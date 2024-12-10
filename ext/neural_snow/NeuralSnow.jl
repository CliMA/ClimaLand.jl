module NeuralSnow

using Flux
using ClimaLand
using ClimaLand.Snow: AbstractDensityModel, SnowModel, SnowParameters
import ClimaLand.Snow:
    density_prog_vars,
    density_prog_types,
    density_prog_names,
    update_density!,
    update_density_prog!,
    snow_depth
import ClimaLand.Parameters as LP
using Thermodynamics

using HTTP, Flux, BSON
include("./ModelTools.jl")
using .ModelTools

export NeuralDepthModel

"""
    NeuralDepthModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density
is calculated from a neural network determining the rate of change of snow depth, dzdt.
"""
struct NeuralDepthModel{FT} <: AbstractDensityModel{FT}
    z_model::Flux.Chain
    α::FT
end

"""
    get_znetwork()
Return the snow-depth neural network from the paper.
"""
function get_znetwork()
    download_link = ClimaLand.Artifacts.neural_snow_znetwork_link()
    z_idx = 1
    p_idx = 7
    nfeatures = 7
    zmodel = ModelTools.make_model(nfeatures, 4, z_idx, p_idx)
    zmodel_state = BSON.load(IOBuffer(HTTP.get(download_link).body))[:zstate]
    Flux.loadmodel!(zmodel, zmodel_state)
    ModelTools.settimescale!(zmodel, 86400.0)
    return zmodel
end

"""
    converted_model_type!(model::Flux.Chain, FT::Type)
A wrapper function to aid conversion of the utilized Flux network weights to the type determined by the simulation.
Currently supports `Float32` and `Float64` types.
"""
function converted_model_type(model::Flux.Chain, FT::DataType)
    if FT == Float32
        return Flux.f32(model)
    elseif FT == Float64
        return Flux.f64(model)
    else
        error(
            "Conversion of Flux Chain weights to the desired float-type is not supported. Please implement the desired conversion
      in NeuralSnow.convert_model_type!() or a custom constructor for the NeuralDepthModel with the desired type.",
        )
    end
end

"""
    NeuralDepthModel(FT::DataType; Δt::Union{Nothing, FT}=nothing; model::Flux.Chain, α::Union{FT, Nothing})
An outer constructor for the `NeuralDepthModel` density parameterization for usage in a snow model.
Can choose to use custom models via the `model` argument, or a custom exponential moving average weight for the network inputs.
The optional Δt argument (model time step, in seconds) can be used to set the lower boundary scaling on the network (when not using a custom model)
in accordance with the paper upon initialization, equivalent to ModelTools.settimescale!(model, Δt, dtype = FT).
If a weight for the exponential moving average is not provided, the default is determined as 1/43200 s⁻¹. This
weight will be unstable for any Δt greater than 43200 seconds, and an alternative value must be supplied.
"""
function NeuralDepthModel(
    FT::DataType;
    Δt::Union{AbstractFloat, Nothing} = nothing,
    model::Flux.Chain = get_znetwork(),
    α::Union{AbstractFloat, Nothing} = nothing,
)
    usemodel = converted_model_type(model, FT)
    weight = !isnothing(α) ? FT(α) : FT(2 / 86400)
    if !isnothing(Δt)
        ModelTools.settimescale!(usemodel, Δt, dtype = FT) #assumes Δt is provided in seconds
        if (Δt > 43200) & (isnothing(α))
            error("Please supply a weight for the
                   exponential moving average, the
                   default value is unstable in this case.")
        end
    end
    return NeuralDepthModel{FT}(usemodel, weight)
end

#Define the additional prognostic variables needed for using this parameterization:
ClimaLand.Snow.density_prog_vars(m::NeuralDepthModel) =
    (:Z, :P_avg, :T_avg, :R_avg, :Qrel_avg, :u_avg)
ClimaLand.Snow.density_prog_types(m::NeuralDepthModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT)
ClimaLand.Snow.density_prog_names(m::NeuralDepthModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)

#Extend/Define the appropriate functions needed for this parameterization:

"""
    snow_depth(m::NeuralDepthModel, Y, p, params)

An extension of the `snow_depth` function to the NeuralDepthModel density parameterization, which includes the prognostic 
depth variable and thus does not need to derive snow depth from SWE and density.
This is sufficient to enable dynamics of the auxillary variable `ρ_snow` without extension of update_density!, and avoids
redundant computations in the computation of runoff.
"""
ClimaLand.Snow.snow_depth(m::NeuralDepthModel, Y, p, params) = Y.snow.Z


"""
    eval_nn(vmodel, z::FT, swe::FT, P::FT, T::FT, R::FT, qrel::FT, u::FT)::FT where {FT}
Helper function for evaluating the neural network in a pointwise manner over a `ClimaCore.Field`
and returning the output in a broadcastable way.
"""
function eval_nn(
    model,
    z::FT,
    swe::FT,
    P::FT,
    T::FT,
    R::FT,
    qrel::FT,
    u::FT,
)::FT where {FT}
    #model() of a Vector{FT} returns a 1-element Matrix, return the internal value: 
    return model([z, swe, qrel, R, u, T, P])[1]
end

"""
    dzdt(density::NeuralDepthModel, model::SnowModel{FT}, Y, p, t) where {FT}
Returns the change in snow depth (rate) given the current model state and the `NeuralDepthModel`
density paramterization, passing the approximate average of the forcings over the last 24 hours instead of
the instantaneous value.
"""
function dzdt(density::NeuralDepthModel, Y)
    return eval_nn.(
        Ref(density.z_model),
        Y.snow.Z,
        Y.snow.S, # When snow-cover-fraction variable is implemented, make sure this value changes to the right input
        Y.snow.P_avg,
        Y.snow.T_avg,
        Y.snow.R_avg,
        Y.snow.Qrel_avg,
        Y.snow.u_avg,
    )
end

"""
    clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT)::FT
A helper function which clips the tendency of Z such that
its behavior is consistent with that of S: if all snow melts
within a timestep, we clip the tendency of S so that it does
not become negative, and here we also clip the tendency of Z
so that depth does not become negative. Additionally, if the
tendencies of Z and S are such that we would encounter Z < S
(rho_snow > rho_liq), we also clip the tendency.
"""
function clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT)::FT where {FT}
    #Case if S is set to zero:
    if (S + dSdt * Δt) <= 0
        return -Z / Δt
        #Case if Z would have been set to Z < S:
    elseif (Z + dZdt * Δt) < (S + dSdt * Δt)
        #more stable form for very small Z, S
        return ((dSdt * Δt + S) - Z) / Δt
    else
        return dZdt
    end
end

"""
    update_density_prog!(density::NeuralDepthModel, model::SnowModel, Y, p)
Updates all prognostic variables associated with density/depth given the current model state and the `NeuralDepthModel`
density paramterization.
"""
function update_density_prog!(
    density::NeuralDepthModel,
    model::SnowModel,
    dY,
    Y,
    p,
)

    dY.snow.Z .=
        clip_dZdt.(
            Y.snow.S,
            Y.snow.Z,
            dY.snow.S, #assumes dY.snow.S is updated (and clipped) before dY.snow.Z
            dzdt(density, Y),
            model.parameters.Δt,
        )

    @. dY.snow.P_avg = density.α * (abs(p.drivers.P_snow) - Y.snow.P_avg)
    @. dY.snow.T_avg =
        density.α *
        (p.drivers.T - model.parameters.earth_param_set.T_freeze - Y.snow.T_avg)
    @. dY.snow.R_avg = density.α * (p.drivers.SW_d - Y.snow.R_avg)
    dY.snow.Qrel_avg .=
        density.α .* (
            Thermodynamics.relative_humidity.(
                model.boundary_conditions.atmos.thermo_params,
                p.drivers.thermal_state,
            ) .- Y.snow.Qrel_avg
        )
    @. dY.snow.u_avg = density.α * (p.drivers.u - Y.snow.u_avg)
end

end

#Todo: Write tests for these functions and a few in Snow.jl or snow_parameterizations.jl

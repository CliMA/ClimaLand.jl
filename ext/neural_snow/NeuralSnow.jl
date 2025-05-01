module NeuralSnow
using StaticArrays
using Flux
using ClimaLand
using ClimaLand.Snow:
    AbstractDensityModel, SnowModel, SnowParameters, snow_bulk_density
import ClimaComms
import ClimaLand.Snow:
    density_prog_vars,
    density_prog_types,
    density_prog_names,
    update_density_and_depth!,
    update_density_prog!
import ClimaLand.Parameters as LP
import Adapt
import ClimaCore
using Thermodynamics


using HTTP, Flux, BSON
include("./ModelTools.jl")
using .ModelTools

export NeuralDepthModel, ModelTools

"""
    NeuralDepthModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density is calculated from a prognostic snow depth variable,
along with the prognostic SWE variable, using a neural network for the rate of change of snow depth, `dz/dt`.
The input to the network are temporally averaged, which we achieve using an exponentially moving average,
with a rate of `\alpha`.
"""
struct NeuralDepthModel{FT, MD <: Flux.Chain} <: AbstractDensityModel{FT}
    "The Flux neural network for compute dz/dt"
    z_model::MD
    "The inverse of the averaging window time (1/s)"
    α::FT
end

function Adapt.adapt_structure(to, from::NeuralDepthModel)
    # Adapt the neural network to the new type
    adapted_model = NeuralDepthModel(
        Adapt.adapt_structure(SArray, from.z_model),
        Adapt.adapt_structure(to, from.α),
    )
    return adapted_model
end
# Adapt.@adapt_structure NeuralDepthModel


# function Adapt.adapt_structure(to, from::Flux.Chain)
#     return Adapt.adapt_structure(to, from)
# end

"""
    get_znetwork()
Return the snow-depth neural network from Charbonneau et al (2025; https://arxiv.org/abs/2412.06819),
and set the network's scaling such that snowpack depth remains nonnegative for a default timestep of 86400 seconds (1 day).
Note that because we are loading a pre-trained model, the number of features input, the size of the output, etc, are hard coded in the function below.
"""
function get_znetwork()
    download_link = ClimaLand.Artifacts.neural_snow_znetwork_link()
    z_idx = 1
    p_idx = 7
    nfeatures = 7
    zmodel = ModelTools.make_model(Float32, nfeatures, 4, z_idx, p_idx)
    zmodel_state = BSON.load(IOBuffer(HTTP.get(download_link).body))[:zstate]
    # Flux.loadmodel!(zmodel, zmodel_state) #Return to this to deterine how to load model for Flux v0.15 and above
    # ModelTools.settimescale!(zmodel, 86400.0)
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
An outer constructor for the NeuralDepthModel density parameterization for usage in a snow model.

By default, the neural network weights trained in Charbonneau et. al. (2025) are used,
but one can also supply their own model by keyword argument. Since the default network
receives inputs representing the average value per feature per day (86400 s), an exponential-weighting
scheme is used to track the moving average of the input fields. A default weighting is supplied
as `\alpha = 2/86400 s⁻¹` following common practice of choosing a weighting of 2/(desired temporal window size),
though the user may specify their own value by keyword argument. Note that this default value is
unstable (whenever `1/\alpha > model timestep`) for model timesteps larger than half a day (43200 seconds).
In the case of using the default model in the Charbonneau et. al. paper, The user should also specify a model timestep
which sets the model scaling in order to guarantee the snow depth remains nonnegative."
"""
function NeuralDepthModel(
    ::Type{FT};
    Δt::Union{AbstractFloat, Nothing} = nothing,
    model::Union{Flux.Chain, Nothing} = nothing,
    α::Union{AbstractFloat, Nothing} = nothing,
) where {FT}
    usemodel = isnothing(model) ? get_znetwork() : model
    usemodel = converted_model_type(usemodel, FT)
    # usemodel = ClimaComms.device() isa ClimaComms.AbstractCPUDevice ? usemodel : Flux.gpu(usemodel)
    weight = !isnothing(α) ? FT(α) : FT(2 / 86400)
    if !isnothing(Δt) & !isnothing(model)
        ModelTools.settimescale!(usemodel, Δt, dtype = FT) #assumes Δt is provided in seconds
        if (Δt > 43200) & (isnothing(α))
            error("Please supply a weight for the
                   exponential moving average, the
                   default value is unstable in this case.")
        end
    end
    return Adapt.adapt(
        ClimaComms.array_type(ClimaComms.device()),
        NeuralDepthModel{FT, typeof(usemodel)}(usemodel, weight),
    )
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
    update_density_and_depth!(ρ_snow, z_snow,density::NeuralDepthModel, Y, p, params::SnowParameters,)

Updates the snow density and depth in place given the current model state. Default for all model types,
can be extended for alternative density parameterizations.
"""
function update_density_and_depth!(
    ρ_snow,
    z_snow,
    density::NeuralDepthModel,
    Y,
    p,
    params::SnowParameters,
)
    @. z_snow = Y.snow.Z
    @. ρ_snow = snow_bulk_density(Y.snow.S, z_snow, params)
end

# TODO: Make this into two flux layers?
"""
    eval_nn(density::NeuralDepthModel, z::FT, swe::FT, P::FT, T::FT, R::FT, qrel::FT, u::FT)::FT where {FT}
Helper function for evaluating the neural network in a pointwise manner over a `ClimaCore.Field`
and returning the output in a broadcastable way.
"""
function eval_nn(
    density::NeuralDepthModel{},
    z::FT,
    swe::FT,
    P::FT,
    T::FT,
    R::FT,
    qrel::FT,
    u::FT,
)::FT where {FT}
    return density.z_model(SVector(z, swe, qrel, R, u, T, P))[]
end

"""
    update_dzdt!(density::NeuralDepthModel, model::SnowModel, Y, p, t)

Updates the dY.snow.Z field in places with the predicted change in snow depth (rate) given the model state `Y` and the `NeuralDepthModel`
density paramterization.
"""
function update_dzdt!(dzdt, density::NeuralDepthModel, Y)
    dzdt .=
        eval_nn.(
            Ref(density),
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
    update_dzdt!(dY.snow.Z, density, Y)

    # Now we clip the tendency so that Z stays within approximately physical bounds.
    @. dY.snow.Z = clip_dZdt(
        Y.snow.S,
        Y.snow.Z,
        dY.snow.S, #assumes dY.snow.S is updated (and clipped) before dY.snow.Z
        dY.snow.Z,
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

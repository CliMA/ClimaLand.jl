module NeuralSnow

using Flux
using ClimaLand
using ClimaLand.Snow: AbstractDensityModel, SnowModel, SnowParameters
import ClimaLand.Snow:
    density_prog_vars,
    density_prog_types,
    density_prog_names,
    update_density!,
    update_density_prog!
import ClimaLand.Parameters as LP
using Thermodynamics

using DataFrames, Dates, CSV, HTTP, Flux, cuDNN, StatsBase, BSON
ModelTools = Base.get_extension(ClimaLand, :NeuralSnowExt).ModelTools
#^^is this correct format since this is within the extension itself, or do i just use include(./ModelTools); using ModelTools; here?

export NeuralDepthModel, snow_depth, update_density_prog!

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
    get_znetwork(; download_link = "https://caltech.box.com/shared/static/ay7cv0rhuiytrqbongpeq2y7m3cimhm4.bson")
A default function for returning the snow-depth neural network from the paper.
"""
#This is the only place that creates the necessity of BSON for the neural extension and the make_model() from ModelTools. Keep or get rid of?
#Can I change this to something CliMA already includes? is BSON vs JLD2 or anything else better? An artifact (can they store BSON/.jld2 files for direct usage without those packages?)?
function get_znetwork(;
    download_link = "https://caltech.box.com/shared/static/ay7cv0rhuiytrqbongpeq2y7m3cimhm4.bson",
)
    z_idx = 1
    p_idx = 7
    nfeatures = 7
    zmodel = ModelTools.make_model(nfeatures, 4, z_idx, p_idx)
    zmodel_state = BSON.load(IOBuffer(HTTP.get(download_link).body))[:zstate]
    Flux.loadmodel!(zmodel, zmodel_state)
    ModelTools.settimescale!(zmodel, 86400.0)
    return zmodel

    #should I write this function to give the ability to also load from local file too?
    #BSON.@load "../SnowResearch/cleancode/testhrmodel.bson" zstate
    #Flux.loadmodel!(zmodel, zstate)
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
If a weight for the exponential moving average is not provided, the weight is determined as the minimum of 2Δt/86400 and 1,
if Δt is also not provided, a default value of 1/10 is used.
"""
function NeuralDepthModel(
    FT::DataType;
    Δt = nothing,
    model::Flux.Chain = get_znetwork(),
    α = nothing,
)
    usemodel = converted_model_type(model, FT)
    weight =
        !isnothing(α) ? FT(α) :
        (isnothing(Δt) ? FT(1 / 10) : FT(min(2 * Δt / 86400, 1)))
    if !isnothing(Δt)
        ModelTools.settimescale!(usemodel, Δt, dtype = FT) #assumes Δt is provided in seconds
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
        Y.snow.S,
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
its behavior is consistent with that of S.
"""
function clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT)::FT where {FT}
    if (S + dSdt * Δt) <= eps(FT)
        return -Z / Δt
    elseif (Z + dZdt * Δt) < (S + dSdt * Δt)
        return (S - Z) / Δt + dSdt
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

    @. dY.snow.P_avg =
        density.α / model.parameters.Δt * (abs(p.drivers.P_snow) - Y.snow.P_avg)
    @. dY.snow.T_avg =
        density.α / model.parameters.Δt *
        (p.drivers.T - model.parameters.earth_param_set.T_freeze - Y.snow.T_avg)
    @. dY.snow.R_avg =
        density.α / model.parameters.Δt * (p.drivers.SW_d - Y.snow.R_avg)
    dY.snow.Qrel_avg .=
        density.α / model.parameters.Δt .* (
            Thermodynamics.relative_humidity.(
                model.boundary_conditions.atmos.thermo_params,
                p.drivers.thermal_state,
            ) .- Y.snow.Qrel_avg
        )
    @. dY.snow.u_avg =
        density.α / model.parameters.Δt * (p.drivers.u - Y.snow.u_avg)
end

end

#Todo: Write tests for these functions
#Todo: address Anderson1976 sensitivity for small precision snowpacks (z = 0 when SWE > 0 or z > 0 when SWE = 0, but on the order of 1e-10 < eps(Float32), doesn't make NaN at present though)

# The current setup requires the following flowchart, if one wants to add a new density parameterization:
# If the paramterization does not require any additional prognostic variables, you only must redefine/extend update_density!() (this is regardless of requiring 0 or multiple new auxillary variables)
# If the parameterization requires additional prognostic variables, you must redefine/extend update_density_prog!()
#   If these prognostic variables do not include depth, you must extend update_density!() in addition to update_density_prog!()
#   If these prognostic variables include depth and no new auxillary variables are needed, you can choose to only redefine snow_depth() (instead of update_density!()) in addition to update_density_prog!().
#   If these prognostic variables include depth and new auxillary variables are needed, must extend update_density_prog!() and update_density!, and extending snow_depth() is recommended but not necessary.
# All these assume someone would not make a parameterization that makes ρ_snow a prognostic variable (and thus have to remove it from the auxillary variables), aka only new diagnostic parameterizations for density are made.

#Is there away to further simplify/collapse this structure so new density parameterizations only require extension of the same functions each time?

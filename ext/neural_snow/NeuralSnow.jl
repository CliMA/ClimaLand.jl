module NeuralSnow

using Flux
using ClimaLand
using ClimaLand.Snow: AbstractDensityModel, SnowModel, SnowParameters
import ClimaLand.Snow:
    density_prog_vars,
    density_prog_types,
    density_prog_names,
    density_aux_vars,
    density_aux_types,
    density_aux_names,
    update_density!,
    update_density_prog!
import ClimaLand.Parameters as LP
using Thermodynamics

using DataFrames, Dates, CSV, HTTP, Flux, cuDNN, StatsBase, BSON
ModelTools = Base.get_extension(ClimaLand, :NeuralSnowExt).ModelTools
#^^is this correct format since this is within the extension itself, or do i just use include(./ModelTools); using ModelTools; here?
#the only reason I am using this is for the make_model() (once) and settimescale!() (twice) functions.

export NeuralDepthModel, update_density!, update_density_prog!

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
The Δt argument can be used to set the lower boundary value on the network (when not using a custom model) in accordance with the paper upon initialization,
equivalent to ModelTools.settimescale!(model, Δt, dtype = FT).
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
        (isnothing(Δt) ? FT(1 / 12) : FT(maximum([2 * Δt / 86400, 1])))
    # should we set the timescale of the model to Δt? getting rid of get_znetwork() and below line mitigates the need to import ModelTools entirely.
    if !isnothing(Δt)
        ModelTools.settimescale!(usemodel, Δt, dtype = FT) #assumes Δt is provided in seconds
    end
    return NeuralDepthModel{FT}(usemodel, weight)
end

#Define the additional prognostic variables needed for using this parameterization:
#do I need to export these for them to be utilized in the simulation?
ClimaLand.Snow.density_prog_vars(m::NeuralDepthModel) =
    (:Z, :P_avg, :T_avg, :R_avg, :Qrel_avg, :u_avg)
ClimaLand.Snow.density_prog_types(m::NeuralDepthModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT)
ClimaLand.Snow.density_prog_names(m::NeuralDepthModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)

#so far, every type is just blank here - is it even worth providing this feature or just get rid of this feature in Snow.jl entirely?
ClimaLand.Snow.density_aux_vars(m::NeuralDepthModel) = ()
ClimaLand.Snow.density_aux_types(m::NeuralDepthModel{FT}) where {FT} = ()
ClimaLand.Snow.density_aux_names(m::NeuralDepthModel) = ()

"""
    snow_density(density::NeuralDepthModel, SWE::FT, z::FT, parameters::SnowParameters{FT}) where {FT}
Returns the snow density given the current model state and the `NeuralDepthModel`
density paramterization.
"""
function snow_density(
    density::NeuralDepthModel{FT},
    SWE::FT,
    z::FT,
    parameters::SnowParameters{FT},
)::FT where {FT}
    ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    if SWE <= eps(FT) #if there is no snowpack, aka avoid NaN
        return FT(ρ_l)
    end
    ρ_new = SWE / z * ρ_l
    return ρ_new
end

"""
    update_density!(density::NeuralDepthModel, params::SnowParameters, Y, p)
Updates the snow density given the current model state and the `NeuralDepthModel`
density paramterization.
"""
function update_density!(
    density::NeuralDepthModel,
    params::SnowParameters,
    Y,
    p,
)
    p.snow.ρ_snow .= snow_density.(Ref(density), Y.snow.S, Y.snow.Z, params)
end

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
    input = FT.([z, swe, qrel, R, u, T, P])
    return model(input)[1]
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
        #do we need to do something about setting q_l = 1 here?
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
    #Get Inputs:
    #is this too many allocations?
    z = Y.snow.Z
    swe = Y.snow.S
    dswedt = p.snow.applied_water_flux #do i need to multiply by -1 here for directions?
    dprecipdt_snow = abs.(p.drivers.P_snow) #need to change downward direction to scalar
    air_temp = p.drivers.T .- model.parameters.earth_param_set.T_freeze
    sol_rad = abs.(p.drivers.SW_d) #need to change downward direction to scalar
    rel_hum =
        Thermodynamics.relative_humidity.(
            Ref(model.boundary_conditions.atmos.thermo_params),
            p.drivers.thermal_state,
        )
    wind_speed = p.drivers.u
    dt = model.parameters.Δt
    β = density.α

    dY.snow.Z .= clip_dZdt.(swe, z, dswedt, dzdt(density, Y), dt)
    dY.snow.P_avg = β .* (dprecipdt_snow .- Y.snow.P_avg)
    dY.snow.T_avg = β .* (air_temp .- Y.snow.T_avg)
    dY.snow.R_avg = β .* (sol_rad .- Y.snow.R_avg)
    dY.snow.Qrel_avg = β .* (rel_hum .- Y.snow.Qrel_avg)
    dY.snow.u_avg = β .* (wind_speed .- Y.snow.u_avg)
end

end

#do tests need to be made for these functions? since it is an extension?

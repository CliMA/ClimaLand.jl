module NeuralSnow
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
import ClimaCore
using Thermodynamics

using Flux
using Adapt, StaticArrays #only needed to circumvent make_static_model() warning
using ..ConstrainedNeuralModels

export NeuralDepthModel

"""
    Snow_Depth_Lower_Bound

Establishes a functor type to use for the lower bound within the 
ConstrainedNeuralModel that will serve as the NeuralDepthModel.
Tracks the index of the input feature associated with the snow depth.
"""
@bound_type struct Snow_Depth_Lower_Bound
    z_idx::Int
end

"""
    (b::Snow_Depth_Lower_Bound)(
        pred::SVector{1, FT},
        input::SVector{7, FT},
    )::FT where {FT <: AbstractFloat}

Establishes a functor method to use as the static lower bound for single-inputs  
within the NeuralDepthModel. Since the NeuralDepthModel will use custom
fixed layers for increased performance, the actual boundary function
`-z/Δt` only needs to be passed here as `z`, meaning the element of
the input vector associated with the snow depth, i.e. the function
is merely `input[b.z_idx]`.
"""
@bound function (b::Snow_Depth_Lower_Bound)(
    pred::SVector{1, FT},
    input::SVector{7, FT},
)::FT where {FT <: AbstractFloat}
    return input[b.z_idx]
end

"""
    Snow_Depth_Upper_Bound

Establishes a functor type to use for the upper bound within the 
ConstrainedNeuralModel that will serve as the NeuralDepthModel.
Tracks the index of the input feature associated with the precipitation.
"""
@bound_type struct Snow_Depth_Upper_Bound
    precip_idx::Int
end

"""
    (b::Snow_Depth_Upper_Bound)(
        pred::SVector{1, FT},
        input::SVector{7, FT},
    )::FT where {FT <: AbstractFloat}

Establishes a functor method to use as the static upper bound for single-inputs  
within the NeuralDepthModel. No upper bound should exist whenever snow
precipitation exists, and otherwise the upper limit is dz/dt = 0. This is
represented by the formula: `(input[b.precip_idx] > 0) * relu(pred[1])`.
"""
@bound function (b::Snow_Depth_Upper_Bound{FT})(
    pred::SVector{1, FT},
    input::SVector{7, FT},
)::FT where {FT <: AbstractFloat}
    return >(input[b.precip_idx], 0) * relu(pred[1])
end

"""
    _set_timescale!(m::ConstrainedNeuralModel, Δt::Real)

A utility method specifically for the ConstrainedNeuralModel used
for the NeuralDepthModel. This should not be used for arbitrary
ConstrainedNeuralModels. This serves to apply the `1/Δt` scaling
to the lower boundary function, allowing the timescale to be set
prior to fixing the model, and takes this specific form from the
custom fixed layers used for this model. It can only be called
prior to the model being made static.
"""
function _set_timescale!(m::ConstrainedNeuralModel, Δt::Real)
    m.fixed_layers[1].weight[2, 2] = 1 * eltype(m.out_scale)(1 / Δt)
    m.initial_fixed_layer[2, 2] = 1 * eltype(m.out_scale)(1 / Δt)
end

#=
#=
This below code shows what building the same model from scratch
would look like:
(Undecided on whether to leave this or delete it by the merge
of this branch, as the same information will be visible in
the model metadata file.)
=#
const depth_model_in_scales::AbstractArray{Float32} =  1 ./ [
            0.68659294 # z (m)
            0.25578135 # SWE (m)
            0.20772743 # relative humidity (0-1)
            76.2825 # solar radiation (W/m²)
            0.63011056 # wind speed (m/s)
            6.3486657 # air temp (⁰C)
            6.9992964f-8 # water-equiv rate of snowfall (m/s)
        ]

const n_features::Int = 7
const k::Int = 4
const z_idx::Int = 1
const p_idx::Int = 7
const Δt::Int = 86400
const snow_depth_out_scale::Float32 = 7.643518e-6

function get_snow_depth_fixed_layers(::Type{FT}) where {FT<:AbstractFloat}
    get_relus = Matrix{FT}([1 0 0; 0 1 0; 1 0 -1])
    get_min = Matrix{FT}([1 1 -1; 0 1 0])
    get_max = Matrix{FT}([1 -1])
    return Chain(
        Dense(get_relus, false, relu),
        Dense(get_min, false, relu),
        CNM.MulLayer(get_max)
    )
end

snow_depth_predmodel = Chain(
    Dense(n_features, nfeatures*k, relu),
    Dense(n_features*k, n_features, elu),
    Dense(n_features, 1),
)

const m = ConstrainedNeuralModel(
    Float32,
    snow_depth_predmodel,
    upper_bound = Snow_Depth_Upper_Bound(p_idx),
    lower_bound = Snow_Depth_Lower_Bound(z_idx),
    in_scales = depth_model_in_scales,
    out_scale = snow_depth_out_scale,
    fixed_layers = get_snow_depth_fixed_layers(Float32, Δt)
)
=#

"""
    NeuralDepthModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density is calculated from a prognostic snow depth variable,
along with the prognostic SWE variable, using a ConstrainedNeuralModel for the rate of change of snow depth, `dz/dt`.
The input to the network are temporally averaged, which we achieve using an exponentially moving average,
with a rate of `\alpha`.

During construction, all Arrays in the provided model are converted to StaticArrays (this is done
for GPU compatibilty). For GPU compatibilty, the model must not allocate or take the adjoint
of arrays.
"""
struct NeuralDepthModel{FT, MD <: ConstrainedNeuralModel} <:
       AbstractDensityModel{FT}
    "The Flux neural network for compute dz/dt"
    z_model::MD
    "The inverse of the averaging window time (1/s)"
    α::FT
    function NeuralDepthModel(z_model::MD, α::FT) where {FT, MD}
        #=
        Usually, we'd just use ConstrainedNeuralModel.make_static_model(z_model)
        here, but as we specifically only want the single-input version
        for modeling here, we circumvent the method to avoid the warning
        it would output about not also enabling batched inputs. On the
        other hand, this is the only reason Adapt/StaticArray need to
        be additionally imported into this module:
        =#
        static_z_model = Adapt.adapt(SArray, z_model)
        new{FT, typeof(static_z_model)}(static_z_model, α)
    end
end

"""
    get_znetwork(FT)
Return the snow-depth neural network from Charbonneau et al (2025; https://doi.org/10.1175/AIES-D-24-0040.1),
with the network's scaling such that snowpack depth remains nonnegative for a default timestep of 86400 seconds (1 day).
"""
function get_znetwork(FT)
    #TODO: FIX THIS METHOD ONCE THE RIGHT ARTIFACTS ARE IN
    #=
    download_link = ClimaLand.Artifacts.neural_snow_znetwork_link()
    z_idx = 1
    p_idx = 7
    nfeatures = 7
    zmodel = ModelTools.make_model(FT, nfeatures, 4, z_idx, p_idx)
    zmodel_state = BSON.load(IOBuffer(HTTP.get(download_link).body))[:zstate]
    Flux.loadmodel!(zmodel, zmodel_state) #Return to this to deterine how to load model for Flux v0.15 and above
    ModelTools.settimescale!(zmodel, 86400.0)
    return zmodel =#
end


"""
    NeuralDepthModel(FT::DataType; Δt::Union{Nothing, FT}=nothing; model::Flux.Chain, α::Union{FT, Nothing})
An outer constructor for the NeuralDepthModel density parameterization for usage in a snow model.

This model utilizes the neural network formulated in Charbonneau et. al. (2025) to predict the rate of change of snow height dz/dt.
Since the default network receives inputs representing the average value per feature per day (86400 s), an exponential-weighting
scheme is used to track the moving average of the input fields. A default weighting is supplied
as `\alpha = 2/86400 s⁻¹` following common practice of choosing a weighting of 2/(desired temporal window size),
though the user may specify their own value by keyword argument. Note that this default value is
unstable (whenever `1/\alpha > model timestep`) for model timesteps larger than half a day (43200 seconds).
The user should also specify a model timestep Δt in seconds which sets the model scaling in order to guarantee the snow depth remains nonnegative.
"""
function NeuralDepthModel(
    ::Type{FT};
    Δt::Union{<:Real, Nothing} = nothing,
    α::Union{<:Real, Nothing} = nothing,
) where {FT <: AbstractFloat}
    usemodel = get_znetwork(FT)
    weight = !isnothing(α) ? FT(α) : FT(2 / 86400)
    if !isnothing(Δt)
        _set_timescale!(usemodel, Δt) #assumes Δt is provided in seconds
        if (Δt > 43200) && isnothing(α)
            error("Please supply a weight for the
                   exponential moving average, the
                   default value is unstable in this case.")
        end
    end
    return NeuralDepthModel(usemodel, weight)
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


"""
    eval_nn(density::NeuralDepthModel, z::FT, swe::FT, P::FT, T::FT, R::FT, qrel::FT, u::FT)::FT where {FT}

Helper function for evaluating the neural network in a pointwise manner over a `ClimaCore.Field`
and returning the output in a broadcastable way.
"""
function eval_nn(
    density::NeuralDepthModel,
    z::FT,
    swe::FT,
    P::FT,
    T::FT,
    R::FT,
    qrel::FT,
    u::FT,
)::FT where {FT}
    # model() of a Vector{FT} returns a 1-element Matrix
    # use SVector for gpu compatibility

    #use "SWE" or use "scf-adjusted SWE"?
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
            Y.snow.S,
            Y.snow.P_avg,
            Y.snow.T_avg,
            Y.snow.R_avg,
            Y.snow.Qrel_avg,
            Y.snow.u_avg,
            #might need to include the snow cover fraction here from p.snow.scf to get the right value
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
#Consider how snow cover fraction will impact this:
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
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
            ) .- Y.snow.Qrel_avg
        )
    @. dY.snow.u_avg = density.α * (p.drivers.u - Y.snow.u_avg)
end

end

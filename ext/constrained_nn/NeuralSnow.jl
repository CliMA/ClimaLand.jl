module NeuralSnow
using ClimaLand
using ClimaLand.Snow:
    AbstractDensityModel, AbstractAlbedoModel, SnowModel, snow_bulk_density
import ClimaComms
import ClimaLand.Snow:
    extra_prog_vars,
    extra_prog_types,
    extra_prog_domain_names,
    update_density_and_depth!,
    update_snow_albedo!,
    compute_extra_prog_tendency!
import ClimaLand.Parameters as LP
import ClimaCore
import ClimaParams as CP
using Thermodynamics

using Flux: relu
using StaticArrays
using ..ConstrainedNeuralModels

export NeuralDepthModel, NeuralAlbedoModel

"""
   AbstractResetAlbedoModel{FT<:AbstractFloat}

Defines the parameterization type for getting the initial value of a new snowpack, to
be combined with the prognostic update dynamics of a NeuralAlbedoModel.

An instance of this type must have at least:
    1) one method defined, `get_initial_albedo(alb_reset_model::T, Y, p)`
       for T <: AbstractResetAlbedoModel that returns the initial albedo value,
       which will be called in `albedo_reset_rate()`
    2) a float field "reset_limit" (m/s) of the critical liquid-equivalent snowfall preciptiation rate,
       which will also get used in resetting the snow albedo.
    3) a float field "no_snowpack_S" (m) of a delimiting height for which a snowpack is considered new if
       S < no_snowpack_S and P > reset_limit
"""
abstract type AbstractResetAlbedoModel{FT <: AbstractFloat} end

"""
    ConstantResetAlbedo{FT} <: AbstractResetAlbedoModel{FT}
Establishes the parameterization for initial values of albedo in a new snowpack,
where the dynamics of the model are begun from a specified constant value.
Default values are a new albedo of 0.8 for 3 mm liquid-equivalent of snow precipitation in 24 hours.
"""
struct ConstantResetAlbedo{FT} <: AbstractResetAlbedoModel{FT}
    "Initial albedo constant (unitless 0-1)"
    init_α::FT
    "Threshold of liquid-equiavlent snow precipitation for a reset (m/s)"
    reset_limit::FT
    "Small/no snowpack limit which S must be under to trigger a reset with suffifient precipitation"
    no_snowpack_S::FT
end
function ConstantResetAlbedo(
    toml_dict::CP.ParamDict{FT};
    init_alpha = toml_dict["snow_initial_albedo"],
    reset_limit = toml_dict["albedo_psnow_reset_rate"],
    no_snowpack_S = toml_dict["new_snowpack_S_threshold"],
) where {FT <: AbstractFloat}
    return ConstantResetAlbedo{FT}(init_alpha, reset_limit, no_snowpack_S)
end

"""
    get_initial_albedo(alb_reset_model::ConstantResetAlbedo, Y, p)

Defines the behavior for a ConstantResetAlbedo parameterization,
which is called within `albedo_reset_rate()` and combined with the
NeuralAlbedoModel output in `update_dαdt!`.
"""
function get_initial_albedo(alb_reset_model::ConstantResetAlbedo, Y, p)
    return alb_reset_model.init_α
end

"""
    Snow_Depth_Lower_Bound

Establishes a functor type to use for the lower bound within the 
ConstrainedNeuralModel that will serve as the NeuralDepthModel.
Tracks the index of the input feature associated with the snow depth.
"""
struct Snow_Depth_Lower_Bound
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
function (b::Snow_Depth_Lower_Bound)(
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
struct Snow_Depth_Upper_Bound
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
function (b::Snow_Depth_Upper_Bound)(
    pred::SVector{1, FT},
    input::SVector{7, FT},
)::FT where {FT <: AbstractFloat}
    return >(input[b.precip_idx], 0) * relu(pred[1])
end

"""
    Snow_Albedo_Bound{FT <: AbstractFloat}

Establishes a functor type to use for the bounds within the 
ConstrainedNeuralModel that will serve as the NeuralAlbedoModel.
Tracks the index of the input feature associated with the albedo,
as well as a limiting value and the model timescale (the inverse
of it, since multiplication is faster than division).
"""
struct Snow_Albedo_Bound{FT <: AbstractFloat, T <: AbstractVector{FT}}
    alb_idx::Int
    lim::T
    t_scale_factor::T
end

"""
    (b::Snow_Albedo_Bound)(
        pred::SVector{1, FT},
        input::SVector{6, FT},
    )::FT where {FT <: AbstractFloat}

Establishes a functor method to use as the static bound for single-inputs  
within the NeuralAlbedoModel. The bound exists to stop the model from going
beyond the limiting value set by the Snow_Albedo_Bound. The formula is:
(b.lim[1] - input[b.alb_idx]) * b.t_scale_factor[1] (b.t_scale_factor = [1/Δt])
"""
function (b::Snow_Albedo_Bound)(
    pred::SVector{1, FT},
    input::SVector{6, FT},
)::FT where {FT <: AbstractFloat}
    return (b.lim[1] - input[b.alb_idx]) * b.t_scale_factor[1]
end

"""
    _set_timescale!(m::ConstrainedNeuralModel, Δt::Real)

A utility method specifically for the ConstrainedNeuralModels used
for the NeuralSnow module. This should not be used for arbitrary
ConstrainedNeuralModels. This serves to apply the `1/Δt` scaling
to the appropriate constraints, allowing the timescale to be set
prior to fixing the model, with the NeuralDepthModel and NeuralAlbedoModel
having slightly different implementations due to their custom structures.
It can only be called prior to the models being made static.
"""
function _set_timescale!(m::ConstrainedNeuralModel, Δt::Real)
    FT = eltype(m.out_scale.sc)
    if typeof(m.constraints.upper_bound) <: Snow_Depth_Upper_Bound
        m.fixed_layers[1].weight[2, 2] = 1 * FT(1 / Δt)
        m.initial_fixed_layer[2, 2] = 1 * eltype(m.out_scale.sc)(1 / Δt)
    elseif typeof(m.constraints.upper_bound) <: Snow_Albedo_Bound
        m.constraints.upper_bound.t_scale_factor[1] = FT(1 / Δt)
        m.constraints.lower_bound.t_scale_factor[1] = FT(1 / Δt)
    else
        error("Unsupported NeuralSnow model type.")
    end
end

"""
    NeuralDepthModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density is calculated from a prognostic snow depth variable,
along with the prognostic SWE variable, using a ConstrainedNeuralModel for the rate of change of snow depth, `dz/dt`.
The input to the network are temporally averaged, which we achieve using an exponentially moving average,
with a rate of `w`.

During construction, all Arrays in the provided model are converted to StaticArrays (this is done
for GPU compatibilty). For GPU compatibilty, the model must not allocate or take the adjoint
of arrays.
"""
struct NeuralDepthModel{FT, MD <: ConstrainedNeuralModel} <:
       AbstractDensityModel{FT}
    "The ConstrainedNeuralModel for computing dz/dt"
    z_model::MD
    "The inverse of the averaging window time (1/s)"
    w::FT
    "The minimum allowable density of the snowpack under this evolution, as a fraction of the liquid water density"
    ρ_min_frac::FT
    function NeuralDepthModel(z_model::MD, w::FT, ρ_min_frac::FT) where {FT, MD}
        static_z_model = make_static_model(z_model, skip_check = true)
        new{FT, typeof(static_z_model)}(static_z_model, w, ρ_min_frac)
    end
end

"""
    NeuralAlbedoModel{FT <: AbstractFloat} <: AbstractAlbedoModel{FT}
Establishes the albedo parameterization where snow albedo is calculated with a ConstrainedNeuralModel
for the rate of change of snow albedo, `dα/dt`. When a new snowpack
appears, the supplied initial albedo model is used to update the value.

During construction, all Arrays in the provided model are converted to StaticArrays (this is done
for GPU compatibilty). For GPU compatibilty, the model must not allocate or take the adjoint
of arrays.
"""
struct NeuralAlbedoModel{
    FT <: AbstractFloat,
    AM <: ConstrainedNeuralModel,
    IA <: AbstractResetAlbedoModel,
    ZAF <: Function,
} <: AbstractAlbedoModel{FT}
    "The ConstrainedNeuralModel for computing dα/dt"
    alb_model::AM
    "The initial albedo parameterization"
    new_alb::IA
    "The function for interacting with the larger model to get zenith angle at solar noon"
    za_solarnoon::ZAF
    function NeuralAlbedoModel(
        alb_model::AM,
        init_model::IA,
        zaf::ZAF,
    ) where {AM, IA, ZAF}
        static_alb_model = make_static_model(alb_model, skip_check = true)
        FT = eltype(static_alb_model.out_scale.sc)
        new{FT, typeof(static_alb_model), typeof(init_model), typeof(zaf)}(
            static_alb_model,
            init_model,
            zaf,
        )
    end
end

"""
    get_znetwork(FT)
Return the snow-depth neural network from Charbonneau et al (2025; https://doi.org/10.1175/AIES-D-24-0040.1),
with the network's scaling such that snowpack depth remains nonnegative for a default timestep of 86400 seconds (1 day).
"""
function get_znetwork(FT; params = nothing)
    structure_path = ClimaLand.Artifacts.neural_depth_model_structure_path()
    if isnothing(params)
        param_path = ClimaLand.Artifacts.neural_depth_model_params_path()
        return convert_model(load_model(param_path, structure_path), FT)
    else
        @assert params isa Vector "Parameters must be supplied as a flattened vector."
        zmodel_data = load_model(structure_path)
        return convert_model(load_model(FT.(params), zmodel_data), FT)
    end
end

"""
    get_albnetwork(FT)
Return the snow-albedo neural network developed by Charbonneau et al (under review),
with the network's scaling such that snowpack albedo remains between 0.3 and 1.

NOTE: the albedo network involves a sinh() activation function, which can make it
more unstable to inputs that are not the correct units. Ensure inputs are consistent
with model design.
"""
function get_albnetwork(FT; params = nothing)
    structure_path = ClimaLand.Artifacts.neural_albedo_model_structure_path()
    if isnothing(params)
        param_path = ClimaLand.Artifacts.neural_albedo_model_params_path()
        return convert_model(load_model(param_path, structure_path), FT)
    else
        @assert params isa Vector "Parameters must be supplied as a flattened vector."
        alb_model_data = load_model(structure_path)
        return convert_model(load_model(FT.(params), alb_model_data), FT)
    end
end

"""
    NeuralDepthModel(
        toml_dict::CP.ParamDict;
        model_params::Union{Vector{<:Real}, Nothing},
        Δt::Union{Nothing, FT}=nothing,
        w::Union{FT, Nothing},
        minimum_density_fraction = toml_dict["minimum_NDM_density"]
    )
An outer constructor for the NeuralDepthModel density parameterization for usage in a snow model.

This model utilizes the neural network formulated in Charbonneau et. al. (2025) to predict the rate of change of snow height dz/dt.
Since the default network receives inputs representing the average value per feature per day (86400 s), an exponential-weighting
scheme is used to track the moving average of the input fields. A default weighting is supplied
as `w = 2/86400 s⁻¹` following common practice of choosing a weighting of 2/(desired temporal window size),
though the user may specify their own value by keyword argument. Note that this default value is
unstable (whenever `1/w > model timestep`) for model timesteps larger than half a day (43200 seconds).
The user should also specify a model timestep Δt in seconds which sets the model scaling in order to guarantee the snow depth remains nonnegative.
With the defaults, this model requires about 1.5 days of burn-in time for the EMA to properly reflect the 24-hour moving average.
For a custom weight `w` and time-step Δt, the burn-in time (in seconds) is roughly 3/w whenever wΔt << 1.
"""
function NeuralDepthModel(
    toml_dict::CP.ParamDict{FT};
    model_params::Union{Vector{<:Real}, Nothing} = nothing,
    Δt::Union{<:Real, Nothing} = nothing,
    w::Union{<:Real, Nothing} = nothing,
    minimum_density_fraction = FT(toml_dict["minimum_NDM_density"]),
) where {FT <: AbstractFloat}
    _DT_ = FT(LP.LandParameters(toml_dict).insol_params.day)
    usemodel = get_znetwork(FT, params = model_params)
    weight = !isnothing(w) ? FT(w) : FT(2 / _DT_)
    if !isnothing(Δt)
        _set_timescale!(usemodel, Δt) #assumes Δt is provided in seconds
        if (Δt > _DT_ / 2) && isnothing(w)
            error("Please supply a weight for the
                   exponential moving average, the
                   default value is unstable in this case.")
        end
    end
    return NeuralDepthModel(usemodel, weight, minimum_density_fraction)
end

"""
    NeuralAlbedoModel(
        toml_dict,
        surface_space;
        initial_albedo = ConstantResetAlbedo(FT),
        model_params::Union{Vector{<:Real}, Nothing} = nothing,
        Δt::Union{<:Real, Nothing} = nothing,
    ) where {FT <: AbstractFloat}

An outer constructor for the NeuralAlbedoModel albedo parameterization for usage in a snow model. It requires model domains that 
contain latitude and longitude information, requiring the passing of the model's surface space to the constructor.

This model utilizes the neural network develped by Charbonneau et. al. (2025) to predict the rate of change of snow albedo dα/dt.
The user should also specify a model timestep Δt in seconds which sets the model scaling in order to guarantee the snow albedo remains physical.
If unspecified, the initial albedo value of a new snowpack is dictated by a ConstantResetAlbedo parameterization of 0.8.
"""
function NeuralAlbedoModel(
    toml_dict::CP.ParamDict{FT},
    surface_space;
    initial_albedo = ConstantResetAlbedo(toml_dict),
    model_params::Union{Vector{<:Real}, Nothing} = nothing,
    Δt::Union{<:Real, Nothing} = nothing,
) where {FT <: AbstractFloat}
    usemodel = get_albnetwork(FT, params = model_params)
    !isnothing(Δt) && _set_timescale!(usemodel, Δt)
    model_lat = ClimaCore.Fields.coordinate_field(surface_space).lat
    model_lon = ClimaCore.Fields.coordinate_field(surface_space).long
    insol_ps = LP.LandParameters(toml_dict).insol_params
    zaf =
        (t, start_date) -> ClimaLand.solar_noon_zenith_cosine(
            t,
            start_date;
            lat = model_lat,
            lon = model_lon,
            insol_params = insol_ps,
        )
    return NeuralAlbedoModel(usemodel, initial_albedo, zaf)
end

#Define the additional prognostic variables needed for using these parameterization:
ClimaLand.Snow.extra_prog_vars(::NeuralDepthModel) =
    (:Z, :P_avg, :T_avg, :R_avg, :Qrel_avg, :u_avg)
ClimaLand.Snow.extra_prog_vars(::NeuralAlbedoModel) = (:A,)

ClimaLand.Snow.extra_prog_types(::NeuralDepthModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT)
ClimaLand.Snow.extra_prog_types(::NeuralAlbedoModel{FT}) where {FT} = (FT,)

ClimaLand.Snow.extra_prog_domain_names(::NeuralDepthModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)
ClimaLand.Snow.extra_prog_domain_names(::NeuralAlbedoModel) = (:surface,)

#Extend/Define the appropriate functions needed for this parameterization:

"""
    update_density_and_depth!(ρ_snow, z_snow,density::NeuralDepthModel, Y, p, earth_param_set)

Updates the snow density and depth in place given the current model state. Dispatched form for the
NeuralDepthModel type.
"""
function ClimaLand.Snow.update_density_and_depth!(
    ρ_snow,
    z_snow,
    density::NeuralDepthModel,
    Y,
    p,
    earth_param_set,
)
    #For now: assume this model was trained to represent Y.snow.Z as z-per-ground-area, just like p.snow.Z is, and take in per-ground-area inputs.
    #Do one extra clamp to deal with machine precision issues:
    @. z_snow = clamp(Y.snow.Z, Y.snow.S, Y.snow.S / density.ρ_min_frac) #p.snow.z_snow is z-per-ground-area, so we set its value directly. Need to change to "Y.snow.Z * scf" if we change Y.snow.Z to represent true snow depth.
    @. ρ_snow = snow_bulk_density(Y.snow.S, z_snow, earth_param_set) #make sure both passed args are both per-snow-area or both per-ground-area.
end

"""
    update_snow_albedo!(α, m::NeuralAlbedoModel, Y, p, t, earth_param_set,)

Updates the snow albedo in place given the current model state. Dispatched form
specifically for the NeuralAlbedoModel type.
"""
function ClimaLand.Snow.update_snow_albedo!(
    α,
    m::NeuralAlbedoModel,
    Y,
    p,
    t,
    earth_param_set,
)
    c = m.alb_model.constraints
    @. α = clamp(Y.snow.A, c.lower_bound.lim[1], c.upper_bound.lim[1]) #try clipping again for machine precision concerns?
end

"""
    eval_z_nn(z_model::ConstrainedNeuralModel, z::FT, swe::FT, P::FT, T::FT, R::FT, qrel::FT, u::FT)::FT where {FT}

Helper function for evaluating the neural network in a pointwise manner over a `ClimaCore.Field`
and returning the output in a broadcastable way.
"""
function eval_z_nn(
    z_model::ConstrainedNeuralModel,
    z::FT,
    swe::FT,
    P::FT,
    T::FT,
    R::FT,
    qrel::FT,
    u::FT,
)::FT where {FT}
    # model() of a Vector{FT} returns a 1-element Matrix
    # use SVector for GPU compatibility

    #See above comment - assume model predicts per-ground-area from per-ground-area inputs, since we'd need to find a `swe_per_snow_area(SWE, scf)` function
    #that is physically consistent in the small snowpack limit AND doesn't either:
    # - make the snowpack start to grow again as scf gets smaller faster than z or SWE gets smaller
    # - make SWE-per-snow-area >> z-per-snow-area within some regime
    return z_model(SVector(z, swe, qrel, R, u, T, P))[]
end


"""
    eval_alb_nn(albmodel::ConstrainedNeuralModel, α::FT, P::FT, T::FT, ZA::FT)::FT where {FT}

Helper function for evaluating the neural network in a pointwise manner over a `ClimaCore.Field`
and returning the output in a broadcastable way.
"""
function eval_alb_nn(
    albmodel::ConstrainedNeuralModel,
    α::FT,
    P::FT,
    ZA::FT,
    T::FT,
)::FT where {FT}
    # model() of a Vector{FT} returns a 1-element Matrix
    # use SVector for GPU compatibility
    return albmodel(SVector(α, P, ZA, T, α * T, α * P))[]
end

"""
    ema_tendency(w::FT, new::FT, old::FT)::FT where {FT}

A utility function to get the derivative value of an exponential moving average (EMA),
in terms of its old value, the incoming value, and the scaling weight.
"""
function ema_tendency(w::FT, new::FT, old::FT)::FT where {FT}
    return w * (new - old)
end

"""
    snow_SW_down(model_type::Val{(:snow,)}, p)

Obtains the correct downward shortwave radiation magnitude for a snow model without an
overlying canopy (standalone or just coupled with soil).
"""
@inline snow_SW_down(model_type::Union{Val{(:snow,)}, Val{(:snow, :soil)}}, p) =
    p.drivers.SW_d

"""
    snow_SW_down(model_type::Union{
        Val{(:canopy, :snow, :soil, :soilco2)},
        Val{(:canopy, :snow, :soil)},
        }, p)

Obtains the correct downward shortwave radiation magnitude for a snow model with
an overlying canopy.
"""
@inline function snow_SW_down(
    model_type::Union{
        Val{(:canopy, :snow, :soil, :soilco2)},
        Val{(:canopy, :snow, :soil)},
    },
    p,
)
    return @. (
        p.canopy.radiative_transfer.nir.trans *
        p.canopy.radiative_transfer.nir_d +
        p.canopy.radiative_transfer.par.trans *
        p.canopy.radiative_transfer.par_d
    ) #network just takes incoming shortwave, not net
end

"""
    reset_alb_tendency(S::FT, P::FT, S_lim::FT, P_lim::FT, reset_rate::FT, prognostic_rate::FT)::FT where {FT}

Handles albedo model prognostic updates, returning the albedo to the initial value
as dicated by an AbstractResetAlbedoModel if a new snowpack is being formed, and
otherwise returning the rate as dicated by the prognostic model.
"""
function reset_alb_tendency(
    S::FT,
    P::FT,
    S_lim::FT,
    P_lim::FT,
    reset_rate::FT,
    prognostic_rate::FT,
)::FT where {FT}
    return (S < S_lim) && (P > P_lim) ? reset_rate : prognostic_rate
end

"""
    albedo_reset_rate(alb_reset_model::AbstractResetAlbedoModel, Δt::FT, Y, p) where {FT<:AbstractFloat}

Returns the albedo reset rate to return the albedo value to that dictated by the
AbstractResetAlbedoModel for the next time step. Requires definition of the method
get_initial_albedo(alb_reset_model::T, Y, p) where T is the desired type of AbstractResetAlbedoModel.
"""
function albedo_reset_rate(
    alb_reset_model::AbstractResetAlbedoModel,
    Δt::FT,
    Y,
    p,
) where {FT <: AbstractFloat}
    return (get_initial_albedo(alb_reset_model, Y, p) .- Y.snow.A) ./ Δt
end

"""
    update_dzdt!(dzdt, density::NeuralDepthModel, model, Y, p)

Updates the dY.snow.Z field in places with the predicted change in snow depth (rate) given the model state `Y` and the `NeuralDepthModel`
density paramterization.
"""
function update_dzdt!(dzdt, density::NeuralDepthModel, model, Y, p)
    #=
    Technically, using the EMA vars for inputs directly like below does not account
    for an "off-by-1" lag error within them; they do not include the "present"
    information by the time they are used as an input. For example, during step N,
    when `Y.snow.P_avg` is called below, it stores the value
        (w*P_snow(N-1) + w*(1-w)*P_snow(N-2) + ... + (1-w)^(N-1) * initial_value_of_driver_P_snow);
    it does not include the value P_snow(N) = p.drivers.P_snow.

    For an Euler-stepping scheme (1 evaluation/update per time-step), we can definitely
    "fix" this "off-by-1" lag by passing a modified input for each EMA variable,
    e.g. instead pass the network the input
        w*p.drivers.P_snow + (1-w)*(Y.snow.P_avg)

    But, since the prognostic variables update for sub-steps within a non-Euler
    stepping scheme, I am unsure if this would work on ARS111 or other schemes
    (that is, does it over-bias/over-weight the contribution of p.drivers.P_snow?).
    This would be worth checking/verifying if Δt was more than an hour or so,
    since beyond that Δt magnitude, the off-by-1 lag would be significant, and 
    Y.snow.Z updates would be effectively delayed compared to the timing of
    snowfall. But, since the soil model requires time-steps on the order ≤ 15 min
    for stability, this lag is numerically negligible, and we save computational
    resources by not doing the extra modification, so whether the modification
    works for sub-stepping has been abandoned in favor of computational efficiency.

    When this was first implemented, using model callbacks inside the simulation
    did not exist - to completely mitigate the above concern, one could shift
    the calculation of the 24-hour average to a model-specific callback variable
    that only updates once per true time-step instead of during sub-timesteps,
    and this would also circumvent the need for modiciation/extra computations.
    Ideally, one would use the moving-window past-24-hour average including the
    present value, but this would require a circular buffer of size dependent on
    Δt; an EMA gives beneficial weighting to immediate values while only requiring
    storage of a single additional field.
    =#
    dzdt .=
        eval_z_nn.(
            Ref(density.z_model),
            p.snow.z_snow,
            Y.snow.S,
            Y.snow.P_avg,
            Y.snow.T_avg,
            Y.snow.R_avg,
            Y.snow.Qrel_avg,
            Y.snow.u_avg,
        )
end

"""
    update_dαdt!(dαdt, albedo::NeuralAlbedoModel, model, Y, p)

Updates the dY.snow.A field in places with the predicted change in snow albedo given the model state `Y` and the `NeuralAlbedoModel`
albedo paramterization.
"""
function update_dαdt!(dαdt, albedo::NeuralAlbedoModel, model, Y, p, t)
    start_date = model.boundary_conditions.radiation.start_date
    dαdt .=
        reset_alb_tendency.(
            Y.snow.S,
            abs.(p.drivers.P_snow),
            albedo.new_alb.no_snowpack_S,
            albedo.new_alb.reset_limit,
            albedo_reset_rate(albedo.new_alb, model.parameters.Δt, Y, p),
            eval_alb_nn.(
                Ref(albedo.alb_model),
                Y.snow.A,
                abs.(p.drivers.P_snow),
                albedo.za_solarnoon(t, start_date),
                p.drivers.T .- model.parameters.earth_param_set.T_freeze,
            ),
        )
end

"""
    clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT, ρ_min_frac::FT)::FT

A helper function which clips the tendency of Z such that
its behavior is consistent with that of S: if all snow melts
within a timestep, we clip the tendency of S so that it does
not become negative, and here we also clip the tendency of Z
so that depth does not become negative. Additionally, if the
tendencies of Z and S are such that we would encounter Z < S
(rho_snow > rho_liq), we also clip the tendency.
"""
function clip_dZdt(
    S::FT,
    Z::FT,
    dSdt::FT,
    dZdt::FT,
    Δt::FT,
    ρ_min_frac::FT,
)::FT where {FT}
    #This function is simple for now since we let Y.snow.Z also be the per-ground-area value, just like Y.snow.S is
    new_S_ground_area = dSdt * Δt + S
    new_Z = dZdt * Δt + Z #also ground-area presently
    #new_scf = (something), if we pivot to let Y.snow.Z be the per-snow-area value

    new_safe_z = clamp(new_Z, new_S_ground_area, new_S_ground_area / ρ_min_frac) #new_S_ground_area -> (new_S_ground_area / new_scf) if if Y.snow.Z becomes per-snow-area, which can be unstable
    return (new_safe_z - Z) / Δt
end

"""
    compute_extra_prog_tendency!(density::NeuralDepthModel, model::SnowModel, dY, Y, p, t)

Updates all prognostic variables associated with density/depth given the current model state and the `NeuralDepthModel`
density paramterization.
"""
function ClimaLand.Snow.compute_extra_prog_tendency!(
    density::NeuralDepthModel,
    model::SnowModel,
    dY,
    Y,
    p,
    t,
)
    prog_comps = model.boundary_conditions.prognostic_land_components
    update_dzdt!(dY.snow.Z, density, model, Y, p)

    # Now clip the tendency so that Z stays within approximately physical bounds.
    @. dY.snow.Z = clip_dZdt(
        Y.snow.S,
        Y.snow.Z,
        dY.snow.S, #assumes dY.snow.S is updated (and clipped) before dY.snow.Z
        dY.snow.Z,
        model.parameters.Δt,
        density.ρ_min_frac,
    )

    @. dY.snow.P_avg =
        ema_tendency(density.w, abs(p.drivers.P_snow), Y.snow.P_avg)
    @. dY.snow.T_avg = ema_tendency(
        density.w,
        p.drivers.T - model.parameters.earth_param_set.T_freeze,
        Y.snow.T_avg,
    )
    dY.snow.R_avg .=
        ema_tendency.(density.w, snow_SW_down(Val(prog_comps), p), Y.snow.R_avg)
    dY.snow.Qrel_avg .=
        ema_tendency.(
            density.w,
            Thermodynamics.relative_humidity.(
                model.boundary_conditions.atmos.thermo_params,
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
            ),
            Y.snow.Qrel_avg,
        )
    @. dY.snow.u_avg = ema_tendency(density.w, p.drivers.u, Y.snow.u_avg)
end

"""
    compute_extra_prog_tendency!(albedo::NeuralAlbedoModel, model::SnowModel, dY, Y, p, t)

Updates all prognostic variables associated with snow albedo given the current model state and the `NeuralAlbedoModel`
albedo paramterization.
"""
function ClimaLand.Snow.compute_extra_prog_tendency!(
    albedo::NeuralAlbedoModel,
    model::SnowModel,
    dY,
    Y,
    p,
    t,
)
    update_dαdt!(dY.snow.A, albedo, model, Y, p, t)
end

end

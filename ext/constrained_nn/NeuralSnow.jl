module NeuralSnow
using ClimaLand
using ClimaLand.Snow:
    AbstractDensityModel,
    AbstractAlbedoModel,
    SnowModel,
    SnowParameters,
    snow_bulk_density
import ClimaComms
import ClimaLand.Snow:
    extra_prog_vars,
    extra_prog_types,
    extra_prog_names,
    update_density_and_depth!,
    update_snow_albedo!,
    update_extra_prog!
import ClimaLand.Parameters as LP
import ClimaCore
using Thermodynamics

using Flux: relu
using StaticArrays
using ..ConstrainedNeuralModels

export NeuralDepthModel, NeuralAlbedoModel

"""
   AlbedoInitialValueModel{FT<:AbstractFloat}

Defines the parameterization type for getting the initial value of a new snowpack, to
be combined with the prognostic update dynamics of a NeuralAlbedoModel.

An instance of this type must have at least:
    1) one functor method defined, (::T)(model, Y, p)
       for T <: AlbedoInitialValueModel that returns the initial albedo value,
       which will be called in `albedo_reset_rate()`
    2) a float "reset_limit" (m/s) of critical preciptiation rate,
       which will also get used in resetting the snow albedo.
"""
abstract type AlbedoInitialValueModel{FT <: AbstractFloat} end

"""
    ConstantInitialAlbedo{FT} <: AlbedoInitialValueModel{FT}
Establishes the parameterization for initial values of albedo in a new snowpack,
where the dynamics of the model are begun from a specified constant value.
Default values are a new albedo of 0.8 for 3 mm liquid-equivalent of snow precipitation in 24 hours.
"""
struct ConstantInitialAlbedo{FT} <: AlbedoInitialValueModel{FT}
    "Initial albedo constant (unitless 0-1)"
    init_α::FT
    "Threshold of liquid-equiavlent snow precipitation for a reset (m/s)"
    reset_limit::FT
end
ConstantInitialAlbedo(::Type{FT}) where {FT <: AbstractFloat} =
    ConstantInitialAlbedo(FT(0.8), FT(0.003 / 86400))

"""
    (m::ConstantInitialAlbedo)(model, Y, p) where {FT}

Defines the behavior for a ConstantInitialAlbedo parameterization,
which is called within `albedo_reset_rate()` and combined with the
NeuralAlbedoModel output in `update_dαdt!`.
"""
function (m::ConstantInitialAlbedo)(model, Y, p)
    return m.init_α
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
    function NeuralDepthModel(z_model::MD, w::FT) where {FT, MD}
        static_z_model = make_static_model(z_model, skip_check = true)
        new{FT, typeof(static_z_model)}(static_z_model, w)
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
    IA <: AlbedoInitialValueModel,
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
    NeuralDepthModel(FT::Type{<:AbstractFloat}; Δt::Union{Nothing, FT}=nothing, w::Union{FT, Nothing})
An outer constructor for the NeuralDepthModel density parameterization for usage in a snow model.

This model utilizes the neural network formulated in Charbonneau et. al. (2025) to predict the rate of change of snow height dz/dt.
Since the default network receives inputs representing the average value per feature per day (86400 s), an exponential-weighting
scheme is used to track the moving average of the input fields. A default weighting is supplied
as `w = 2/86400 s⁻¹` following common practice of choosing a weighting of 2/(desired temporal window size),
though the user may specify their own value by keyword argument. Note that this default value is
unstable (whenever `1/w > model timestep`) for model timesteps larger than half a day (43200 seconds).
The user should also specify a model timestep Δt in seconds which sets the model scaling in order to guarantee the snow depth remains nonnegative.
"""
function NeuralDepthModel(
    ::Type{FT};
    model_params::Union{Vector{<:Real}, Nothing} = nothing,
    Δt::Union{<:Real, Nothing} = nothing,
    w::Union{<:Real, Nothing} = nothing,
) where {FT <: AbstractFloat}
    usemodel = get_znetwork(FT, params = model_params)
    weight = !isnothing(w) ? FT(w) : FT(2 / 86400)
    if !isnothing(Δt)
        _set_timescale!(usemodel, Δt) #assumes Δt is provided in seconds
        if (Δt > 43200) && isnothing(w)
            error("Please supply a weight for the
                   exponential moving average, the
                   default value is unstable in this case.")
        end
    end
    return NeuralDepthModel(usemodel, weight)
end

"""
    NeuralAlbedoModel(
        ::Type{FT},
        surface_space,
        earth_param_set;
        initial_albedo = ConstantInitialAlbedo(FT),
        model_params::Union{Vector{<:Real}, Nothing} = nothing,
        Δt::Union{<:Real, Nothing} = nothing,
    ) where {FT <: AbstractFloat}

An outer constructor for the NeuralAlbedoModel albedo parameterization for usage in a snow model. It requires model domains that 
contain latitude and longitude information, requiring the passing of the model's surface space to the constructor.

This model utilizes the neural network develped by Charbonneau et. al. (2025) to predict the rate of change of snow albedo dα/dt.
The user should also specify a model timestep Δt in seconds which sets the model scaling in order to guarantee the snow albedo remains physical.
If unspecified, the initial albedo value of a new snowpack is dictated by a ConstantInitialAlbedo parameterization of 0.8.
"""
function NeuralAlbedoModel(
    ::Type{FT},
    surface_space,
    earth_param_set;
    initial_albedo = ConstantInitialAlbedo(FT),
    model_params::Union{Vector{<:Real}, Nothing} = nothing,
    Δt::Union{<:Real, Nothing} = nothing,
) where {FT <: AbstractFloat}
    usemodel = get_albnetwork(FT, params = model_params)
    !isnothing(Δt) && _set_timescale!(usemodel, Δt)
    model_lat = ClimaCore.Fields.coordinate_field(surface_space).lat
    model_lon = ClimaCore.Fields.coordinate_field(surface_space).long
    insol_ps = earth_param_set.insol_params
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

ClimaLand.Snow.extra_prog_names(::NeuralDepthModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)
ClimaLand.Snow.extra_prog_names(::NeuralAlbedoModel) = (:surface,)

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
    # ^this is not technically true, but creates a lot of issues trying to rectify it, as there is not a stable
    # or physically consistent formula to get swe-per-snow-area now that we've let Y.snow.S be the per-ground-area.
    @. z_snow = Y.snow.Z #p.snow.z_snow is z-per-ground-area, so we set its value directly. Need to change to "Y.snow.Z * scf" if we change Y.snow.Z to represent true snow depth.
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
    @. α = Y.snow.A
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
    #that is physically consistent in the small snowpack limit AND doesn't:
    # - either make the snowpack start to grow again as scf gets smaller faster than z or SWE gets smaller
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
    ema_step(w::FT, new::FT, old::FT)::FT where {FT}

A utility function to get the derivative value of an exponential moving average (EMA),
in terms of its old value, the incoming value, and the scaling weight.
"""
function ema_step(w::FT, new::FT, old::FT)::FT where {FT}
    return w * (new - old)
end

"""
    ema_value(w::FT, new::FT, old::FT)::FT where {FT}

A utility function to get the value of an exponential moving average (EMA),
in terms of its old value, the incoming value, and the scaling weight.
"""
function ema_value(w::FT, new::FT, old::FT)::FT where {FT}
    return w * new + (1 - w) * old
end

"""
    snow_surf_SW_down(mtype::Val{(:snow,)}, p)

Obtains the correct downward shortwave radiation magnitude for a snow model without an
overlying canopy (standalone or just coupled with soil).
"""
@inline snow_surf_SW_down(mtype::Union{Val{(:snow,)}, Val{(:snow, :soil)}}, p) =
    p.drivers.SW_d

"""
    snow_surf_SW_down(mtype::Union{
        Val{(:canopy, :snow, :soil, :soilco2)},
        Val{(:canopy, :snow, :soil)},
        }, p)

Obtains the correct downward shortwave radiation magnitude for a snow model with
an overlying canopy.
"""
@inline function snow_surf_SW_down(
    mtype::Union{
        Val{(:canopy, :snow, :soil, :soilco2)},
        Val{(:canopy, :snow, :soil)},
    },
    p,
)
    return p.drivers.SW_d
    #=@. -(
        p.canopy.radiative_transfer.nir.trans *
        p.canopy.radiative_transfer.nir_d *
        (1 - p.snow.α_snow) +
        p.canopy.radiative_transfer.par.trans *
        p.canopy.radiative_transfer.par_d *
        (1 - p.snow.α_snow)
    )
        =#
end

"""
    reset_condition(S::FT, P::FT, dt::FT, lim::FT, reset_rate::FT, prognostic_rate::FT)::FT where {FT}

Handles albedo model prognostic updates, returning the albedo to the initial value
as dicated by an AlbedoInitialValueModel if a new snowpack is being formed, and
otherwise returning the rate as dicated by the prognostic model.
"""
function reset_condition(
    S::FT,
    P::FT,
    lim::FT,
    reset_rate::FT,
    prognostic_rate::FT,
)::FT where {FT}
    #Do I need to make 0.001 a parameter too? could use parameters.ΔS, but that's pretty big comparatively?
    return (S < 0.001) && (P > lim) ? reset_rate : prognostic_rate
end

"""
    albedo_reset_rate(init_alb::AlbedoInitialValueModel, Δt::FT, Y, p) where {FT<:AbstractFloat}

Returns the albedo reset rate to return the albedo value to that dictated by the
AlbedoInitialValueModel for the next time step. Requires definition of the functor
method (::T)(model, Y, p) where T is the desired type of AlbedoInitialValueModel.
"""
function albedo_reset_rate(
    init_alb::AlbedoInitialValueModel,
    model,
    Δt::FT,
    Y,
    p,
) where {FT <: AbstractFloat}
    return (init_alb(model, Y, p) .- Y.snow.A) ./ Δt
end

"""
    update_dzdt!(dzdt, density::NeuralDepthModel, model, Y, p)

Updates the dY.snow.Z field in places with the predicted change in snow depth (rate) given the model state `Y` and the `NeuralDepthModel`
density paramterization.
"""
function update_dzdt!(dzdt, density::NeuralDepthModel, model, Y, p)
    #=
    the input for the step from time t to t+1 should include the variable information
    from time t (especially considering precipiation), even though the EMA prognostic vars will not have been
    updated with the dY information from time t yet. This means we need
    to get that value here:
    =#
    prog_comp = model.boundary_conditions.prognostic_land_components
    dzdt .=
        eval_z_nn.(
            Ref(density.z_model),
            Y.snow.Z,
            Y.snow.S,
            ema_value.(density.w, abs.(p.drivers.P_snow), Y.snow.P_avg),
            ema_value.(
                density.w,
                p.drivers.T .- model.parameters.earth_param_set.T_freeze,
                Y.snow.T_avg,
            ),
            ema_value.(
                density.w,
                snow_surf_SW_down(Val(prog_comp), p),
                Y.snow.R_avg,
            ),
            ema_value.(
                density.w,
                Thermodynamics.relative_humidity.(
                    model.boundary_conditions.atmos.thermo_params,
                    p.drivers.T,
                    p.drivers.P,
                    p.drivers.q,
                ),
                Y.snow.Qrel_avg,
            ),
            ema_value.(density.w, p.drivers.u, Y.snow.u_avg),
        )
end

"""
    update_dαdt!(dαdt, albedo::NeuralAlbedoModel, model, Y, p)

Updates the dY.snow.A field in places with the predicted change in snow albedo given the model state `Y` and the `NeuralAlbedoModel`
albedo paramterization.
"""
function update_dαdt!(dαdt, albedo::NeuralAlbedoModel, model, Y, p, t)
    start_date = model.boundary_conditions.radiation.start_date #is there a better way to reference this?
    dαdt .=
        reset_condition.(
            Y.snow.S,
            abs.(p.drivers.P_snow),
            albedo.new_alb.reset_limit,
            albedo_reset_rate(albedo.new_alb, model, model.parameters.Δt, Y, p),
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
    #This function is simple for now since we let Y.snow.Z also be the per-ground-area value, just like Y.snow.S is
    #Unsure how to edit this function if Y.snow.Z becomes per-snow-area:
    new_S_ground_area = dSdt * Δt + S
    new_Z = dZdt * Δt + Z #also ground-area
    #new_scf = (something), if we pivot to let Y.snow.Z be the per-snow-area value

    #Handle cases if S is set to 0, or if Z would have been set to Z < S:
    if new_S_ground_area <= 0 #0, or eps(FT)?
        return -Z / Δt
    elseif new_Z < new_S_ground_area # new_Z * new_scf < new_S, if Y.snow.Z becomes per-snow-area:
        return (new_S_ground_area - Z) / Δt  # (new_S / new_scf - Z) / Δt, if Y.snow.Z becomes per-snow-area, which can be unstable
    else
        return dZdt
    end
end

"""
    update_extra_prog!(density::NeuralDepthModel, model::SnowModel, dY, Y, p, t)

Updates all prognostic variables associated with density/depth given the current model state and the `NeuralDepthModel`
density paramterization.
"""
function ClimaLand.Snow.update_extra_prog!(
    density::NeuralDepthModel,
    model::SnowModel,
    dY,
    Y,
    p,
    t,
)
    update_dzdt!(dY.snow.Z, density, model, Y, p)

    # Now we clip the tendency so that Z stays within approximately physical bounds.
    @. dY.snow.Z = clip_dZdt(
        Y.snow.S,
        Y.snow.Z,
        dY.snow.S, #assumes dY.snow.S is updated (and clipped) before dY.snow.Z
        dY.snow.Z,
        model.parameters.Δt,
    )

    @. dY.snow.P_avg = ema_step(density.w, abs(p.drivers.P_snow), Y.snow.P_avg)
    @. dY.snow.T_avg = ema_step(
        density.w,
        p.drivers.T - model.parameters.earth_param_set.T_freeze,
        Y.snow.T_avg,
    )
    @. dY.snow.R_avg = ema_step(density.w, p.drivers.SW_d, Y.snow.R_avg)
    dY.snow.Qrel_avg .=
        ema_step.(
            density.w,
            Thermodynamics.relative_humidity.(
                model.boundary_conditions.atmos.thermo_params,
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
            ),
            Y.snow.Qrel_avg,
        )
    @. dY.snow.u_avg = ema_step(density.w, p.drivers.u, Y.snow.u_avg)
end

"""
    update_extra_prog!(albedo::NeuralAlbedoModel, model::SnowModel, dY, Y, p, t)

Updates all prognostic variables associated with snow albedo given the current model state and the `NeuralAlbedoModel`
albedo paramterization.
"""
function ClimaLand.Snow.update_extra_prog!(
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

#fix tests for lat/lon constructor and add tests for snow_surf_SW_down
# try clamping without complex canopy value (fix it after that test is done), see if NaNs go away

# sodankyla latitude is wrong in the artifact
# how to make the sinh()/model NaN safe? / does model not work at really-high latitude sites where solar zenith angles don't do what's expected in training data?

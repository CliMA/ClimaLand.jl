module InlandWater

using ClimaLand
using LazyBroadcast: lazy
using DocStringExtensions
using ClimaCore
import ClimaCore: Fields, Spaces
import ..Parameters as LP
import ClimaParams as CP
using Thermodynamics
using SurfaceFluxes
import ClimaLand:
    AbstractExpModel,
    make_update_aux,
    make_compute_exp_tendency,
    make_update_boundary_fluxes,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    prognostic_domain_names,
    auxiliary_types,
    auxiliary_domain_names,
    component_temperature,
    component_specific_humidity,
    get_update_surface_humidity_function,
    get_update_surface_temperature_function,
    surface_height,
    surface_displacement_height,
    surface_roughness_model,
    get_∂T_sfc∂T_function,
    get_∂q_sfc∂T_function,
    surface_albedo,
    surface_emissivity,
    turbulent_fluxes!,
    get_drivers,
    get_earth_param_set,
    total_liq_water_vol_per_area!,
    total_energy_per_area!,
    return_momentum_fluxes
using ClimaLand:
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    AbstractBC
import ClimaLand.Domains: surface_to_subsurface
import ClimaLand.Domains:
    Point, Column, SphericalShell, SphericalSurface, HybridBox, Plane
using Interpolations
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput

export SlabLakeParameters,
    SlabLakeModel,
    AbstractInlandWaterModel,
    AtmosDrivenLakeBC,
    lake_boundary_fluxes!,
    inland_water_mask

"""
    AbstractInlandWaterModel{FT} <: ClimaLand.AbstractExpModel{FT}

Abstract type for inland water models.
"""
abstract type AbstractInlandWaterModel{FT} <: ClimaLand.AbstractExpModel{FT} end

"""
    SlabLakeParameters{FT, PSE}

Fixed-depth slab lake parameters. The slab lake keeps a constant depth
while evolving its energy prognostically.

The `conductance` G (W m⁻² K⁻¹) represents the lake-side thermal
conductance, i.e. G = κ_lake / depth, combining lake thermal conductivity
and depth into a single tunable parameter. The sediment heat flux is:
    Q_sed = -G_eff * (T_lake - T_soil)
where
    G_eff = 1 / (1/G + Δz_soil/(2*κ_soil))

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SlabLakeParameters{FT <: AbstractFloat, PSE}
    "Mixed-layer depth of the slab lake (m)"
    depth::FT = FT(10)
    "Lake–sediment conductance (W m⁻² K⁻¹)"
    conductance::FT = FT(0.1)
    "Open-water albedo"
    liquid_albedo::FT = FT(0.08)
    "Lake-ice albedo"
    ice_albedo::FT = FT(0.6)
    "Lake emissivity"
    emissivity::FT = FT(0.97)
    "Roughness length for momentum over lake water/ice (m)"
    z_0m::FT = FT(0.001)
    "Roughness length for scalars over lake water/ice (m)"
    z_0b::FT = FT(0.0001)
    "Physical constants and clima-wide parameters"
    earth_param_set::PSE
end

Base.broadcastable(ps::SlabLakeParameters) = tuple(ps)

"""
    SlabLakeParameters{FT}(; earth_param_set, kwargs...)

Convenience constructor that infers the PSE type parameter.
"""
function SlabLakeParameters{FT}(;
    earth_param_set::PSE,
    kwargs...,
) where {FT, PSE}
    return SlabLakeParameters{FT, PSE}(; earth_param_set, kwargs...)
end

"""
    SlabLakeParameters(
        toml_dict::CP.ParamDict;
        depth = toml_dict["slab_lake_depth"],
        conductance = toml_dict["slab_lake_conductance"],
        liquid_albedo = toml_dict["slab_lake_liquid_albedo"],
        ice_albedo = toml_dict["slab_lake_ice_albedo"],
        emissivity = toml_dict["slab_lake_emissivity"],
        z_0m = toml_dict["slab_lake_z_0m"],
        z_0b = toml_dict["slab_lake_z_0b"],
    )

TOML dictionary constructor for `SlabLakeParameters`.
"""
function SlabLakeParameters(
    toml_dict::CP.ParamDict;
    depth = toml_dict["slab_lake_depth"],
    conductance = toml_dict["slab_lake_conductance"],
    liquid_albedo = toml_dict["slab_lake_liquid_albedo"],
    ice_albedo = toml_dict["slab_lake_ice_albedo"],
    emissivity = toml_dict["slab_lake_emissivity"],
    z_0m = toml_dict["slab_lake_z_0m"],
    z_0b = toml_dict["slab_lake_z_0b"],
)
    FT = CP.float_type(toml_dict)
    earth_param_set = LP.LandParameters(toml_dict)
    return SlabLakeParameters{FT, typeof(earth_param_set)}(;
        depth,
        conductance,
        liquid_albedo,
        ice_albedo,
        emissivity,
        z_0m,
        z_0b,
        earth_param_set,
    )
end

"""
    AtmosDrivenLakeBC{A, B, C} <: ClimaLand.AbstractBC

Boundary condition for the slab lake driven by atmospheric forcing.
Same pattern as `AtmosDrivenSnowBC`.

$(DocStringExtensions.FIELDS)
"""
struct AtmosDrivenLakeBC{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
    C <: Tuple,
} <: ClimaLand.AbstractBC
    "Atmospheric forcing"
    atmos::A
    "Radiative forcing"
    radiation::B
    "Prognostic land components present"
    prognostic_land_components::C
end

"""
    SlabLakeModel{FT, PS, D, BC, M} <: AbstractInlandWaterModel{FT}

A slab lake model for inland water points.

$(DocStringExtensions.FIELDS)
"""
struct SlabLakeModel{FT, PS, D, BC, M, DZ} <: AbstractInlandWaterModel{FT}
    "Lake parameters"
    parameters::PS
    "The lake domain (surface only)"
    domain::D
    "Boundary conditions"
    boundary_conditions::BC
    "Inland water mask on the surface space (fraction 0–1)"
    inland_water_mask::M
    "Top soil layer center-to-surface distance (from soil domain)"
    Δz_top::DZ
end

"""
    SlabLakeModel{FT}(;
        parameters,
        domain,
        boundary_conditions,
        inland_water_mask,
        Δz_top,
    )

Construct a `SlabLakeModel`. The `Δz_top` field should be the center-to-surface
distance of the top soil layer (from the soil domain), used for computing
the sediment heat flux.
"""
function SlabLakeModel{FT}(;
    parameters::PS,
    domain::D,
    boundary_conditions::BC,
    inland_water_mask::M,
    Δz_top::DZ,
) where {FT, PS, D, BC, M, DZ}
    return SlabLakeModel{FT, PS, D, BC, M, DZ}(
        parameters,
        domain,
        boundary_conditions,
        inland_water_mask,
        Δz_top,
    )
end

"""
    SlabLakeModel(
        FT,
        domain,
        forcing,
        toml_dict::CP.ParamDict,
        Δz_top;
        inland_water_mask,
        prognostic_land_components = (:lake,),
        parameters = SlabLakeParameters(toml_dict),
    )

Convenience constructor for `SlabLakeModel` that builds the model from a
TOML parameter dictionary, following the same pattern as `SnowModel`.
"""
function SlabLakeModel(
    FT,
    domain,
    forcing,
    toml_dict::CP.ParamDict,
    Δz_top;
    inland_water_mask,
    prognostic_land_components = (:lake,),
    parameters = SlabLakeParameters(toml_dict),
)
    boundary_conditions = AtmosDrivenLakeBC(
        forcing.atmos,
        forcing.radiation,
        prognostic_land_components,
    )
    return SlabLakeModel{
        FT,
        typeof(parameters),
        typeof(domain),
        typeof(boundary_conditions),
        typeof(inland_water_mask),
        typeof(Δz_top),
    }(
        parameters,
        domain,
        boundary_conditions,
        inland_water_mask,
        Δz_top,
    )
end

ClimaLand.name(::SlabLakeModel) = :lake
ClimaLand.get_earth_param_set(model::SlabLakeModel) =
    model.parameters.earth_param_set

# ─── Lake physics ─────────────────────────────────────────────────────────────

function lake_energy_at_freezing(
    q_l::FT,
    params::SlabLakeParameters{FT},
    earth_param_set,
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _T_freeze = FT(LP.T_freeze(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    _cp_i = FT(LP.cp_i(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))
    cp = q_l * _cp_l + (FT(1) - q_l) * _cp_i
    return _ρ_l * params.depth * cp * (_T_freeze - _T_ref) -
           _ρ_l * params.depth * (FT(1) - q_l) * _LH_f0
end

function lake_liquid_fraction(
    U::FT,
    params::SlabLakeParameters{FT},
    earth_param_set,
) where {FT}
    U_ice = lake_energy_at_freezing(FT(0), params, earth_param_set)
    U_liq = lake_energy_at_freezing(FT(1), params, earth_param_set)
    if U <= U_ice
        return FT(0)
    elseif U >= U_liq
        return FT(1)
    else
        return (U - U_ice) / (U_liq - U_ice)
    end
end

function lake_energy_from_temperature(
    T::FT,
    params::SlabLakeParameters{FT},
    earth_param_set,
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _T_freeze = FT(LP.T_freeze(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    _cp_i = FT(LP.cp_i(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))
    if T <= _T_freeze
        return _ρ_l * params.depth * _cp_i * (T - _T_ref) -
               _ρ_l * params.depth * _LH_f0
    else
        return _ρ_l * params.depth * _cp_l * (T - _T_ref)
    end
end

function lake_temperature(
    U::FT,
    q_l::FT,
    params::SlabLakeParameters{FT},
    earth_param_set,
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _T_freeze = FT(LP.T_freeze(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    _cp_i = FT(LP.cp_i(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))
    depth = params.depth
    if q_l < eps(FT)
        return _T_ref + (U + _ρ_l * depth * _LH_f0) / (_ρ_l * depth * _cp_i)
    elseif q_l > FT(1) - eps(FT)
        return _T_ref + U / (_ρ_l * depth * _cp_l)
    else
        return _T_freeze
    end
end

function lake_volumetric_internal_energy(
    T::FT,
    q_l::FT,
    earth_param_set,
) where {FT}
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    _cp_i = FT(LP.cp_i(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))
    cp = q_l * _cp_l + (FT(1) - q_l) * _cp_i
    return _ρ_l * (cp * (T - _T_ref) - (FT(1) - q_l) * _LH_f0)
end

function lake_runoff_energy_flux(
    runoff::FT,
    T::FT,
    q_l::FT,
    earth_param_set,
) where {FT}
    return runoff * lake_volumetric_internal_energy(T, q_l, earth_param_set)
end

lake_surface_albedo(q_l::FT, params::SlabLakeParameters{FT}) where {FT} =
    q_l * params.liquid_albedo + (FT(1) - q_l) * params.ice_albedo

"""
    lake_sediment_heat_flux(T_lake, T_soil, κ_soil, Δz_soil, params)

Compute sediment heat flux using series-conductance model:
    Q_sed = -G_eff * (T_lake - T_soil)
where G_eff = 1 / (1/G + Δz_soil/(2*κ_soil)).
"""
function lake_sediment_heat_flux(
    T_lake::FT,
    T_soil::FT,
    κ_soil::FT,
    Δz_soil::FT,
    params::SlabLakeParameters{FT},
) where {FT}
    G = params.conductance
    G_eff = FT(1) / (FT(1) / G + Δz_soil / (FT(2) * κ_soil))
    return -G_eff * (T_lake - T_soil)
end

function lake_specific_humidity(
    T_sfc::FT,
    q_l::FT,
    ρ_sfc::FT,
    earth_param_set,
) where {FT}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    qsat_ice = Thermodynamics.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Ice(),
    )
    qsat_liq = Thermodynamics.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Liquid(),
    )
    return (FT(1) - q_l) * qsat_ice + q_l * qsat_liq
end

# ─── Interface: prognostic / auxiliary variable declarations ──────────────────

ClimaLand.prognostic_vars(::SlabLakeModel) = (:U,)
ClimaLand.prognostic_types(::SlabLakeModel{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(::SlabLakeModel) = (:surface,)

function ClimaLand.auxiliary_vars(::SlabLakeModel)
    return (
        :T,
        :q_l,
        :sediment_heat_flux,
        :surface_water_flux,
        :surface_energy_flux,
        :runoff_energy_flux,
        :runoff,
        :lake_fraction,
        :lake_albedo,
        :turbulent_fluxes,
        :R_n,
    )
end

function ClimaLand.auxiliary_types(::SlabLakeModel{FT}) where {FT}
    return (
        FT, # T
        FT, # q_l
        FT, # sediment_heat_flux
        FT, # surface_water_flux
        FT, # surface_energy_flux
        FT, # runoff_energy_flux
        FT, # runoff
        FT, # lake_fraction
        FT, # lake_albedo
        NamedTuple{
            (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T),
            Tuple{FT, FT, FT, FT, FT},
        }, # turbulent_fluxes
        FT, # R_n
    )
end

function ClimaLand.auxiliary_domain_names(::SlabLakeModel)
    return (
        :surface, # T
        :surface, # q_l
        :surface, # sediment_heat_flux
        :surface, # surface_water_flux
        :surface, # surface_energy_flux
        :surface, # runoff_energy_flux
        :surface, # runoff
        :surface, # lake_fraction
        :surface, # lake_albedo
        :surface, # turbulent_fluxes
        :surface, # R_n
    )
end

# ─── Surface interface methods ────────────────────────────────────────────────

ClimaLand.component_temperature(::SlabLakeModel, Y, p) = p.lake.T

function ClimaLand.surface_emissivity(model::SlabLakeModel{FT}, Y, p) where {FT}
    return model.parameters.emissivity
end

function ClimaLand.surface_albedo(model::SlabLakeModel, Y, p)
    return p.lake.lake_albedo
end

function ClimaLand.surface_roughness_model(
    model::SlabLakeModel{FT},
    Y,
    p,
) where {FT}
    ps = model.parameters
    return SurfaceFluxes.ConstantRoughnessParams{FT}(ps.z_0m, ps.z_0b)
end

function ClimaLand.surface_height(model::SlabLakeModel{FT}, Y, p) where {FT}
    return FT(0)
end

function ClimaLand.surface_displacement_height(
    model::SlabLakeModel{FT},
    Y,
    p,
) where {FT}
    return FT(0)
end

function ClimaLand.component_specific_humidity(
    model::SlabLakeModel{FT},
    Y,
    p,
) where {FT}
    earth_param_set = model.parameters.earth_param_set
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    T_sfc = p.lake.T
    h_sfc = surface_height(model, Y, p)
    atmos = model.boundary_conditions.atmos
    ρ_sfc = @. lazy(
        ClimaLand.compute_ρ_sfc(
            surface_flux_params,
            p.drivers.T,
            p.drivers.P,
            p.drivers.q,
            atmos.h - h_sfc,
            T_sfc,
        ),
    )
    return @. lazy(
        lake_specific_humidity(T_sfc, p.lake.q_l, ρ_sfc, earth_param_set),
    )
end

# ─── Update aux ───────────────────────────────────────────────────────────────

function ClimaLand.make_update_aux(model::SlabLakeModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        params = model.parameters
        earth_param_set = params.earth_param_set
        mask = model.inland_water_mask

        # Lake fraction from mask
        @. p.lake.lake_fraction = mask

        # Thermodynamic state from prognostic U
        @. p.lake.q_l = lake_liquid_fraction(Y.lake.U, params, earth_param_set)
        @. p.lake.T =
            lake_temperature(Y.lake.U, p.lake.q_l, params, earth_param_set)
        @. p.lake.lake_albedo = lake_surface_albedo(p.lake.q_l, params)

        # Sediment heat flux (needs soil κ and T from p.soil)
        T_soil_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.T)
        κ_soil_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.κ)
        Δz_top = model.Δz_top
        @. p.lake.sediment_heat_flux =
            p.lake.lake_fraction * lake_sediment_heat_flux(
                p.lake.T,
                T_soil_sfc,
                κ_soil_sfc,
                Δz_top,
                params,
            )
    end
    return update_aux!
end

# ─── Boundary fluxes ──────────────────────────────────────────────────────────

function ClimaLand.make_update_boundary_fluxes(
    model::SlabLakeModel{FT},
) where {FT}
    function update_boundary_fluxes!(p, Y, t)
        lake_boundary_fluxes!(
            model.boundary_conditions,
            Val(model.boundary_conditions.prognostic_land_components),
            model,
            Y,
            p,
            t,
        )
    end
    return update_boundary_fluxes!
end

"""
    lake_boundary_fluxes!(bc, prognostic_land_components, model, Y, p, t)

Compute the lake surface boundary fluxes.
In standalone mode, the lake computes its own turbulent fluxes and radiation.
In integrated mode, this is overridden by the LandModel.
"""
function lake_boundary_fluxes!(
    bc::AtmosDrivenLakeBC,
    prognostic_land_components,
    model::SlabLakeModel{FT},
    Y,
    p,
    t,
) where {FT}
    # Compute turbulent fluxes using the generic interface
    turbulent_fluxes!(p.lake.turbulent_fluxes, bc.atmos, model, Y, p, t)

    # Compute net radiation in standalone mode
    ClimaLand.net_radiation!(p.lake.R_n, bc.radiation, model, Y, p, t)

    earth_param_set = model.parameters.earth_param_set
    f = p.lake.lake_fraction

    # Surface water flux: precip + evaporation
    @. p.lake.surface_water_flux =
        f * (p.drivers.P_liq + p.lake.turbulent_fluxes.vapor_flux)

    # Runoff = negative of net water flux (excess leaves as runoff)
    @. p.lake.runoff = -p.lake.surface_water_flux

    # Runoff energy flux
    @. p.lake.runoff_energy_flux =
        f * lake_runoff_energy_flux(
            p.lake.runoff / max(f, sqrt(eps(FT))),
            p.lake.T,
            p.lake.q_l,
            earth_param_set,
        )

    # Surface energy flux
    # Precipitation enthalpy: rain arrives as liquid water at air temperature
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    @. p.lake.surface_energy_flux =
        f * (
            p.lake.R_n +
            p.lake.turbulent_fluxes.lhf +
            p.lake.turbulent_fluxes.shf +
            p.drivers.P_liq * _ρ_l * _cp_l * (p.drivers.T - _T_ref)
        )
    return nothing
end

# ─── Tendency ─────────────────────────────────────────────────────────────────

function ClimaLand.make_compute_exp_tendency(
    model::SlabLakeModel{FT},
) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        @. dY.lake.U =
            -p.lake.surface_energy_flux - p.lake.runoff_energy_flux +
            p.lake.sediment_heat_flux
    end
    return compute_exp_tendency!
end

# ─── Conservation diagnostics ─────────────────────────────────────────────────

function ClimaLand.total_energy_per_area!(sfc_field, ::SlabLakeModel, Y, p, t)
    @. sfc_field = Y.lake.U * p.lake.lake_fraction
end

function ClimaLand.total_liq_water_vol_per_area!(
    sfc_field,
    model::SlabLakeModel{FT},
    Y,
    p,
    t,
) where {FT}
    @. sfc_field = model.parameters.depth * p.lake.lake_fraction
end

# ─── Drivers ──────────────────────────────────────────────────────────────────

function ClimaLand.get_drivers(model::SlabLakeModel)
    bc = model.boundary_conditions
    if bc isa
       AtmosDrivenLakeBC{<:PrescribedAtmosphere, <:PrescribedRadiativeFluxes}
        return (bc.atmos, bc.radiation)
    else
        return ()
    end
end

# ─── Domain accessor ──────────────────────────────────────────────────────────

# The lake is a surface-only model. `get_domain` is used internally.
ClimaLand.get_domain(model::SlabLakeModel) = model.domain

# ─── Inland Water Mask ────────────────────────────────────────────────────────

"""
    inland_water_mask(
        surface_space;
        filepath = nothing,
        varname = "landseamask",
        threshold = 80.0,
        landsea_mask = nothing,
        latitude_bounds = (-60.0, 60.0),
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
       interpolation_method = Interpolations.Constant()
    )

Reads a water-fraction dataset from `filepath` (defaults to the IMERG
land-sea mask artifact), regrids to the `surface_space`, and identifies
inland water points.

A point is classified as inland water if its water fraction is above
`threshold` and its latitude is within `latitude_bounds`. Points
outside the latitude bounds (e.g. polar ice shelves) are never
classified as inland water, since the slab lake treatment
is not appropriate for ice sheets.

If a `landsea_mask` is provided (binary field: 1 = land,
0 = ocean), only points that are classified as land in the `landsea_mask`
but have a water fraction above `threshold` are marked as inland water.
This ensures ocean points are not double-counted.

The IMERG land-sea mask (`landseamask` variable) encodes water fraction:
0 = pure land, 100 = pure water/ocean. A `threshold` of 80 means any
grid point with > 80% water fraction that is inside the simulation's
land domain and within the latitude bounds is treated as inland water.

Returns a binary ClimaCore Field: 1 = inland water, 0 = land/ocean.
Returns `nothing` for Point/Column domains (no horizontal extent).
"""
function inland_water_mask(
    surface_space;
    filepath = nothing,
    varname = "landseamask",
    threshold = 80.0,
    landsea_mask = nothing,
    latitude_bounds = (-60.0, 60.0),
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Constant(),
)
    if isnothing(filepath)
        filepath = ClimaLand.Artifacts.imerg_landsea_mask_path()
    end
    water_frac = SpaceVaryingInput(
        filepath,
        varname,
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    # Points with high water fraction are water-like (inland lakes, rivers)
    is_high_water_frac = _apply_threshold_above.(water_frac, threshold)

    # Exclude polar regions where the slab lake treatment is inappropriate
    # (e.g. Antarctic ice shelves have high IMERG water fraction but are not lakes)
    lat_field = ClimaCore.Fields.coordinate_field(surface_space).lat
    lat_min, lat_max = latitude_bounds
    in_bounds = _within_latitude_bounds.(lat_field, lat_min, lat_max)
    is_high_water_frac = is_high_water_frac .* in_bounds

    if isnothing(landsea_mask)
        return is_high_water_frac
    else
        # Only mark as inland water if the simulation treats it as land
        # (landsea_mask == 1) but the water fraction data says it's watery
        return is_high_water_frac .* landsea_mask
    end
end

_apply_threshold_above(field, value) =
    field > value ? eltype(field)(1) : eltype(field)(0)

_within_latitude_bounds(lat, lat_min, lat_max) =
    (lat >= lat_min && lat <= lat_max) ? typeof(lat)(1) : typeof(lat)(0)

# Points and Columns do not have a horizontal dim, so a horizontal mask cannot be applied
function inland_water_mask(domain::Union{Point, Column}; kwargs...)
    error(
        "inland_water_mask is not supported for Point or Column domains (no horizontal extent)",
    )
end

function inland_water_mask(
    domain::Union{SphericalShell, SphericalSurface, HybridBox, Plane};
    filepath = nothing,
    kwargs...,
)
    # HybridBox and Plane domains might not have longlat, which is needed for the mask
    if hasproperty(domain, :longlat) && isnothing(domain.longlat)
        error(
            "inland_water_mask requires a domain with longlat coordinates to regrid the NetCDF mask",
        )
    end
    return inland_water_mask(domain.space.surface; filepath, kwargs...)
end

end # module InlandWater

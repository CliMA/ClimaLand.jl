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
    surface_albedo,
    surface_emissivity,
    surface_displacement_height,
    surface_roughness_model,
    get_∂T_sfc∂T_function,
    get_∂q_sfc∂T_function,
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
    CoupledAtmosphere,
    CoupledRadiativeFluxes,
    AbstractBC
import ClimaLand.Domains:
    Point, Column, SphericalShell, SphericalSurface, HybridBox, Plane
using Interpolations
using Dates
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.DataHandling

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
struct SlabLakeParameters{FT <: AbstractFloat, PSE}
    "Mixed-layer depth of the slab lake (m)"
    depth::FT
    "Lake–sediment conductance (W m⁻² K⁻¹)"
    conductance::FT
    "Open-water albedo"
    liquid_albedo::FT
    "Lake-ice albedo"
    ice_albedo::FT
    "Lake emissivity"
    emissivity::FT
    "Roughness length for momentum over lake water/ice (m)"
    z_0m::FT
    "Roughness length for scalars over lake water/ice (m)"
    z_0b::FT
    "Physical constants and clima-wide parameters"
    earth_param_set::PSE
end

Base.broadcastable(ps::SlabLakeParameters) = tuple(ps)

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
    return SlabLakeParameters{FT, typeof(earth_param_set)}(
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
struct SlabLakeModel{FT, PS, D, BC, M} <: AbstractInlandWaterModel{FT}
    "Lake parameters"
    parameters::PS
    "The lake domain (surface only)"
    domain::D
    "Boundary conditions"
    boundary_conditions::BC
    "Inland water mask on the surface space (fraction 0–1)"
    inland_water_mask::M
end

"""
    SlabLakeModel{FT}(;
        parameters,
        domain,
        boundary_conditions,
        inland_water_mask,
    )

Construct a `SlabLakeModel`.
"""
function SlabLakeModel{FT}(;
    parameters::PS,
    domain::D,
    boundary_conditions::BC,
    inland_water_mask::M,
) where {FT, PS, D, BC, M}
    return SlabLakeModel{FT, PS, D, BC, M}(
        parameters,
        domain,
        boundary_conditions,
        inland_water_mask,
    )
end

"""
    SlabLakeModel(
        FT,
        domain,
        forcing,
        toml_dict::CP.ParamDict;
        inland_water_mask = inland_water_mask(domain.space.surface),
        prognostic_land_components = (:lake,),
        parameters = SlabLakeParameters(toml_dict),
    )

Convenience constructor for `SlabLakeModel` that builds the model from a
TOML parameter dictionary, following the same pattern as for other
component models. The default inland water mask is used.
"""
function SlabLakeModel(
    FT,
    domain,
    forcing,
    toml_dict::CP.ParamDict;
    inland_water_mask = inland_water_mask(domain.space.surface),
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
    }(
        parameters,
        domain,
        boundary_conditions,
        inland_water_mask,
    )
end

ClimaLand.name(::SlabLakeModel) = :lake
ClimaLand.get_earth_param_set(model::SlabLakeModel) =
    model.parameters.earth_param_set

# ─── Interface: prognostic / auxiliary variable declarations ──────────────────

ClimaLand.prognostic_vars(::SlabLakeModel) = (:U,) # per lake area
ClimaLand.prognostic_types(::SlabLakeModel{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(::SlabLakeModel) = (:surface,)

function ClimaLand.auxiliary_vars(::SlabLakeModel)
    return (
        :T,
        :q_l,
        :sediment_heat_flux,
        :runoff_energy_flux,
        :runoff,
        :albedo,
        :turbulent_fluxes,
        :R_n,
    )
end
get_turb_fluxes_type(FT, ::AtmosDrivenLakeBC) = NamedTuple{
    (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T),
    Tuple{FT, FT, FT, FT, FT},
}
get_turb_fluxes_type(
    FT,
    ::AtmosDrivenLakeBC{<:CoupledAtmosphere, <:CoupledRadiativeFluxes},
) = NamedTuple{
    (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T, :ρτxz, :ρτyz, :buoyancy_flux),
    Tuple{FT, FT, FT, FT, FT, FT, FT, FT},
}

function ClimaLand.auxiliary_types(m::SlabLakeModel{FT}) where {FT}
    turb_fluxes_type = get_turb_fluxes_type(FT, m.boundary_conditions)
    return (
        FT, # T
        FT, # q_l
        FT, # sediment_heat_flux
        FT, # runoff_energy_flux
        FT, # runoff
        FT, # albedo
        turb_fluxes_type,
        FT, # R_n
    )
end

function ClimaLand.auxiliary_domain_names(::SlabLakeModel)
    return (
        :surface, # T
        :surface, # q_l
        :surface, # sediment_heat_flux
        :surface, # runoff_energy_flux
        :surface, # runoff
        :surface, # albedo
        :surface, # turbulent_fluxes
        :surface, # R_n
    )
end

# ─── Surface interface methods ────────────────────────────────────────────────

ClimaLand.component_temperature(::SlabLakeModel, Y, p) = p.lake.T
ClimaLand.surface_albedo(::SlabLakeModel, Y, p) = p.lake.albedo
ClimaLand.surface_emissivity(m::SlabLakeModel, Y, p) = m.parameters.emissivity

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
        # Thermodynamic state from prognostic U
        @. p.lake.q_l = lake_liquid_fraction(Y.lake.U, params)
        @. p.lake.T = lake_temperature(Y.lake.U, p.lake.q_l, params)
        @. p.lake.albedo = lake_surface_albedo(p.lake.q_l, params)
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
    lake_boundary_fluxes!(
        bc::AtmosDrivenLakeBC,
        prognostic_land_components::Val{(:lake,)},
        model::SlabLakeModel{FT},
        Y,
        p,
        t,
    ) where {FT}

Compute the lake surface boundary fluxes, with the option to dispatch off
of the value of the `prognostic_land_components`. 

For integrated models, the only change is that a sediment heat flux
 is added. It is assumed to be zero here.
Future interactions with the snow component will also be handled there.
"""
function lake_boundary_fluxes!(
    bc::AtmosDrivenLakeBC,
    prognostic_land_components::Val{(:lake,)},
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

    # Runoff = negative of net water flux of precip + evap
    # (excess leaves as runoff)
    @. p.lake.runoff = -(
        p.drivers.P_liq + p.drivers.P_snow + p.lake.turbulent_fluxes.vapor_flux
    )

    @. p.lake.runoff_energy_flux = lake_runoff_energy_flux(
        p.lake.runoff,
        p.lake.T,
        p.lake.q_l,
        earth_param_set,
    )

    # In the future, if we wish to model this for standalone runs, we should add a PrescribedSoil
    # to the AtmosDrivenLakeBC struct and extract what is needed (κ_soil, Δz_soil, T_soil) from there.
    @. p.lake.sediment_heat_flux = 0
    return nothing
end

# ─── Tendency ─────────────────────────────────────────────────────────────────

function ClimaLand.make_compute_exp_tendency(
    model::SlabLakeModel{FT},
) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        # Precipitation enthalpy: rain arrives as liquid water at air temperature
        earth_param_set = model.parameters.earth_param_set
        _ρ_l = LP.ρ_cloud_liq(earth_param_set)
        _cp_l = LP.cp_l(earth_param_set)
        _cp_i = LP.cp_i(earth_param_set)
        _T_ref = LP.T_0(earth_param_set)
        _LH_f0 = LP.LH_f0(earth_param_set)

        @. dY.lake.U = -(
            p.lake.R_n +
            p.lake.turbulent_fluxes.lhf +
            p.lake.turbulent_fluxes.shf +
            p.drivers.P_liq * _ρ_l * _cp_l * (p.drivers.T - _T_ref) +
            p.drivers.P_snow *
            _ρ_l *
            (_cp_i * (p.drivers.T - _T_ref) - _LH_f0) +
            p.lake.runoff_energy_flux - p.lake.sediment_heat_flux
        )
    end
    return compute_exp_tendency!
end

# ─── Conservation diagnostics ─────────────────────────────────────────────────

function ClimaLand.total_energy_per_area!(
    sfc_field,
    model::SlabLakeModel,
    Y,
    p,
    t,
)
    mask = model.inland_water_mask
    @. sfc_field = Y.lake.U * mask
end

function ClimaLand.total_liq_water_vol_per_area!(
    sfc_field,
    model::SlabLakeModel{FT},
    Y,
    p,
    t,
) where {FT}
    mask = model.inland_water_mask
    @. sfc_field = model.parameters.depth * mask
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
        filepath = ClimaLand.Artifacts.imerg_landsea_mask_path(),
        varname = "landseamask",
        threshold = 80.0,
        landsea_mask = nothing,
        latitude_bounds = (-60.0, 60.0),
        use_era5_lakes = true,
        era5_lake_cover_threshold = 0.5,
        era5_lake_depth_threshold = 5.0,
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = ...,
        interpolation_method = Interpolations.Constant(),
    )

Builds the inland water mask by combining two data sources:

1. **IMERG** water-fraction mask (default): points with water fraction
   above `threshold` are classified as inland water.
2. **ERA5 lakes** (when `use_era5_lakes = true`, the default): points
   where ERA5 lake cover `cl > era5_lake_cover_threshold` (default 0.5)
   **and** ERA5 lake depth `dl > era5_lake_depth_threshold` (default 5 m)
   are also classified as inland water, even if IMERG does not flag them.

The final mask is the **union** of both sources, filtered by
`latitude_bounds` (to exclude polar ice shelves) and optionally by
`landsea_mask` (to exclude ocean points).

If a `landsea_mask` is provided (binary field: 1 = land,
0 = ocean), only points that are classified as land in the `landsea_mask`
are eligible for inland water classification.

Returns a binary ClimaCore Field: 1 = inland water, 0 = land/ocean.
"""
function inland_water_mask(
    surface_space;
    filepath = ClimaLand.Artifacts.imerg_landsea_mask_path(),
    varname = "landseamask",
    threshold = 80.0,
    landsea_mask = nothing,
    latitude_bounds = (-60.0, 60.0),
    use_era5_lakes = true,
    era5_lake_cover_threshold = 0.5,
    era5_lake_depth_threshold = 5.0,
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Constant(),
)
    regridder_kwargs = (; extrapolation_bc, interpolation_method)

    # --- IMERG-based mask ---
    water_frac = SpaceVaryingInput(
        filepath,
        varname,
        surface_space;
        regridder_type,
        regridder_kwargs,
    )
    is_water = _apply_threshold_above.(water_frac, threshold)

    # Exclude polar regions (e.g. Antarctic ice shelves)
    lat_field = ClimaCore.Fields.coordinate_field(surface_space).lat
    lat_min, lat_max = latitude_bounds
    in_bounds = _within_latitude_bounds.(lat_field, lat_min, lat_max)
    is_water = is_water .* in_bounds

    # --- ERA5 lake cover + depth mask ---
    if use_era5_lakes
        # Read ERA5 lake cover
        dh_cover = DataHandling.DataHandler(
            ClimaLand.Artifacts.era5_lake_cover_path(),
            "cl",
            surface_space;
            regridder_type,
            regridder_kwargs,
        )
        dates_cover = DataHandling.available_dates(dh_cover)
        lake_cover =
            DataHandling.regridded_snapshot(dh_cover, first(dates_cover))
        close(dh_cover)

        # Read ERA5 lake depth
        dh_depth = DataHandling.DataHandler(
            ClimaLand.Artifacts.era5_lake_depth_path(),
            "dl",
            surface_space;
            regridder_type,
            regridder_kwargs,
        )
        dates_depth = DataHandling.available_dates(dh_depth)
        lake_depth =
            DataHandling.regridded_snapshot(dh_depth, first(dates_depth))
        close(dh_depth)

        FT = eltype(lake_cover)
        cl_thresh = FT(era5_lake_cover_threshold)
        dl_thresh = FT(era5_lake_depth_threshold)
        era5_is_lake = @. ifelse(
            lake_cover > cl_thresh && lake_depth > dl_thresh,
            FT(1),
            FT(0),
        )
        era5_is_lake = era5_is_lake .* in_bounds
        # Union: pixel is inland water if IMERG OR ERA5 flags it
        is_water =
            @. ifelse(is_water > FT(0) || era5_is_lake > FT(0), FT(1), FT(0))
    end

    if isnothing(landsea_mask)
        return is_water
    else
        return is_water .* landsea_mask
    end
end

_apply_threshold_above(field, value) =
    field > value ? eltype(field)(1) : eltype(field)(0)

_within_latitude_bounds(lat, lat_min, lat_max) =
    (lat >= lat_min && lat <= lat_max) ? typeof(lat)(1) : typeof(lat)(0)


# ─── Parameterizations ───────────────────────────────────────────────────────
include("./lake_parameterizations.jl")

end # module InlandWater

"""
Forward-model engine for the DK-Sor CalLMIP pipeline.

Builds the integrated ClimaLand `LandModel` (canopy + snow + soil + soil-CO2) at the
DK-Sor site and runs it for one or more years. The main entry point `run_prior_year`
returns either monthly (NEE, LHF, SHF) fluxes for the calibration objective or, with
`callmip=true`, the full set of daily CalLMIP diagnostics; it also supports multi-year
continuous runs, spin-up, and parameter/soil/radiation overrides. This file is included
by the calibration, posterior, and output scripts so they share one model definition.
Helpers to load FLUXNET observations and compute skill metrics are included.
"""

import ClimaComms
ClimaComms.@import_required_backends

using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.TimeManager: date
import EnsembleKalmanProcesses as EKP     # needed before loading observations.jld2
using NCDatasets
using Dates
using Statistics
using LinearAlgebra
using Printf
using JLD2
using CairoMakie

# ── Configuration ──────────────────────────────────────────────────────────────
const FT          = Float64
const DT          = Float64(900)
const SITE_ID_VAL = :DK_Sor
const CLIMALAND_DIR = pkgdir(ClimaLand)
# Forcing + flux observations come from the CalLMIP ClimaArtifacts (see src/Artifacts.jl):
#   callmip_phase1_forcing → per-site *_Met.nc  (PLUMBER2/FLUXNET2015, from ME-org)
#   callmip_phase1         → per-site *_Flux.nc (callmip-org/Phase1; "1a" = DK-Sor test)
const MET_NC_PATH   = ClimaLand.Artifacts.callmip_phase1_forcing_path("DK-Sor")
const FLUX_NC_PATH  = ClimaLand.Artifacts.callmip_phase1_flux_path("DK-Sor"; phase = "1a")
const OUTDIR        = CLIMALAND_DIR

# Years to check. Override with CALLMIP_CHECK_YEARS="1997,1998,..." to test the
# full span (e.g. the prior-crash years 2003 & 2011 that calibration flagged).
const CHECK_YEARS = haskey(ENV, "CALLMIP_CHECK_YEARS") ?
    parse.(Int, split(ENV["CALLMIP_CHECK_YEARS"], ",")) : [2005, 2008, 2010]

# CalLMIP daily diagnostic variables (run_prior_year(...; callmip=true)).
# soilrn/soillhf/soilshf (canonical soil-surface fluxes, exposed for LandModel in
# src/diagnostics) → hfg = −soilrn−soillhf−soilshf, evspsblsoi = soillhf/Lv.
const CALLMIP_SURFACE_VARS = ["nee","lhf","shf","gpp","er","hr","trans","ct","lai","cveg",
                              "soilrn","soillhf","soilshf"]
const CALLMIP_COLUMN_VARS  = ["swc","tsoil","soc"]

const _SITE_LOC    = FluxnetSimulations.get_location(FT, Val(SITE_ID_VAL))
const _SITE_HEIGHT = FluxnetSimulations.get_fluxtower_height(FT, Val(SITE_ID_VAL))
const _SITE_PARAMS = FluxnetSimulations.get_parameters(FT, Val(SITE_ID_VAL))

const lat         = _SITE_LOC.lat
const long        = _SITE_LOC.long
const time_offset = _SITE_LOC.time_offset
const atmos_h     = _SITE_HEIGHT.atmos_h

"""
Run one year, return (nee_monthly, lhf_monthly, shf_monthly).
`param_overrides` is an optional Dict(name => value); empty ⇒ default (prior)
parameters. Pass the posterior mean to evaluate the calibrated model.
"""
function run_prior_year(year; param_overrides::Dict = Dict{String,Float64}(),
                        dt = DT, callmip::Bool = false, with_co2::Bool = true,
                        return_sim::Bool = false, domain_override = nothing,
                        free_drainage::Bool = false, soil_overrides = (;),
                        rad_overrides = (;),
                        T_init_override = nothing, met_path_override = nothing,
                        spinup_days::Int = 0, stop_year::Int = year,
                        init_state = nothing, return_state::Bool = false,
                        multi_col::Bool = false,
                        multi_col_sites::Vector{String} =
                            [replace(String(SITE_ID_VAL), "_" => "-")],
                        multi_col_time_offsets = nothing,
                        multi_col_atmos_h = nothing)
    met_path = isnothing(met_path_override) ? MET_NC_PATH : met_path_override
    # spinup_days>0 starts the run earlier (e.g. a 60-day spin-up).
    # stop_year>year runs CONTINUOUSLY through (stop_year, 12, 31) in a single
    # integration (preserves inter-annual soil/carbon memory); default = one year.
    # NB: the :daily reduction can't emit the final day's mean (its window closes at
    # stop_date 00:00, beyond the forcing extent), so the last calendar day is fill.
    start_date = DateTime(year, 1, 1) - Day(spinup_days)
    stop_date  = DateTime(stop_year + 1, 1, 1)

    # Non-TOML params (site/radiation, e.g. α_NIR_leaf) are pulled out of
    # param_overrides and applied to rad_overrides below, not written to the TOML.
    param_overrides = Dict(param_overrides)
    for nm in ("α_NIR_leaf", "α_PAR_leaf", "ϵ_canopy")
        if haskey(param_overrides, nm)
            rad_overrides = merge((; Symbol(nm) => pop!(param_overrides, nm)), rad_overrides)
        end
    end

    if isempty(param_overrides)
        local_toml = LP.create_toml_dict(FT)          # default (prior) parameters
    else
        tmpf = tempname() * ".toml"
        open(tmpf, "w") do io
            for (nm, v) in param_overrides
                println(io, "[\"$nm\"]"); println(io, "value = $(Float64(v))")
                println(io, "type  = \"float\""); println(io)
            end
        end
        local_toml = LP.create_toml_dict(FT; override_files = [tmpf])
        rm(tmpf; force = true)
    end

    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(SITE_ID_VAL))
    # Optional domain override (e.g. test the global model's 50 m / (10,0.5) /
    # 15-layer column, whose deep base stays below the seasonal freezing depth).
    if !isnothing(domain_override)
        haskey(domain_override, :zmin)     && (zmin = domain_override.zmin)
        haskey(domain_override, :nelements)&& (nelements = domain_override.nelements)
        haskey(domain_override, :dz_tuple) && (dz_tuple = FT.(domain_override.dz_tuple))
    end

    # Multi-column run: gather per-site forcing paths, coordinates, and UTC offsets.
    # Any `callmip_phase1_forcing` site works; coordinates come from `get_location`
    # when a site file exists, else from the met file's `latitude`/`longitude`. The
    # shared vertical mesh and parameters are the reference site's (SITE_ID_VAL).
    if multi_col
        multi_col_met_paths =
            [ClimaLand.Artifacts.callmip_phase1_forcing_path(s) for s in multi_col_sites]
        _has_loc(sym) =
            hasmethod(FluxnetSimulations.get_location, Tuple{typeof(FT), Val{sym}})
        locs = map(enumerate(multi_col_sites)) do (i, s)
            sym = Symbol(replace(s, "-" => "_"))
            if _has_loc(sym)
                loc = FluxnetSimulations.get_location(FT, Val(sym))
                (; long = FT(loc.long), lat = FT(loc.lat), offset = loc.time_offset)
            else
                lo, la = NCDataset(multi_col_met_paths[i], "r") do ds
                    (Float64(ds["longitude"][1, 1]), Float64(ds["latitude"][1, 1]))
                end
                (; long = FT(lo), lat = FT(la), offset = nothing)
            end
        end
        multi_col_longlats = [(l.long, l.lat) for l in locs]
        multi_col_offsets = if !isnothing(multi_col_time_offsets)
            collect(multi_col_time_offsets)
        else
            map(zip(multi_col_sites, locs)) do (s, l)
                isnothing(l.offset) &&
                    error("No get_location for multi_col site $s; pass \
                           multi_col_time_offsets to specify its UTC offset.")
                l.offset
            end
        end
        # Per-column atmospheric reference height (tower measurement height): each
        # site's get_fluxtower_height when defined, else a required override vector.
        _has_height(sym) = hasmethod(
            FluxnetSimulations.get_fluxtower_height,
            Tuple{typeof(FT), Val{sym}},
        )
        multi_col_heights = if !isnothing(multi_col_atmos_h)
            collect(multi_col_atmos_h)
        else
            map(multi_col_sites) do s
                sym = Symbol(replace(s, "-" => "_"))
                _has_height(sym) || error(
                    "No get_fluxtower_height for multi_col site $s; pass \
                     multi_col_atmos_h to specify its atmos_h.",
                )
                FT(FluxnetSimulations.get_fluxtower_height(FT, Val(sym)).atmos_h)
            end
        end
        # Per-column site geometry (SAI, rooting_depth, canopy height): each site's
        # get_parameters when defined, else the reference site's value. These are
        # site-intrinsic physical attributes, not shared calibratable parameters.
        multi_col_geom = map(multi_col_sites) do s
            sym = Symbol(replace(s, "-" => "_"))
            if hasmethod(
                FluxnetSimulations.get_parameters,
                Tuple{typeof(FT), Val{sym}},
            )
                p = FluxnetSimulations.get_parameters(FT, Val(sym))
                (;
                    SAI = FT(p.SAI),
                    rooting_depth = FT(p.rooting_depth),
                    h_canopy = FT(p.h_canopy),
                )
            else
                (;
                    SAI = FT(_SITE_PARAMS.SAI),
                    rooting_depth = FT(_SITE_PARAMS.rooting_depth),
                    h_canopy = FT(_SITE_PARAMS.h_canopy),
                )
            end
        end
    end
    # multi_col=true runs on a MultiColumnFiniteDifferenceSpace (one column per site);
    # it returns a Column-typed domain so everything downstream is unchanged. See
    # ColumnEnsemble in src/shared_utilities/Domains.jl.
    land_domain = if multi_col
        ColumnEnsemble(;
            zlim      = (FT(zmin), FT(zmax)),
            nelements, dz_tuple,
            longlat   = multi_col_longlats,
        )
    else
        Column(;
            zlim      = (FT(zmin), FT(zmax)),
            nelements, dz_tuple,
            longlat   = (long, lat),
        )
    end
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    # multi_col=true remaps each site's forcing onto its column of the N-point surface
    # via a MultiColumnDataHandler, with per-column UTC offsets and a per-column
    # atmos_h (a surface Field of tower heights scattered in multi_col_sites order).
    forcing = if multi_col
        atmos_h_field = ClimaCore.Fields.array2field(
            FT.(multi_col_heights),
            land_domain.space.surface,
        )
        FluxnetSimulations.prescribed_forcing_netcdf_multicol(
            multi_col_met_paths, land_domain.space.surface, atmos_h_field,
            start_date, local_toml, FT;
            hour_offset_from_UTC = multi_col_offsets,
        )
    else
        FluxnetSimulations.prescribed_forcing_netcdf(
            met_path, lat, long, time_offset, atmos_h,
            start_date, local_toml, FT,
        )
    end

    # LAI: single-col reads one met file into a scalar-in-space TimeVaryingInput;
    # multi_col builds a per-column (space+time) LAI from each site's LAI_alternative
    # via the same MultiColumnDataHandler machinery as the forcing, plus a per-site
    # `maxLAI` vector (same max-of-valid reduction) used for the per-column RAI below.
    if multi_col
        LAI = FluxnetSimulations.prescribed_LAI_netcdf_multicol(
            multi_col_met_paths, land_domain.space.surface, start_date, FT;
            hour_offset_from_UTC = multi_col_offsets,
        )
        maxLAI = map(multi_col_met_paths) do p
            NCDataset(p, "r") do ds
                lai_data = Float64.(coalesce.(ds["LAI_alternative"][1, 1, :], NaN))
                maximum(lai_data[.!isnan.(lai_data)])
            end
        end
    else
        LAI, maxLAI = NCDataset(met_path, "r") do ds
            time_vals = ds["time"][:]
            lai_data  = Float64.(coalesce.(ds["LAI_alternative"][1, 1, :], NaN))
            lai_secs  = Float64[
                Second(t - Hour(time_offset) - start_date).value
                for t in time_vals
            ]
            valid = .!isnan.(lai_data)
            TimeVaryingInput(lai_secs[valid], lai_data[valid]), maximum(lai_data[valid])
        end
    end

    (; soil_ν, θ_r, soil_K_sat, soil_vg_n, soil_vg_α,
       ν_ss_om, ν_ss_quartz, ν_ss_gravel,
       z_0m_soil, z_0b_soil, soil_α_PAR, soil_α_NIR, soil_ϵ,
       Ω, χl, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
       ϵ_canopy, ac_canopy, Drel, g0,
       SAI, f_root_to_shoot, rooting_depth, h_canopy,
       conductivity_model, retention_model, plant_ν, plant_S_s) = _SITE_PARAMS

    # Optional soil-parameter overrides (e.g. test alternative values:
    # soil_K_sat=1e-5, ν_ss_om=0.03 — better-drained, lower-organic → higher κ).
    haskey(soil_overrides, :soil_K_sat) && (soil_K_sat = FT(soil_overrides.soil_K_sat))
    haskey(soil_overrides, :ν_ss_om)    && (ν_ss_om    = FT(soil_overrides.ν_ss_om))
    haskey(soil_overrides, :ν_ss_quartz)&& (ν_ss_quartz= FT(soil_overrides.ν_ss_quartz))
    haskey(soil_overrides, :ν_ss_gravel)&& (ν_ss_gravel= FT(soil_overrides.ν_ss_gravel))

    # Optional radiation/leaf-optics overrides (site params, not in the TOML) to
    # probe net-radiation control of SHF (α_NIR dominates SW_n absorbed energy).
    haskey(rad_overrides, :α_NIR_leaf) && (α_NIR_leaf = FT(rad_overrides.α_NIR_leaf))
    haskey(rad_overrides, :τ_NIR_leaf) && (τ_NIR_leaf = FT(rad_overrides.τ_NIR_leaf))
    haskey(rad_overrides, :α_PAR_leaf) && (α_PAR_leaf = FT(rad_overrides.α_PAR_leaf))
    haskey(rad_overrides, :ϵ_canopy)   && (ϵ_canopy   = FT(rad_overrides.ϵ_canopy))

    # RAI = maxLAI * f_root_to_shoot. Single-col: scalars. Multi-col: scatter each
    # per-column quantity into a surface Field. RAI/SAI share one concrete Field type
    # (PrescribedAreaIndices{FS} requires it). Canopy height and rooting depth are
    # site geometry and vary per column too (biomass height/rooting_depth accept a
    # Field); making height a Field also gives per-column canopy roughness.
    if multi_col
        _sfc = land_domain.space.surface
        RAI = ClimaCore.Fields.array2field(FT.(maxLAI) .* FT(f_root_to_shoot), _sfc)
        SAI_bio = ClimaCore.Fields.array2field(
            FT[g.SAI for g in multi_col_geom],
            _sfc,
        )
        rooting_depth_bio = ClimaCore.Fields.array2field(
            FT[g.rooting_depth for g in multi_col_geom],
            _sfc,
        )
        h_canopy_bio = ClimaCore.Fields.array2field(
            FT[g.h_canopy for g in multi_col_geom],
            _sfc,
        )
    else
        RAI = maxLAI * f_root_to_shoot
        SAI_bio = SAI
        rooting_depth_bio = rooting_depth
        h_canopy_bio = h_canopy
    end

    prognostic_land_components =
        with_co2 ? (:canopy, :snow, :soil, :soilco2) : (:canopy, :snow, :soil)

    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = soil_α_PAR, NIR_albedo = soil_α_NIR)
    retention_parameters = (;
        ν      = soil_ν,
        θ_r,
        K_sat  = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    )
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

    # Bottom BC: default is the convenience-ctor zero-flux (water+heat); with
    # free_drainage=true use EnergyWaterFreeDrainage (water drains out the base,
    # the physical case for the permeable Sorø subsurface).
    bottom_bc_kw = free_drainage ?
        (; bottom_bc = Soil.EnergyWaterFreeDrainage()) : (;)
    soil = Soil.EnergyHydrology{FT}(
        land_domain, forcing, local_toml;
        prognostic_land_components,
        additional_sources   = (ClimaLand.RootExtraction{FT}(),),
        albedo               = soil_albedo,
        runoff               = ClimaLand.Soil.Runoff.SurfaceRunoff(),
        retention_parameters, composition_parameters,
        S_s = FT(1e-3), z_0m = z_0m_soil, z_0b = z_0b_soil,
        bottom_bc_kw...,
    )

    soilco2 = if with_co2
        co2_drivers = Soil.Biogeochemistry.SoilDrivers(
            Soil.Biogeochemistry.PrognosticMet(soil.parameters), forcing.atmos,
        )
        Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, co2_drivers, local_toml)
    else
        nothing
    end

    radiation_parameters = (;
        Ω, G_Function = CLMGFunction(χl),
        α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
    )
    radiative_transfer = Canopy.TwoStreamModel{FT}(
        canopy_domain, local_toml; radiation_parameters, ϵ_canopy,
    )

    biomass = Canopy.PrescribedBiomassModel{FT}(;
        LAI, SAI = SAI_bio, RAI,
        rooting_depth = rooting_depth_bio, height = h_canopy_bio,
    )
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        (; atmos = forcing.atmos, radiation = forcing.radiation,
           ground = ClimaLand.PrognosticGroundConditions{FT}()),
        LAI, local_toml;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis       = Canopy.PModel{FT}(canopy_domain, local_toml),
        conductance          = Canopy.PModelConductance{FT}(local_toml; Drel),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, local_toml; soil_params = (; ν = soil_ν, θ_r)),
        hydraulics = Canopy.PlantHydraulicsModel{FT}(
            canopy_domain, local_toml;
            ν = plant_ν, S_s = plant_S_s, conductivity_model, retention_model),
        energy  = Canopy.BigLeafEnergyModel{FT}(local_toml; ac_canopy),
        biomass,
    )
    snow = Snow.SnowModel(
        FT, canopy_domain, forcing, local_toml, dt;
        prognostic_land_components,
    )
    land = LandModel{FT}(canopy, snow, soil, soilco2, nothing)

    # Default IC temperature = met Tair at the start of the year (Jan → near
    # freezing), applied uniformly to the whole column. T_init_override lets us
    # test a realistic warm deep-soil temperature (annual mean ~281 K) — the
    # real deep soil is NOT at January's surface temperature.
    T_init_K = if !isnothing(T_init_override)
        Float64(T_init_override)
    else
        NCDataset(met_path, "r") do ds
            t_dates = DateTime.(ds["time"][:])
            idx_yr  = findfirst(t -> Dates.year(t) == year, t_dates)
            isnothing(idx_yr) ? 283.15 : Float64(coalesce(ds["Tair"][1, 1, idx_yr], 283.15))
        end
    end

    # Per-column initial temperature for multi_col: each column starts at its own
    # site's Tair at the start of `year` (same read as T_init_K), scattered into a
    # surface Field (canopy energy) and a subsurface Field (soil internal energy).
    # Single-col keeps the scalar T_init_K path unchanged.
    T_init_sfc = nothing
    T_init_sub = nothing
    if multi_col
        Tinit_vec = if !isnothing(T_init_override)
            fill(FT(T_init_override), length(multi_col_met_paths))
        else
            map(multi_col_met_paths) do p
                NCDataset(p, "r") do ds
                    t_dates = DateTime.(ds["time"][:])
                    idx_yr  = findfirst(t -> Dates.year(t) == year, t_dates)
                    FT(
                        isnothing(idx_yr) ? 283.15 :
                        Float64(coalesce(ds["Tair"][1, 1, idx_yr], 283.15)),
                    )
                end
            end
        end
        _sfc = land_domain.space.surface
        _sub = land_domain.space.subsurface
        _Nv = nelements   # cell-center vertical levels of the shared mesh
        _N = length(Tinit_vec)
        T_init_sfc = ClimaCore.Fields.array2field(collect(Tinit_vec), _sfc)
        T_init_sub = ClimaCore.Fields.array2field(
            repeat(reshape(collect(Tinit_vec), 1, _N), _Nv, 1), _sub)
    end

    function set_ic!(Y, p, t, model)
        FT_l = eltype(Y.soil.ρe_int)
        if isnothing(init_state)
            # Cold start from arbitrary IC (soil near saturation, T = met Tair).
            ν_l  = FT_l(soil_ν)
            θr_l = FT_l(θ_r)
            Y.soil.ϑ_l .= θr_l + (ν_l - θr_l) * FT_l(0.95)
            Y.soil.θ_i .= FT_l(0)
            ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l, Y.soil.θ_i,
                soil.parameters.ρc_ds, soil.parameters.earth_param_set)
            Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
                Y.soil.θ_i, ρc_s, multi_col ? T_init_sub : FT_l(T_init_K),
                soil.parameters.earth_param_set)
            Y.snow.S .= FT_l(0); Y.snow.S_l .= FT_l(0); Y.snow.U .= FT_l(0)
            Y.canopy.energy.T .= multi_col ? T_init_sfc : FT_l(T_init_K)
            Y.canopy.hydraulics.ϑ_l .= model.canopy.hydraulics.parameters.ν
        else
            # Restart the WATER + ENERGY (temperature) prognostic states from a
            # spun-up state (parent() copies are robust to space-identity). Carbon
            # is NOT carried — it is prescribed below (SOC/CO2/O2 from files).
            parent(Y.soil.ϑ_l)              .= parent(init_state.soil.ϑ_l)
            parent(Y.soil.θ_i)              .= parent(init_state.soil.θ_i)
            parent(Y.soil.ρe_int)           .= parent(init_state.soil.ρe_int)
            parent(Y.snow.S)                .= parent(init_state.snow.S)
            parent(Y.snow.S_l)              .= parent(init_state.snow.S_l)
            parent(Y.snow.U)                .= parent(init_state.snow.U)
            parent(Y.canopy.energy.T)       .= parent(init_state.canopy.energy.T)
            parent(Y.canopy.hydraulics.ϑ_l) .= parent(init_state.canopy.hydraulics.ϑ_l)
        end
        # Carbon pools are PRESCRIBED (from ClimaLand defaults), never spun up.
        if with_co2
            Y.soilco2.CO2 .= FT_l(6e-5)
            Y.soilco2.O2  .= FT_l(0.21)
            SOC_top = FT_l(15.0); SOC_bot = FT_l(0.5)
            τ_soc   = FT_l(1.0 / log(SOC_top / SOC_bot))
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
            @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
        end
    end

    ow = ClimaDiagnostics.Writers.DictWriter()
    if callmip
        diags = vcat(
            ClimaLand.default_diagnostics(land, start_date, "";
                output_writer = ow, output_vars = CALLMIP_SURFACE_VARS,
                reduction_period = :daily),
            ClimaLand.default_diagnostics(land, start_date, "";
                output_writer = ow, output_vars = CALLMIP_COLUMN_VARS,
                reduction_period = :daily),
        )
    else
        diags = ClimaLand.default_diagnostics(land, start_date, "";
            output_writer = ow, output_vars = ["lhf", "shf", "nee"],
            reduction_period = :daily)
    end
    simulation = LandSimulation(
        start_date, stop_date, dt, land;
        set_ic!, updateat = Second(dt), diagnostics = diags,
    )

    # Debug hook: return the (un-solved) simulation so a driver can step it
    # manually and inspect the soil state at the divergence onset.
    return_sim && return simulation

    solve!(simulation)

    # Spin-up hook: return the final prognostic state so a driver can restart the
    # output run from a spun-up (equilibrated water+temp) state.
    return_state && return simulation._integrator.u

    if callmip
        # Daily CalLMIP output: surface vars → Vector, column vars → (n_z × n_days).
        surface_data = Dict{String, Vector{Float64}}()
        ref_dates    = Date[]
        for var in CALLMIP_SURFACE_VARS
            try
                times, data =
                    ClimaLand.Diagnostics.diagnostic_as_vectors(ow, "$(var)_1d_average")
                surface_data[var] = Float64.(data)
                # ClimaDiagnostics timestamps each :daily mean at the END of its
                # averaging window, so the value written at 00:00 of day D+1 is the
                # mean over day D. Subtract 1 day so each daily mean is labelled by
                # the calendar day it represents (verified empirically: the model
                # otherwise leads the FLUXNET obs by exactly 1 day). This also fills
                # 1997-01-01 and drops the spurious trailing 2015-01-01.
                isempty(ref_dates) && (ref_dates = Date.(date.(times)) .- Day(1))
            catch
                @warn "  [$year] surface '$var' not found"
            end
        end
        column_data = Dict{String, Matrix{Float64}}()
        for var in CALLMIP_COLUMN_VARS
            key = "$(var)_1d_average"
            if haskey(ow.dict, key)
                raw = ow.dict[key]
                tt  = sort(collect(keys(raw)))
                column_data[var] =
                    reduce(hcat, [Float64.(vec(parent(raw[t]))) for t in tt])
            else
                @warn "  [$year] column '$var' not found"
            end
        end
        z_soil = try
            Float64.(vec(parent(ClimaCore.Fields.coordinate_field(
                axes(land.soil.domain.fields.z)).z)))
        catch
            Float64[]
        end
        return surface_data, column_data, z_soil, ref_dates
    end

    function get_monthly(diag_name, scale = 1.0)
        times, data = ClimaLand.Diagnostics.diagnostic_as_vectors(ow, diag_name)
        months = Dates.month.(date.(times))
        return [mean(data[months .== m]) * scale for m in 1:12]
    end

    nee_m = get_monthly("nee_1d_average", 12.0 * 86400.0)  # mol/m²/s → gC/m²/d
    lhf_m = get_monthly("lhf_1d_average")
    shf_m = get_monthly("shf_1d_average")
    return nee_m, lhf_m, shf_m
end

"""Load monthly obs for a given year from the FLUXNET daily NC file."""
function load_obs_year(year)
    NCDataset(FLUX_NC_PATH, "r") do ds
        dates   = DateTime.(ds["time"][:])
        nee_all = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
        qle_all = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
        qh_all  = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
        for arr in (nee_all, qle_all, qh_all)
            arr[arr .>= 1.0e19] .= NaN
        end

        nee_m = Float64[]; lhf_m = Float64[]; shf_m = Float64[]
        for mon in 1:12
            mask = (Dates.year.(dates) .== year) .& (Dates.month.(dates) .== mon)
            nee_v = nee_all[mask]; qle_v = qle_all[mask]; qh_v = qh_all[mask]
            nee_f = filter(isfinite, nee_v); push!(nee_m, isempty(nee_f) ? NaN : mean(nee_f))
            lhf_f = filter(isfinite, qle_v); push!(lhf_m, isempty(lhf_f) ? NaN : mean(lhf_f))
            shf_f = filter(isfinite, qh_v);  push!(shf_m, isempty(shf_f) ? NaN : mean(shf_f))
        end
        return nee_m, lhf_m, shf_m
    end
end

function rmse_r2(mod, obs)
    mask = isfinite.(mod) .& isfinite.(obs)
    sum(mask) < 2 && return NaN, NaN
    m = mod[mask]; o = obs[mask]
    rmse = sqrt(mean((m .- o).^2))
    r2   = 1.0 - sum((m .- o).^2) / sum((o .- mean(o)).^2)
    return rmse, r2
end

# ── Run simulations ────────────────────────────────────────────────────────────
# Only when run directly (not when `include`d by compute_prior_post_stats.jl,
# which reuses run_prior_year / load_obs_year / rmse_r2).

# results_mod = Dict{Int, NamedTuple}()
# results_obs = Dict{Int, NamedTuple}()

# ok_years = Int[]
# for yr in CHECK_YEARS
#     @info "Running year $yr with default parameters..."
#     try
#         nee_m, lhf_m, shf_m = run_prior_year(yr; multi_col = false)
#         # Sanity: physical fluxes, no NaN/blowup.
#         finite = all(isfinite, nee_m) && all(isfinite, lhf_m) && all(isfinite, shf_m)
#         @info @sprintf("  %s year %d  NEE[%.2f,%.2f] LHF[%.1f,%.1f] SHF[%.1f,%.1f]",
#             finite ? "OK  " : "BAD ", yr,
#             minimum(nee_m), maximum(nee_m),
#             minimum(lhf_m), maximum(lhf_m),
#             minimum(shf_m), maximum(shf_m))
#         if finite
#             results_mod[yr] = (; nee = nee_m, lhf = lhf_m, shf = shf_m)
#             nee_o, lhf_o, shf_o = load_obs_year(yr)
#             results_obs[yr] = (; nee = nee_o, lhf = lhf_o, shf = shf_o)
#             push!(ok_years, yr)
#         end
#     catch e
#         @warn "  FAIL year $yr: $(typeof(e)): $(sprint(showerror, e))"
#     end
# end

# @info "PRIOR FORWARD-MODEL CHECK: $(length(ok_years))/$(length(CHECK_YEARS)) years OK: $ok_years"
# failed = setdiff(CHECK_YEARS, ok_years)
# isempty(failed) || @warn "FAILED years: $failed"

# # Only plot if every requested year ran (keeps the multi-year metrics honest).
# if length(ok_years) != length(CHECK_YEARS)
#     @info "Skipping plot (not all years succeeded); foundational check is the result."
#     exit(0)
# end

# # ── Aggregate metrics across all 3 years ─────────────────────────────────────
# function aggregate_metric(var_key)
#     mod_all = vcat([results_mod[yr][var_key] for yr in CHECK_YEARS]...)
#     obs_all = vcat([results_obs[yr][var_key] for yr in CHECK_YEARS]...)
#     return rmse_r2(mod_all, obs_all)
# end

# nee_rmse, nee_r2 = aggregate_metric(:nee)
# lhf_rmse, lhf_r2 = aggregate_metric(:lhf)
# shf_rmse, shf_r2 = aggregate_metric(:shf)

# @info "Prior performance (2005/2008/2010 combined):"
# @info "  NEE: RMSE=$(round(nee_rmse; digits=2)) gC/m²/d, R²=$(round(nee_r2; digits=3))"
# @info "  LHF: RMSE=$(round(lhf_rmse; digits=2)) W/m², R²=$(round(lhf_r2; digits=3))"
# @info "  SHF: RMSE=$(round(shf_rmse; digits=2)) W/m², R²=$(round(shf_r2; digits=3))"

# # ── Plot ───────────────────────────────────────────────────────────────────────
# months = 1:12
# var_labels = ["NEE (gC m⁻² d⁻¹)", "LHF (W m⁻²)", "SHF (W m⁻²)"]
# var_keys   = [:nee, :lhf, :shf]
# metrics    = [(nee_rmse, nee_r2), (lhf_rmse, lhf_r2), (shf_rmse, shf_r2)]

# fig = Figure(size = (1000, 800))

# for (row, (vk, ylabel, (rmse, r2))) in enumerate(zip(var_keys, var_labels, metrics))
#     for (col, yr) in enumerate(CHECK_YEARS)
#         ax = Axis(fig[row, col];
#             xlabel = col == 2 ? "Month" : "",
#             ylabel = col == 1 ? ylabel : "",
#             title  = row == 1 ? string(yr) : "",
#             xticks = (1:12, string.(1:12)),
#         )

#         obs_v = results_obs[yr][vk]
#         mod_v = results_mod[yr][vk]

#         # Obs: black dots for valid months
#         valid_months = findall(isfinite, obs_v)
#         scatter!(ax, valid_months, obs_v[valid_months];
#             color = :black, markersize = 8, label = "FLUXNET")

#         # Model: colored line
#         lines!(ax, months, mod_v; color = Makie.wong_colors()[row], linewidth = 2,
#             label = "ClimaLand prior")

#         if col == 1 && row == 1
#             axislegend(ax; position = :rt, labelsize = 10)
#         end
#     end
# end

# # Overall title with metrics
# Label(fig[0, :],
#     "DK-Sor prior (default params) vs FLUXNET — 2005/2008/2010\n" *
#     "NEE: RMSE=$(round(nee_rmse; digits=2)) gC/m²/d, R²=$(round(nee_r2; digits=3))  |  " *
#     "LHF: RMSE=$(round(lhf_rmse; digits=2)) W/m², R²=$(round(lhf_r2; digits=3))  |  " *
#     "SHF: RMSE=$(round(shf_rmse; digits=2)) W/m², R²=$(round(shf_r2; digits=3))",
#     fontsize = 13, font = :bold,
# )

# save(joinpath(OUTDIR, "prior_vs_obs.png"), fig; px_per_unit = 2)
# @info "Plot saved → $(joinpath(OUTDIR, "prior_vs_obs.png"))"

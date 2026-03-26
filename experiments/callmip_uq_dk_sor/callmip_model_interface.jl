"""
ClimaCalibrate model interface for CalLMIP Phase 1 output at DK-Sor.

Identical land model setup to model_interface.jl, but requests the full
set of CalLMIP-required diagnostic variables and saves them to JLD2.

Diagnostics saved and their physical meaning:
  Surface diagnostics (scalar timeseries)
  ─────────────────────────────────────────
  nee      Net Ecosystem Exchange          mol CO₂ m⁻² s⁻¹   (→ NEE)
  lhf      Total Latent Heat Flux          W m⁻²              (→ Qle)
  shf      Total Sensible Heat Flux        W m⁻²              (→ Qh)
  gpp      Gross Primary Production        mol CO₂ m⁻² s⁻¹   (→ GPP)
  er       Total Ecosystem Respiration     mol CO₂ m⁻² s⁻¹   (→ Reco)
  trans    Canopy Transpiration            kg m⁻² s⁻¹         (→ TVeg)
  soillhf  Soil Surface Latent Heat Flux   W m⁻²              (→ ESoil via /Lv)
  soilrn   Soil Net Radiation              W m⁻²              (for Qg)
  soilshf  Soil Surface Sensible Heat Flux W m⁻²              (for Qg)
  rn       Total Net Radiation             W m⁻²              (for Qg cross-check)
  ct       Canopy Temperature              K                  (→ AvgSurfT)
  lai      Leaf Area Index                 m² m⁻²             (→ LAI)
  cveg     Vegetation Carbon               kg C m⁻²           (→ TotAbovBioMass)

  Column diagnostics (profile timeseries → integrated in write_callmip_netcdf.jl)
  ─────────────────────────────────────────────────────────────────────────────────
  swc      Soil Water Content per layer    m³ m⁻³             (→ SoilMoist)
  tsoil    Soil Temperature per layer      K                  (→ AvgSurfT fallback)
  soc      Soil Organic Carbon per layer   kg C m⁻³           (→ TotSoilCarb)

Scalar unit conversions applied at write time (write_callmip_netcdf.jl):
  Carbon fluxes: × 12e−3   mol CO₂/m²/s → kg C/m²/s  (ALMA standard)
  ESoil:         ÷ 2.5e6   W/m²         → kg/m²/s     (÷ latent heat of vap.)
  Qg:            soilrn − soillhf − soilshf             (energy balance)
  SoilMoist:     Σ(swc × dz × ρ_water)                 (full column, kg/m²)
  TotSoilCarb:   Σ(soc × dz)                            (full column, kg C/m²)
"""

import ClimaCalibrate
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date
using Insolation

import EnsembleKalmanProcesses as EKP
import JLD2
using Dates
using NCDatasets
using Statistics

# CALLMIP_VARS: all diagnostics requested for CalLMIP output.
# Surface diagnostics (→ scalar per day after daily averaging).
# NOTE: if a variable is unavailable in the model, it will be silently absent.
const CALLMIP_SURFACE_VARS = [
    "nee",      # Net Ecosystem Exchange
    "lhf",      # Total Latent Heat Flux (Qle)
    "shf",      # Total Sensible Heat Flux (Qh)
    "gpp",      # Gross Primary Production
    "er",       # Total Ecosystem Respiration (Ra + Rh)
    "trans",    # Transpiration
    "soillhf",  # Soil LHF (for ESoil = soillhf / Lv)
    "soilrn",   # Soil net radiation (for Qg)
    "soilshf",  # Soil SHF (for Qg = soilrn - soillhf - soilshf)
    "rn",       # Total net radiation (cross-check)
    "ct",       # Canopy temperature (AvgSurfT primary)
    "lai",      # Leaf Area Index
    "cveg",     # Vegetation carbon (TotAbovBioMass)
]
# Column diagnostics (→ vertical profile per day; integrated to total-column in writer).
const CALLMIP_COLUMN_VARS = [
    "swc",      # Soil water content profile  (m³/m³) → SoilMoist
    "tsoil",    # Soil temperature profile    (K)     → AvgSurfT fallback
    "soc",      # Soil organic carbon profile (kg C/m³) → TotSoilCarb
]

# ── Forward Model ──────────────────────────────────────────────────────────────

function ClimaCalibrate.forward_model(iteration, member)
    FT = Float64
    site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)
    climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))

    sim_start = DateTime(2003, 1, 1)
    sim_stop  = DateTime(2014, 1, 1)
    @info "Member $member: simulating $sim_start to $sim_stop"

    calibrate_params_path =
        ClimaCalibrate.parameter_path(OUTPUT_DIR, iteration, member)
    toml_dict =
        LP.create_toml_dict(FT; override_files = [calibrate_params_path])

    # Domain
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID_val))
    (; atmos_h) =
        FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

    land_domain = Column(;
        zlim      = (zmin, zmax),
        nelements = nelements,
        dz_tuple  = dz_tuple,
        longlat   = (long, lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    # Meteorological forcing
    met_nc_path = joinpath(climaland_dir, "DK_Sor",
                           "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
    (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
        met_nc_path, lat, long, time_offset, atmos_h, sim_start, toml_dict, FT)

    # LAI
    surface_space = canopy_domain.space.surface
    met_ds     = NCDataset(met_nc_path, "r")
    lai_data   = Float64.(coalesce.(met_ds["LAI_alternative"][1, 1, :], NaN))
    lai_times  = met_ds["time"][:]
    close(met_ds)
    lai_seconds = [Float64(Second(t - Hour(time_offset) - sim_start).value)
                   for t in lai_times]
    valid_lai = .!isnan.(lai_data)
    LAI = TimeVaryingInput(lai_seconds[valid_lai], lai_data[valid_lai])

    # Site parameters (same as model_interface.jl)
    χl            = FT(0.25)
    α_PAR_leaf    = FT(0.1)
    α_NIR_leaf    = FT(0.45)
    τ_PAR_leaf    = FT(0.05)
    τ_NIR_leaf    = FT(0.25)
    Ω             = FT(1)
    rooting_depth = FT(0.3)
    ν             = FT(0.45)
    θ_r           = FT(0.07)
    K_sat         = FT(1e-5)
    vg_n          = FT(1.6)
    vg_α          = FT(1.6)

    hydrology_cm = ClimaLand.Soil.vanGenuchten{FT}(; α = vg_α, n = vg_n)
    retention_parameters  = (; ν, hydrology_cm, θ_r, K_sat)
    composition_parameters = (; ν_ss_om = FT(0.03),
                                ν_ss_quartz = FT(0.47),
                                ν_ss_gravel = FT(0.12))

    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
    forcing_nt = (;
        atmos     = atmos,
        radiation = radiation,
        ground    = ClimaLand.PrognosticGroundConditions{FT}(),
    )

    biomass = Canopy.PrescribedBiomassModel{FT}(
        land_domain, LAI, toml_dict;
        rooting_depth, height = FT(25), SAI = FT(1.0), RAI = FT(17.5),
    )
    radiation_parameters = (;
        Ω,
        G_Function    = CLMGFunction(χl),
        α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
    )
    radiative_transfer = Canopy.TwoStreamModel{FT}(canopy_domain, toml_dict;
                                                    radiation_parameters)
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        canopy_domain, toml_dict;
        n_stem = 1, h_stem = FT(24), h_leaf = FT(1),
    )
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain, forcing_nt, LAI, toml_dict;
        prognostic_land_components,
        photosynthesis       = Canopy.PModel{FT}(canopy_domain, toml_dict),
        conductance          = Canopy.PModelConductance{FT}(toml_dict),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, toml_dict; soil_params = (; ν, θ_r),
        ),
        biomass,
        radiative_transfer,
        hydraulics,
    )
    soil = ClimaLand.Soil.EnergyHydrology{FT}(
        land_domain, (; atmos, radiation), toml_dict;
        prognostic_land_components, retention_parameters, composition_parameters,
        S_s                 = FT(1e-3),
        additional_sources  = (ClimaLand.RootExtraction{FT}(),),
    )
    land = LandModel{FT}(
        (; atmos, radiation), LAI, toml_dict, land_domain, DT;
        prognostic_land_components, canopy, soil,
    )

    # Initial conditions
    function custom_set_ic!(Y, p, t, model)
        earth_param_set = ClimaLand.get_earth_param_set(model.soil)
        evaluate!(p.drivers.T, atmos.T, t)
        (; θ_r, ν, ρc_ds) = model.soil.parameters
        @. Y.soil.ϑ_l = θ_r + (ν - θ_r) .* FT(0.95)
        Y.soil.θ_i .= FT(0.0)
        ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l, Y.soil.θ_i, ρc_ds, earth_param_set)
        Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i, ρc_s, p.drivers.T, earth_param_set)
        Y.snow.S   .= FT(0)
        Y.snow.S_l .= FT(0)
        Y.snow.U   .= FT(0)
        if model.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
            Y.canopy.energy.T .= p.drivers.T
        end
        n_stem = model.canopy.hydraulics.n_stem
        n_leaf = model.canopy.hydraulics.n_leaf
        for i in 1:(n_stem + n_leaf)
            Y.canopy.hydraulics.ϑ_l.:($i) .=
                model.canopy.hydraulics.parameters.ν
        end
        if !isnothing(model.soilco2)
            Y.soilco2.CO2 .= FT(0.000412)
            Y.soilco2.O2_f .= FT(0.21)
            SOC_top  = FT(15.0)
            SOC_bot  = FT(0.5)
            τ_soc    = FT(1.0 / log(SOC_top / SOC_bot))
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
            @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
        end
    end

    # Diagnostics: surface (scalar) + column (profile) variables
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    # Filter to available vars only — unknown names will cause an error, so we
    # verify against get_possible_diagnostics before requesting.
    possible = ClimaLand.Diagnostics.get_possible_diagnostics(land)
    requested_surface = filter(v -> v in possible, CALLMIP_SURFACE_VARS)
    requested_column  = filter(v -> v in possible, CALLMIP_COLUMN_VARS)
    requested_all = vcat(requested_surface, requested_column)

    skipped = setdiff(vcat(CALLMIP_SURFACE_VARS, CALLMIP_COLUMN_VARS), requested_all)
    isempty(skipped) ||
        @warn "Diagnostics not available in this model build and will be skipped: $skipped"

    diags = ClimaLand.default_diagnostics(
        land, sim_start;
        output_writer    = output_writer,
        output_vars      = requested_all,
        reduction_period = :daily,
    )

    simulation = LandSimulation(
        sim_start, sim_stop, DT, land;
        set_ic!     = custom_set_ic!,
        updateat    = Second(DT),
        diagnostics = diags,
    )
    solve!(simulation)

    # Save diagnostics
    member_path = ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, iteration, member)
    save_callmip_diagnostics(simulation, member_path, requested_surface,
                             requested_column, land)
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Save function
# ─────────────────────────────────────────────────────────────────────────────

"""
    save_callmip_diagnostics(simulation, member_path, surface_vars, column_vars, land)

Extract daily diagnostics from the DictWriter and save to JLD2.
Surface diagnostic data is stored as scalar timeseries.
Column (multi-layer) diagnostic data is stored as (n_z × n_days) matrices,
together with the soil z-coordinate vector (m) so that `write_callmip_netcdf.jl`
can integrate columns to produce total-column quantities.
"""
function save_callmip_diagnostics(simulation, member_path, surface_vars,
                                   column_vars, land)
    # ── Extract scalar (surface) diagnostics ─────────────────────────────────
    surface_data = Dict{String, Any}()
    ref_dates    = nothing
    for var in surface_vars
        try
            dates, vals = extract_daily_diag(simulation, "$(var)_1d_average")
            surface_data[var] = Float64.(vals)
            ref_dates = dates
        catch e
            @warn "Could not extract surface diagnostic '$var'" exception = e
            surface_data[var] = Float64[]
        end
    end

    # ── Extract column (multi-layer) diagnostics ──────────────────────────────
    # For column fields the DictWriter stores a Field (ClimaCore) at each
    # output step.  We unwrap using parent() to recover the raw Float64 array.
    column_data = Dict{String, Any}()
    for var in column_vars
        try
            dates, vals = extract_daily_column_diag(simulation,
                                                     "$(var)_1d_average")
            column_data[var] = vals   # Matrix{Float64}: [n_z × n_days]
        catch e
            @warn "Could not extract column diagnostic '$var'" exception = e
            column_data[var] = Matrix{Float64}(undef, 0, 0)
        end
    end

    # ── Extract soil z-coordinates for column integration ─────────────────────
    z_soil = try
        Float64.(vec(parent(
            ClimaCore.Fields.coordinate_field(
                axes(simulation.model.soil.domain.fields.z)
            ).z
        )))
    catch
        Float64[]
    end

    JLD2.jldsave(
        joinpath(member_path, "callmip_diagnostics.jld2");
        dates        = ref_dates,
        surface_data = surface_data,
        column_data  = column_data,
        z_soil       = z_soil,
    )
end

"""
    extract_daily_diag(simulation, diag_name)

Extract a scalar daily diagnostic from the DictWriter.
Returns (dates::Vector{Date}, values::Vector{Float64}).
"""
function extract_daily_diag(simulation, diag_name)
    writer = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            writer = d.output_writer
            break
        end
    end
    isnothing(writer) && error("Diagnostic '$(diag_name)' not found.")
    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, diag_name)
    model_dates = Date.(times isa Vector{DateTime} ? times : date.(times))
    return model_dates, Float64.(data)
end

"""
    extract_daily_column_diag(simulation, diag_name)

Extract a multi-layer daily diagnostic from the DictWriter.
Returns (dates::Vector{Date}, data::Matrix{Float64}) where
`data` has shape [n_z × n_days] with z ordered from deepest to shallowest.
"""
function extract_daily_column_diag(simulation, diag_name)
    writer = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            writer = d.output_writer
            break
        end
    end
    isnothing(writer) && error("Diagnostic '$(diag_name)' not found.")
    (times, raw_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        writer, diag_name)
    dates = Date.(times isa Vector{DateTime} ? times : date.(times))
    n_days = length(dates)

    # Unwrap: each time-step entry is either a scalar or a ClimaCore Field.
    # Use parent() to extract the underlying Float64 array.
    first_val = raw_data[1]
    if isa(first_val, Number)
        # Unexpectedly scalar — shape as (1 × n_days)
        return dates, reshape(Float64.(raw_data), 1, n_days)
    else
        # Column field — unwrap with parent()
        n_z = length(parent(first_val))
        mat = Matrix{Float64}(undef, n_z, n_days)
        for t in 1:n_days
            mat[:, t] = Float64.(vec(parent(raw_data[t])))
        end
        return dates, mat
    end
end

# ── Observation Map (not used for CalLMIP sims, provided for compatibility) ───
function ClimaCalibrate.observation_map(iteration)
    return zeros(1, 1)   # placeholder — CalLMIP runs do not need EKP
end

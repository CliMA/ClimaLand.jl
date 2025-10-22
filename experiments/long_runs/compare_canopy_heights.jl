# # Comparison: Default vs Spatially-Varying Canopy Height

# This script runs two simulations to compare:
# 1. Default configuration: canopy_height = 1m everywhere
# 2. New configuration: spatially-varying canopy height from CLM data (capped at 8m)

# Run on clima with:
# srun --gpus=1 --mpi=none -t 03:00:00 --pty bash -l
# module load climacommon
# export CLIMACOMMS_DEVICE=CUDA
# export CLIMACOMMS_CONTEXT=SINGLETON
# cd /home/renatob/ClimaLand.jl
# julia --color=yes --project=.buildkite experiments/long_runs/compare_canopy_heights.jl

import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
using ClimaLand.Canopy: clm_canopy_height, effective_canopy_height, PrescribedBiomassModel
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

using Dates

const FT = Float64

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

function setup_model(FT, start_date, stop_date, Δt, domain, toml_dict, use_spatially_varying_height::Bool)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    
    # Forcing data
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    # Plant hydraulics
    conductivity_model = Canopy.PlantHydraulics.Weibull(toml_dict)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve(toml_dict)
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        surface_domain,
        toml_dict;
        conductivity_model,
        retention_model,
    )

    # Roughness lengths
    h_canopy = hydraulics.compartment_surfaces[end]
    z_0m = FT(0.13) * h_canopy
    z_0b = FT(0.1) * z_0m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    # Construct the P model manually
    photosynthesis = PModel{FT}(domain, toml_dict)
    conductance = PModelConductance{FT}(toml_dict)
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(domain, toml_dict)

    # Choose biomass model based on configuration
    if use_spatially_varying_height
        @info "Using spatially-varying canopy height from CLM data"
        # Read spatially-varying canopy height from CLM data
        raw_canopy_height = clm_canopy_height(surface_space)
        # Cap height to stay below atmospheric reference height (ERA5 at 10m)
        canopy_height = effective_canopy_height(raw_canopy_height, FT(10.0); buffer=FT(2.0))
        
        # Create biomass model with spatially-varying height
        biomass = PrescribedBiomassModel{FT}(surface_domain, LAI, toml_dict; height=canopy_height)
    else
        @info "Using default canopy height = 1m everywhere"
        # Use default biomass model (height = 1m from toml_dict)
        biomass = PrescribedBiomassModel{FT}(surface_domain, LAI, toml_dict)
    end

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        hydraulics,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        z_0m,
        z_0b,
        biomass,
    )

    # Snow model setup
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
    horz_degree_res =
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2
    scf = Snow.WuWuSnowCoverFractionModel(toml_dict, horz_degree_res)
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        α_snow,
        scf,
    )

    # Construct the land model
    land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt; snow, canopy)
    return land
end

# Simulation configuration: 1 year from March 2008 to March 2009
start_date = DateTime("2008-03-01")
stop_date = DateTime("2009-03-01")
Δt = 450.0
nelements = (101, 15)

domain = ClimaLand.Domains.global_domain(
    FT;
    context,
    nelements,
    mask_threshold = FT(0.99),
)
toml_dict = LP.create_toml_dict(FT)

@info "=========================================="
@info "COMPARISON RUN: Default vs Spatially-Varying Canopy Height"
@info "Duration: $(start_date) to $(stop_date)"
@info "Resolution: $nelements"
@info "Timestep: $Δt s"
@info "=========================================="

# Run 1: Default configuration (height = 1m)
@info "\n=== RUN 1: Default Configuration (height = 1m) ==="
root_path_default = "comparison_default_height_$(device_suffix)"
diagnostics_outdir_default = joinpath(root_path_default, "global_diagnostics")
outdir_default = ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir_default)

model_default = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict, false)
simulation_default = LandSimulation(start_date, stop_date, Δt, model_default; outdir = outdir_default)

CP.log_parameter_information(toml_dict, joinpath(root_path_default, "parameters.toml"))
@time ClimaLand.Simulations.solve!(simulation_default)

@info "Default run completed. Output saved to: $root_path_default"

# Run 2: Spatially-varying height
@info "\n=== RUN 2: Spatially-Varying Canopy Height ==="
root_path_varying = "comparison_varying_height_$(device_suffix)"
diagnostics_outdir_varying = joinpath(root_path_varying, "global_diagnostics")
outdir_varying = ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir_varying)

model_varying = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict, true)
simulation_varying = LandSimulation(start_date, stop_date, Δt, model_varying; outdir = outdir_varying)

CP.log_parameter_information(toml_dict, joinpath(root_path_varying, "parameters.toml"))
@time ClimaLand.Simulations.solve!(simulation_varying)

@info "Spatially-varying height run completed. Output saved to: $root_path_varying"

@info "\n=========================================="
@info "COMPARISON COMPLETED"
@info "Default results: $root_path_default"
@info "Varying height results: $root_path_varying"
@info "=========================================="

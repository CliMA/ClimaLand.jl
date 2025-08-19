# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, for 6 hours on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# This code also assesses performance, either via Nsight or by running the
# model multiple times and collecting statistics for execution time and allocations
# to make flame graphs. You can choose between the two
# by requestion --profiler nsight or --profiler flamegraph. If not provided,
# flamegraphs are created.

# When run with buildkite on clima, without Nisght, this code also compares with the previous best time
# saved at the bottom of this file

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 6 hours
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every Newton iteration
# Atmos forcing update: every 3 hours
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaTimeSteppers as CTS
import ClimaCore
@show pkgversion(ClimaCore)
using ClimaUtilities.ClimaArtifacts

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using Dates
import Profile, ProfileCanvas
using Test
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--profiler"
        help = "Profiler option: nsight or flamegraph"
        arg_type = String
        default = "flamegraph"
    end
    return parse_args(s)
end

const FT = Float64

const PREVIOUS_BEST_GPU_TIME = 0.67

include("benchmark_sim.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "snowy_land_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)
parsed_args = parse_commandline()
profiler = parsed_args["profiler"]
duration = profiler == "nsight" ? Dates.Minute(90) : Dates.Hour(6)
time_interpolation_method = LinearInterpolation(PeriodicCalendar())

function setup_simulation(duration)
    start_date = DateTime(2008)
    Δt = 450.0
    end_date = start_date + duration
    updateat = Array(start_date:Second(3600 * 3):end_date)
    nelements = (101, 15)
    earth_param_set = LP.LandParameters(FT)
    domain = ClimaLand.Domains.global_domain(FT; nelements = nelements)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    start_date = DateTime(2008)

    # Forcing data
    era5_ncdata_path = ClimaLand.Artifacts.era5_land_forcing_data2008_path(;
        context,
        lowres = true,
    )
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
        context = nothing,
        start_date = start_date,
        end_date = start_date,
    )
    LAI = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = time_interpolation_method,
    )

    # Overwrite some defaults for the canopy model
    # Energy model
    ac_canopy = FT(2.5e3)
    energy = Canopy.BigLeafEnergyModel{FT}(; ac_canopy)

    # Roughness lengths
    hydraulics = Canopy.PlantHydraulicsModel{FT}(surface_domain, LAI)
    h_canopy = hydraulics.compartment_surfaces[end]
    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        earth_param_set;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        energy,
        hydraulics,
        z_0m = z0_m,
        z_0b = z0_b,
    )

    # Construct land model with all default components
    land = LandModel{FT}(forcing, LAI, earth_param_set, domain, Δt; canopy)


    # Set initial conditions
    function set_ic!(Y, p, t0, land)
        (; θ_r, ν, ρc_ds) = land.soil.parameters
        @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
        Y.soil.θ_i .= 0
        T = FT(276.85)
        ρc_s =
            Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                T,
                earth_param_set,
            )
        Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air

        plant_ν = land.canopy.hydraulics.parameters.ν
        Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
        evaluate!(Y.canopy.energy.T, atmos.T, t0)

        Y.snow.S .= 0
        Y.snow.S_l .= 0
        Y.snow.U .= 0
    end
    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        end_date,
        Δt,
        land;
        updateat = updateat,
        set_ic!,
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end
run_benchmarks(
    device,
    () -> setup_simulation(duration),
    "integrated",
    PREVIOUS_BEST_GPU_TIME,
    outdir,
)

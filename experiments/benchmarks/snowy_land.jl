# # Global run of land model

# The code sets up and profiles ClimaLand v1, which
# includes soil, canopy, and snow, for 1 day (2 hours without gpu) on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# See "benchmark_sim.jl" for details on what is benchmarked

# When run with buildkite on clima, without Nisght, this code compares with the previous
# best time saved at the top of this file

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 1 day or 2 hours
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every Newton iteration
# Atmos forcing update: every 3 hours
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
using Dates
using CSV
const FT = Float64

######################################################################
## This result is from a benchmark ran on an A100 on the clima cluster
const PREVIOUS_GPU_TIME = 2.67
######################################################################

include("benchmark_sim.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "snowy_land_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)
parsed_args = parse_commandline()
profiler = parsed_args["profiler"]

function setup_simulation()
    default_params_filepath =
        joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
    toml_dict = LP.create_toml_dict(FT, default_params_filepath)
    start_date = DateTime(2008)
    duration =
        device isa ClimaComms.CPUSingleThreaded ? Dates.Hour(2) : Dates.Day(1)
    Δt = 450.0
    stop_date = start_date + duration
    updateat = Array(start_date:Second(3600 * 3):stop_date)
    time_interpolation_method = LinearInterpolation(PeriodicCalendar())
    nelements = (101, 15)
    earth_param_set = LP.LandParameters(FT)
    domain = ClimaLand.Domains.global_domain(FT; nelements = nelements)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface

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
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    # Overwrite some defaults for the canopy model
    # Energy model
    ac_canopy = FT(2.5e3)
    energy = Canopy.BigLeafEnergyModel{FT}(; ac_canopy)

    # Roughness lengths
    hydraulics = Canopy.PlantHydraulicsModel{FT}(surface_domain, LAI, toml_dict)
    h_canopy = hydraulics.compartment_surfaces[end]
    z_0m = FT(0.13) * h_canopy
    z_0b = FT(0.1) * z_0m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        energy,
        hydraulics,
        z_0m,
        z_0b,
    )

    # Construct land model with all default components
    land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt; canopy)


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
        stop_date,
        Δt,
        land;
        updateat = updateat,
        set_ic!,
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end
run_benchmarks(device, setup_simulation, profiler, PREVIOUS_GPU_TIME, outdir)

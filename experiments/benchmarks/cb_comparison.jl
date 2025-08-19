delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaTimeSteppers as CTS
import ClimaCore
@show pkgversion(ClimaCore)
using ClimaUtilities.ClimaArtifacts
import ClimaDiagnostics
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities.TimeManager: ITime, date
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using Dates
using Test
using ArgParse

include("benchmark_sim.jl")



const FT = Float64


include("benchmark_sim.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()

function generic_setup_simulation(;
    diagnostic_frequency = nothing,
    nan_check_frequency = nothing,
)
    start_date = DateTime(2008)
    duration = Dates.Hour(45)
    Δt = 450.0
    end_date = start_date + duration
    updateat = Array(start_date:Second(3600 * 3):end_date)
    time_interpolation_method = LinearInterpolation(PeriodicCalendar())
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
    output_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        domain.space.subsurface,
        tempdir();
        start_date,
    )
    diagnostics =
        isnothing(diagnostic_frequency) ? nothing :
        ClimaLand.default_diagnostics(
            land,
            start_date;
            output_writer,
            average_period = diagnostic_frequency,
        )
    user_callbacks =
        isnothing(nan_check_frequency) ? () :
        (
            ClimaLand.NaNCheckCallback(
                nan_check_frequency,
                ITime(0; epoch = start_date),
                nothing,
            ),
        )
    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        end_date,
        Δt,
        land;
        updateat,
        set_ic!,
        user_callbacks,
        diagnostics,
    )
    ClimaLand.Simulations.step!(simulation)
    return simulation
end

setup_with_daily() = generic_setup_simulation(; diagnostic_frequency = :daily)
setup_with_monthly() =
    generic_setup_simulation(; diagnostic_frequency = :monthly)
setup_with_half_hourly() =
    generic_setup_simulation(; diagnostic_frequency = :halfhourly)
setup_with_infrequent_nan_check() =
    generic_setup_simulation(; nan_check_frequency = Day(1))
setup_with_frequent_nan_check() =
    generic_setup_simulation(; nan_check_frequency = Hour(10))

# generic_setup_simulation()
# setup_with_daily()
# setup_with_half_hourly()
# setup_with_monthly()
# setup_with_infrequent_nan_check()
# setup_with_frequent_nan_check()

run_timing_benchmarks(
    device,
    generic_setup_simulation;
    MAX_PROFILING_TIME_SECONDS = 280,
)

run_timing_benchmarks(
    device,
    setup_with_monthly;
    MAX_PROFILING_TIME_SECONDS = 280,
)


run_timing_benchmarks(
    device,
    setup_with_daily;
    MAX_PROFILING_TIME_SECONDS = 280,
)


run_timing_benchmarks(
    device,
    setup_with_infrequent_nan_check;
    MAX_PROFILING_TIME_SECONDS = 280,
)

# run_timing_benchmarks(device, setup_with_frequent_nan_check)

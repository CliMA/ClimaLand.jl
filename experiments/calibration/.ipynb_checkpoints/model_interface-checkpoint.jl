# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 730 d
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours

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
using ClimaLand.Canopy: uSPACPiParameters, uSPACConductancePi
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

using Dates

import ClimaCalibrate
import JLD2
import EnsembleKalmanProcesses as EKP

const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "snowy_land_longrun_$(device_suffix)"

function setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data - always use high resolution for calibration runs
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        use_lowres_forcing = true, # change to true to use just 2008, false for all years
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

    # Overwrite some defaults for the canopy model
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

    # ---------- uSPAC Π-group conductance wired to α/β params ----------
    # These keys must match your EKP prior names written into the override TOMLs.
    αR  = FT(CP.get_parameter(toml_dict, "alpha_R"))
    βR  = FT(CP.get_parameter(toml_dict, "beta_R"))
    αF  = FT(CP.get_parameter(toml_dict, "alpha_F"))
    βF  = FT(CP.get_parameter(toml_dict, "beta_F"))
    αT  = FT(CP.get_parameter(toml_dict, "alpha_T"))
    βTs = FT(CP.get_parameter(toml_dict, "beta_Ts"))
    αS  = FT(CP.get_parameter(toml_dict, "alpha_S"))
    βSs = FT(CP.get_parameter(toml_dict, "beta_Ss"))
    
    ΓR = (αR,  βR)    # aridity → ΠR
    ΓF = (αF,  βF)    # aridity → ΠF
    ΓT = (αT,  βTs)   # sand    → ΠT
    ΓS = (αS,  βSs)   # sand    → ΠS

    uspac_pars = uSPACPiParameters{FT}(; ΓR, ΓF, ΓT, ΓS)
    conductance = uSPACConductancePi{FT}(uspac_pars)

    @info "uSPAC conductance active" ΓR ΓF ΓT ΓS
        
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        hydraulics,
        z_0m,
        z_0b,
        conductance = conductance, # for uspac
    )

    # Snow model setup
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
    horz_degree_res =
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2 # mean of resolution in latitude and longitude, in degrees
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

    # Construct the land model with all default components except for snow
    land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt; snow, canopy)
    return land
end

function ClimaCalibrate.forward_model(iteration, member)
    output_dir = CALIBRATE_CONFIG.output_dir
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    nelements = CALIBRATE_CONFIG.nelements
    spinup = CALIBRATE_CONFIG.spinup
    extend = CALIBRATE_CONFIG.extend
    ensemble_member_path =
        ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, member)

    eki = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    minibatch = EKP.get_current_minibatch(eki)

    # Determine start date and end date from the sample date ranges
    start_date = first(sample_date_ranges[minimum(minibatch)]) - spinup
    stop_date = last(sample_date_ranges[maximum(minibatch)]) + extend
    Δt = 450.0

    # Convert to ITimes
    t0 = ITime(0, Dates.Second(1), start_date)
    tf = ITime(
        Dates.value(convert(Dates.Second, stop_date - start_date)),
        epoch = start_date,
    )
    Δt = ITime(Δt, epoch = start_date)

    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    outdir =
        ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_dir)

    domain = ClimaLand.Domains.global_domain(
        FT;
        context,
        nelements,
        mask_threshold = FT(0.99),
    )

    calibrate_params_path =
        ClimaCalibrate.parameter_path(output_dir, iteration, member)
    toml_dict =
        LP.create_toml_dict(FT, override_files = [calibrate_params_path])

    model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)

    # Set up diagnostics
    domain = ClimaLand.get_domain(model)
    diagnostic_domain =
        haskey(domain.space, :subsurface) ? domain.space.subsurface :
        domain.space.surface
    output_writer =
        ClimaDiagnostics.NetCDFWriter(diagnostic_domain, outdir; start_date)

    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date;
        output_writer,
        output_vars = ["lhf", "shf", "lwu", "swu"], # add soil surface e.g. "soil_moisture_surface" when smapping
        reduction_period = :monthly,
        reduction_type = :average,
    )

    simulation = LandSimulation(
        t0,
        tf,
        Δt,
        model;
        outdir,
        user_callbacks = (ClimaLand.ReportCallback(div((tf - t0), 10), t0),),
        diagnostics = diagnostics,
    )
    @info "Run: Global Soil-Canopy-Snow Model"
    @info "Resolution: $nelements"
    @info "Timestep: $Δt s"
    @info "Start Date: $start_date"
    @info "Stop Date: $stop_date"
    CP.log_parameter_information(
        toml_dict,
        joinpath(ensemble_member_path, "log_params_$member.toml"),
    )
    ClimaLand.Simulations.solve!(simulation)
    return nothing
end

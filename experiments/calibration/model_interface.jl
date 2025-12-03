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
import ClimaUtilities.TimeVaryingInputs
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
using Statistics

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

    # ========== Load CWD field ==========
    # cwd_file = joinpath(@__DIR__, "cwd_era5.jld2")


    # cwd_data = JLD2.load(cwd_file)
    # lons_cwd = cwd_data["lons"]
    # lats_cwd = cwd_data["lats"]
    # CWD_spatial = cwd_data["CWD_field"]
    
    # CWD_field = regrid_cwd_to_model_space(lons_cwd, lats_cwd, CWD_spatial, surface_space, FT)


    """
        load_cwd_field(surface_space, FT)

    Load cumulative water deficit (CWD) from a JLD2 file and regrid to model surface space.
    CWD should be a spatial field (2D: lon × lat) representing long-term water stress.
    """
    function load_cwd_field(surface_space, FT)
        cwd_file = joinpath(
            pkgdir(ClimaLand),
            "experiments",
            "calibration",
            "cwd_era5.jld2"
        )
        
        if !isfile(cwd_file)
            @warn "CWD file not found at $cwd_file, using zero CWD"
            return ClimaCore.Fields.zeros(FT, surface_space)
        end
        
        @info "Loading CWD from $cwd_file"
        
        # Load CWD data from JLD2
        cwd_data_dict = JLD2.load(cwd_file)
        
        # file structure:
        # "CWD_field" => Float64 array (lon, lat)
        # "lons" => Float64 vector
        # "lats" => Float64 vector
        
        if !haskey(cwd_data_dict, "CWD_field")
            @error "CWD file missing 'CWD_field' key. Available keys: $(keys(cwd_data_dict))"
            return ClimaCore.Fields.zeros(FT, surface_space)
        end

        cwd_field = cwd_data_dict["CWD_field"]

        # Transfer to GPU if needed
        # Check if the model's surface_space is on GPU
        if parent(ClimaCore.Fields.zeros(FT, surface_space)) isa CUDA.CuArray
            @info "Transferring CWD field to GPU"
            # Create a new field on GPU and copy data
            cwd_field_gpu = similar(ClimaCore.Fields.zeros(FT, surface_space))
            parent(cwd_field_gpu) .= CUDA.CuArray(parent(cwd_field))
            cwd_field = cwd_field_gpu
        end

        @info "CWD field loaded successfully"
        
        return cwd_field
    end
    
    # Load CWD (Cumulative Water Deficit) map
    CWD_field = load_cwd_field(surface_space, FT)

    # ========== Load or calculate SM variability field ==========
    # Use coefficient of variation (CV) for normalized variability across climates
    # Note: SMAP (9km) is coarser than ERA5 (0.25°), but CV captures temporal
    # dynamics within each pixel, serving as a landscape-scale proxy for soil
    # moisture heterogeneity. High CV indicates dynamic conditions (e.g., shallow
    # soils, high drainage), while low CV indicates buffered conditions.
    # sm_metric = :cv
    
    # # SM variability filename
    # if start_year == stop_year
    #     sm_var_file = joinpath(@__DIR__, "sm_variability_smap_$(start_year).jld2")
    # else
    #     sm_var_file = joinpath(@__DIR__, "sm_variability_smap_$(start_year)_$(stop_year).jld2")
    # end
    
    # if !isfile(sm_var_file)
    #     @warn "SM variability file not found, calculating from SMAP data" sm_var_file
        
    #     # Get SMAP data directory from environment or use default
    #     smap_data_dir = get(ENV, "SMAP_DATA_PATH", "/resnick/scratch/egreich/SMAP/SMAP_L3_SM_P_E")
        
    #     # Find SMAP files in the date range
    #     smap_files = find_smap_files(start_date, stop_date; data_dir=smap_data_dir)
        
    #     if isempty(smap_files)
    #         @warn "No SMAP data available in $smap_data_dir, setting default SM variability values"
    #         # Use a reasonable default CV (e.g., 0.3 = 30% variability)
    #         sm_var_field = ClimaCore.Fields.Field(FT(0.3), surface_space)
    #     else
    #         @info "Found $(length(smap_files)) SMAP files for SM variability calculation"
            
    #         # Calculate SM variability from SMAP files
    #         # This will create gridded maps of std, CV, and mean
    #         lons_sm, lats_sm, sm_std_spatial, sm_cv_spatial, sm_mean_spatial = 
    #             calculate_sm_variability_from_smap(
    #                 smap_files,
    #                 start_date,
    #                 stop_date,
    #             )
            
    #         # Use CV as primary metric
    #         sm_var_spatial = sm_cv_spatial
            
    #         jldsave(
    #             sm_var_file; 
    #             lons=lons_sm, 
    #             lats=lats_sm, 
    #             sm_std_spatial,
    #             sm_cv_spatial,
    #             sm_mean_spatial,
    #             sm_var_spatial,
    #             metric=sm_metric,
    #             start_year, 
    #             stop_year
    #         )
    #         @info "Saved SM variability (CV) to" sm_var_file
            
    #         sm_var_field = regrid_sm_variability_to_model_space(lons_sm, lats_sm, sm_var_spatial, surface_space, FT)
    #     end
    # else
    #     @info "Loading pre-calculated SM variability data" sm_var_file
    #     sm_var_data = JLD2.load(sm_var_file)
    #     lons_sm = sm_var_data["lons"]
    #     lats_sm = sm_var_data["lats"]
        
    #     # Prefer CV if available, otherwise use stored primary metric
    #     if haskey(sm_var_data, "sm_cv_spatial")
    #         sm_var_spatial = sm_var_data["sm_cv_spatial"]
    #     elseif haskey(sm_var_data, "sm_var_spatial")
    #         sm_var_spatial = sm_var_data["sm_var_spatial"]
    #     else
    #         sm_var_spatial = sm_var_data["sm_std_spatial"]
    #     end
        
    #     sm_var_field = regrid_sm_variability_to_model_space(lons_sm, lats_sm, sm_var_spatial, surface_space, FT)
    # end
    
    # Now you have both CWD_field and sm_var_field (CV) as covariates
    # CWD: spatial gradient in long-term water deficit (climate-driven)
    # SM_CV: spatial gradient in soil moisture variability (soil/landscape-driven)
    @info "Covariate fields loaded"

    # Forcing data - always use high resolution for calibration runs
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        use_lowres_forcing = false, # change to true to use just 2008, false for all years
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
    # Helper function to extract scalar value from ClimaParams output
    function extract_param(toml_dict, name)
        val = CP.get_parameter_values(toml_dict, name)
        # Handle both NamedTuple and direct numeric returns
        if val isa NamedTuple
            # If it's a NamedTuple, try to get the field with the same name
            return FT(getproperty(val, Symbol(name)))
        else
            return FT(val)
        end
    end
    
    # These keys must match your EKP prior names written into the override TOMLs.
    # NOTE: Prior means should NOT be 0.0, as this causes division-by-zero in stomatalconductance.jl
    # when all π-groups (ΠR, ΠF, ΠT, ΠS) = 0. Safeguards have been added to those functions,
    # but starting with non-zero priors (e.g., 0.1) is better for numerical stability.
    αR  = extract_param(toml_dict, "alpha_R")
    βR  = extract_param(toml_dict, "beta_R")
    αF  = extract_param(toml_dict, "alpha_F")
    βF  = extract_param(toml_dict, "beta_F")
    αT  = extract_param(toml_dict, "alpha_T")
    βTs = extract_param(toml_dict, "beta_Ts")
    αS  = extract_param(toml_dict, "alpha_S")
    βSs = extract_param(toml_dict, "beta_Ss")

    ΓR = (αR,  βR)    # aridity → ΠR
    ΓF = (αF,  βF)    # aridity → ΠF
    ΓT = (αT,  βTs)   # sand    → ΠT
    ΓS = (αS,  βSs)   # sand    → ΠS

    # Pass covariates directly to uSPAC parameters
    # CWD_field from forcing (calculated from ERA5 data)
    # sand will be extracted from soil parameters internally in stomatalconductance.jl
    aridity_idx = CWD_field  # Climatic Water Deficit in mm
    sand = nothing  # Will be extracted from p.soil in update_canopy_conductance!
    
    uspac_pars = uSPACPiParameters{FT}(; ΓR, ΓF, ΓT, ΓS, aridity_idx, sand)
    conductance = uSPACConductancePi{FT}(uspac_pars)

    @info "uSPAC conductance active" ΓR ΓF ΓT ΓS aridity_idx
        
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

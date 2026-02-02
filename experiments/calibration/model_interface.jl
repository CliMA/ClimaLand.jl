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
    # ========== Load Aridity Field ==========
    # Aridity is loaded from pre-processed Global Aridity Index (aridity.jld2)

    """
        load_aridity_field(surface_space, domain, FT)

    Load aridity index from a JLD2 file and regrid to model surface space.
    Aridity is P/ET0 where higher values indicate wetter conditions.
    Ocean points are masked to NaN using the topographic land/sea mask.
    """
    function load_aridity_field(surface_space, domain, FT)
        aridity_file = joinpath(
            pkgdir(ClimaLand),
            "experiments",
            "calibration",
            "aridity.jld2"
        )
        
        if !isfile(aridity_file)
            @warn "Aridity file not found at $aridity_file, using zero aridity"
            return ClimaCore.Fields.zeros(FT, surface_space)
        end
        
        @info "Loading aridity from $aridity_file"
        
        # Load aridity data from JLD2
        aridity_data_dict = JLD2.load(aridity_file)
        
        # file structure:
        # "CWD_field" => Float64 array (lon, lat) [kept for compatibility]
        # "lons" => Float64 vector
        # "lats" => Float64 vector
        
        if !haskey(aridity_data_dict, "CWD_field")
            @error "Aridity file missing 'CWD_field' key. Available keys: $(keys(aridity_data_dict))"
            return ClimaCore.Fields.zeros(FT, surface_space)
        end

        aridity_field = aridity_data_dict["CWD_field"]  # Note: still called CWD_field for compatibility

        # Transfer to GPU if needed
        # Check if the model's surface_space is on GPU
        if parent(ClimaCore.Fields.zeros(FT, surface_space)) isa CUDA.CuArray
            @info "Transferring aridity field to GPU"
            # Create a new field on GPU and copy data
            aridity_field_gpu = similar(ClimaCore.Fields.zeros(FT, surface_space))
            parent(aridity_field_gpu) .= CUDA.CuArray(parent(aridity_field))
            aridity_field = aridity_field_gpu
        end

        # Apply topographic land/sea mask to set ocean points to NaN
        # This ensures ocean masking is consistent with ClimaLand's standard approach
        @info "Applying land/sea mask to aridity field..."
        landsea_mask_field = ClimaLand.Domains.landsea_mask(domain; threshold = 0.5)
        
        # Mask ocean points (where landsea_mask == 0) to NaN
        # This is the authoritative way to identify ocean vs land
        aridity_field = @. ifelse(landsea_mask_field == FT(0), FT(NaN), aridity_field)
        
        # Count ocean vs land points
        n_total = length(parent(aridity_field))
        n_ocean = sum(parent(landsea_mask_field) .== FT(0))
        n_land = n_total - n_ocean
        @info "Ocean masking complete: $(n_ocean)/$(n_total) points are ocean ($(round(100*n_ocean/n_total, digits=1))%)"
        @info "Land points: $(n_land) ($(round(100*n_land/n_total, digits=1))%)"
        
        return aridity_field
    end
    
    # Load aridity field (Global Aridity Index: higher = wetter)
    # Ocean points are masked to NaN using topographic land/sea mask
    aridity_field = load_aridity_field(surface_space, domain, FT)

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

    # -------------------- Extract calibrated parameters --------------------
    
    # Extract from TOML
    function extract_param(toml_dict, name)
        val = CP.get_parameter_values(toml_dict, name)
        if val isa NamedTuple
            return FT(getproperty(val, Symbol(name)))
        else
            return FT(val)
        end
    end
    
    βkx_base  = extract_param(toml_dict, "βkx_base")
    βkx_coord  = extract_param(toml_dict, "βkx_coord")
    βψx50_base  = extract_param(toml_dict, "βψx50_base")
    βψx50_slope  = extract_param(toml_dict, "βψx50_slope")
    βΠR_base  = extract_param(toml_dict, "βΠR_base")
    βΠR_slope  = extract_param(toml_dict, "βΠR_slope")

    # -------------------- Build uSPAC conductance model --------------------
    # Note: Soil parameters (K_sat, ψ_sat, RAI, Zr) are extracted at runtime
    # from p.soil within pi_groups_from_calibrated_traits()
    
    uspac_pars = uSPACPiParameters{FT}(;
        # Calibrated (aridity-dependent):
        βkx_base = βkx_base,
        βkx_coord = βkx_coord,
        βψx50_base = βψx50_base,
        βψx50_slope = βψx50_slope,
        βΠR_base = βΠR_base,
        βΠR_slope = βΠR_slope,
        
        # Covariates (only aridity needed now):
        aridity_idx = aridity_field,
        
        # Trait distribution parameters (Phase 1.5: Climate-dependent variance):
        # Enable trait heterogeneity with climate-dependent variance
        use_trait_distribution = false,
        n_quad = 3,  # Gaussian-Hermite quadrature order (3 points per dimension, efficient)
        
        # Base variance (wet climate reference from literature)
        σ_logkx_base = FT(0.5),   # Liu et al. 2019
        σ_P50_base = FT(1.5),     # Choat et al. 2012
        σ_ΠR_base = FT(0.15),     # Strategy variation
        
        # Climate-dependent variance modifiers (dry → more variable)
        α_σ_logkx = FT(0.3),  # Aridity increases kx variance (niche partitioning)
        α_σ_P50 = FT(1.0),    # Aridity increases P50 variance (drought adaptation spectrum)
        α_σ_ΠR = FT(0.1),     # Aridity increases strategy variance (iso/anisohydric diversity)
        
        # Trait correlations (climate-dependent)
        ρ_kx_P50_base = FT(0.7),   # Safety-efficiency tradeoff (Manzoni et al. 2013)
        ρ_kx_ΠR_base = FT(-0.3),   # High kx requires tight regulation
        ρ_P50_ΠR_base = FT(-0.2),  # Resistant plants can be anisohydric
        α_ρ_kx_P50 = FT(-0.2),     # Weaker coordination in dry climates
        
        # Fixed parameters:
        β_star_frac = FT(0.95),
        β_w_frac    = FT(0.05),
        gsw_max = FT(Inf)
    )
    
    conductance = uSPACConductancePi{FT}(uspac_pars)

    @info "uSPAC conductance active with trait distribution integration" use_trait_distribution=uspac_pars.use_trait_distribution n_quad=uspac_pars.n_quad
        
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil),
        hydraulics,
        z_0m,
        z_0b,
        conductance = conductance,
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
        prognostic_land_components = (:canopy, :snow, :soil),
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
        output_vars = ["lhf", "shf", "lwu", "swu", "swc"], # swc = soil water content (depth-resolved)
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

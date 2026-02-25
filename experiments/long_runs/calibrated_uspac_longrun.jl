# # Global run of land model with calibrated uSPAC hydraulic traits

# This code sets up and runs ClimaLand v1, which includes soil, canopy, and snow,
# on a spherical domain, using ERA5 data as forcing. This version uses the
# calibrated hydraulic trait coordination parameters (iteration_002) that are
# now stored in toml/default_parameters.toml.

# Key differences from standard long runs:
# - Uses uSPACConductancePi instead of PModelConductance or standard PlantHydraulics
# - Loads Global Aridity Index as a spatial covariate for trait coordination
# - Calibrated β parameters automatically loaded from default_parameters.toml
# - Validated against TRY database (P50 median difference: 0.07 MPa)

# Simulation Setup
# Number of spatial elements: determined by global_box_domain()
# Soil depth: 50 m
# Simulation duration: configurable (2 years default, 19 years for LONGER_RUN)
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
using JLD2

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64;

# Set LONGER_RUN=true for 19-year simulation, false for 2-year test
const LONGER_RUN = haskey(ENV, "LONGER_RUN") ? true : false
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "calibrated_uspac_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    
    # ========== Load Aridity Field (required for trait coordination) ==========
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
        
        if !haskey(aridity_data_dict, "CWD_field")
            @error "Aridity file missing 'CWD_field' key. Available keys: $(keys(aridity_data_dict))"
            return ClimaCore.Fields.zeros(FT, surface_space)
        end

        aridity_field = aridity_data_dict["CWD_field"]

        # Apply topographic land/sea mask to set ocean points to NaN
        @info "Applying land/sea mask to aridity field..."
        landsea_mask_field = ClimaLand.Domains.landsea_mask(domain; threshold = 0.5)
        
        # Mask ocean points (where landsea_mask == 0) to NaN
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
    aridity_field = load_aridity_field(surface_space, domain, FT)

    # ========== Forcing Data ==========
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

    # ========== Plant Hydraulics Setup ==========
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

    # ========== Extract Calibrated Parameters from TOML ==========
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
    
    @info "Loaded calibrated hydraulic trait parameters:" βkx_base βkx_coord βψx50_base βψx50_slope βΠR_base βΠR_slope

    # ========== Build uSPAC Conductance Model ==========
    uspac_pars = uSPACPiParameters{FT}(;
        # Calibrated parameters (aridity-dependent):
        βkx_base = βkx_base,
        βkx_coord = βkx_coord,
        βψx50_base = βψx50_base,
        βψx50_slope = βψx50_slope,
        βΠR_base = βΠR_base,
        βΠR_slope = βΠR_slope,
        
        # Covariates (required for trait coordination):
        aridity_idx = aridity_field,
        
        # Trait distribution parameters
        # Set use_trait_distribution=true to enable trait heterogeneity
        use_trait_distribution = false,  # Set to true for trait distribution integration
        n_quad = 3,  # Gaussian-Hermite quadrature order
        
        # Base variance (wet climate reference from literature)
        σ_logkx_base = FT(0.5),   # Liu et al. 2019
        σ_P50_base = FT(1.5),     # Choat et al. 2012
        σ_ΠR_base = FT(0.15),     # Strategy variation
        
        # Climate-dependent variance modifiers (dry → more variable)
        α_σ_logkx = FT(0.3),
        α_σ_P50 = FT(1.0),
        α_σ_ΠR = FT(0.1),
        
        # Trait correlations (climate-dependent)
        ρ_kx_P50_base = FT(0.7),   # Safety-efficiency tradeoff
        ρ_kx_ΠR_base = FT(-0.3),
        ρ_P50_ΠR_base = FT(-0.2),
        α_ρ_kx_P50 = FT(-0.2),
        
        # Fixed parameters:
        β_star_frac = FT(0.95),
        β_w_frac    = FT(0.05),
        gsw_max = FT(Inf)
    )
    
    conductance = uSPACConductancePi{FT}(uspac_pars)
    @info "uSPAC conductance model initialized" use_trait_distribution=uspac_pars.use_trait_distribution

    # ========== Canopy Model ==========
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil),
        hydraulics,
        z_0m,
        z_0b,
        conductance = conductance,  # Use calibrated uSPAC conductance
    )

    # ========== Snow Model ==========
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
        prognostic_land_components = (:canopy, :snow, :soil),
        α_snow,
        scf,
    )

    # ========== Complete Land Model ==========
    land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt; snow, canopy)
    return land
end

# ========== Simulation Configuration ==========
# If not LONGER_RUN, run for 2 years; note that the forcing from 2008 is repeated.
# If LONGER run, run for 19 years, with the correct forcing each year.
start_date = LONGER_RUN ? DateTime("2000-03-01") : DateTime("2008-03-01")
stop_date = LONGER_RUN ? DateTime("2019-03-01") : DateTime("2010-03-01")
Δt = 450.0

# ========== Domain and Model Setup ==========
domain = ClimaLand.Domains.global_box_domain(
    FT;
    context,
    mask_threshold = FT(0.99),
)

# Load parameters from default TOML (includes calibrated β parameters)
toml_dict = LP.create_toml_dict(FT)

# Build the model
model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)

# Create simulation
simulation = LandSimulation(start_date, stop_date, Δt, model; outdir)

# ========== Run Information ==========
@info "="^70
@info "Global Soil-Canopy-Snow Model with Calibrated uSPAC Traits"
@info "="^70
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
@info "Duration: $(Dates.value(stop_date - start_date) / (365.25 * 86400)) years"
@info "Output Directory: $root_path"
@info "="^70
@info "Calibration Info:"
@info "  - Parameters from: iteration_002 (SMAP calibration)"
@info "  - TRY validation: P50 median diff = 0.07 MPa"
@info "  - Conductance model: uSPACConductancePi (trait coordination)"
@info "  - Aridity covariate: Global Aridity Index (P/ET0)"
@info "="^70

CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))

# ========== Run Simulation ==========
ClimaLand.Simulations.solve!(simulation)

# ========== Post-processing and Visualization ==========
@info "Creating visualizations..."
LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation; savedir = root_path, date = stop_date)
LandSimVis.make_leaderboard_plots(simulation; savedir = root_path)

@info "="^70
@info "Simulation complete!"
@info "Results saved to: $root_path"
@info "="^70

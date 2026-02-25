# Template for running ClimaLand with calibrated hydraulic trait coordination parameters
# This script shows how to set up a long run using the calibrated parameters
# from iteration_002 that are now stored in toml/default_parameters.toml

# The key changes from standard long run scripts are:
# 1. Use uSPACConductancePi instead of PModelConductance or other conductance models
# 2. Load aridity field as a covariate (required for trait coordination)
# 3. The calibrated β parameters are automatically loaded from default_parameters.toml

import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.Canopy: uSPACPiParameters, uSPACConductancePi
import ClimaLand.Parameters as LP

const FT = Float64
context = ClimaComms.context()

# ============================================================================
# STEP 1: Load parameters from default_parameters.toml
# ============================================================================
# This automatically loads the calibrated parameters:
# βkx_base, βkx_coord, βψx50_base, βψx50_slope, βΠR_base, βΠR_slope
toml_dict = LP.create_toml_dict(FT)

# ============================================================================
# STEP 2: Load aridity field (required covariate for trait coordination)
# ============================================================================
using JLD2

# Load global aridity index from file
# Note: The file contains "CWD_field" which is actually the aridity data
aridity_file = joinpath(
    @__DIR__,
    "aridity.jld2"
)

# Uncomment these lines when ready to interpolate aridity onto your grid:
# aridity_data = JLD2.load(aridity_file)
# cwd_field = aridity_data["CWD_field"]  # This is the aridity field on ClimaCore space
# # You'll need to extract the parent array and interpolate onto grid
# # See model_interface.jl lines 130-145 for interpolation example

println("Note: Aridity field loading is commented out in this template.")
println("See model_interface.jl for full implementation details.")

# ============================================================================
# STEP 3: Extract calibrated parameters and create uSPAC conductance model
# ============================================================================

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

println("Loaded calibrated parameters:")
println("  βkx_base   = $βkx_base")
println("  βkx_coord  = $βkx_coord")
println("  βψx50_base = $βψx50_base")
println("  βψx50_slope = $βψx50_slope")
println("  βΠR_base   = $βΠR_base")
println("  βΠR_slope  = $βΠR_slope")

# Create uSPAC parameters structure
uspac_pars = uSPACPiParameters{FT}(;
    # Calibrated (aridity-dependent):
    βkx_base = βkx_base,
    βkx_coord = βkx_coord,
    βψx50_base = βψx50_base,
    βψx50_slope = βψx50_slope,
    βΠR_base = βΠR_base,
    βΠR_slope = βΠR_slope,
    
    # Covariates (REQUIRED: set aridity_field to your interpolated aridity data):
    # aridity_idx = aridity_field,
    
    # Trait distribution parameters (optional - set to false for deterministic):
    use_trait_distribution = false,
    n_quad = 3,
    
    # Base variance (if using trait distributions):
    σ_logkx_base = FT(0.5),
    σ_P50_base = FT(1.5),
    σ_ΠR_base = FT(0.15),
    
    # Climate-dependent variance modifiers:
    α_σ_logkx = FT(0.3),
    α_σ_P50 = FT(0.5),
    α_σ_ΠR = FT(0.2),
)

# Create conductance model
# conductance = uSPACConductancePi{FT}(uspac_pars)

# ============================================================================
# STEP 4: Use in canopy model setup
# ============================================================================
# When setting up your CanopyModel, use the uSPAC conductance model:
#
# canopy = CanopyModel{FT}(;
#     conductance = conductance,
#     # ... other canopy parameters
# )
#
# Then include this canopy in your land model as usual.

println("\n" * "="^70)
println("Template complete!")
println("="^70)
println("\nNext steps:")
println("1. Set up your domain and spatial grid")
println("2. Interpolate aridity_raw onto your grid to create aridity_field")
println("3. Complete the uSPACPiParameters by uncommenting aridity_idx line")
println("4. Create the conductance model and use it in your CanopyModel")
println("5. See experiments/calibration/model_interface.jl for full example")
println("\nThe calibrated parameters are now part of the default configuration")
println("and will be automatically used whenever you set up a uSPAC conductance model.")

# Run 2 only: Spatially-varying canopy height
# This script runs only the spatially-varying height simulation
# Assumes Run 1 (default height) has already completed

import SciMLBase
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaParams as CP
using Dates
using ClimaLand
using ClimaLand.Canopy
import ClimaLand.Parameters as LP
using ClimaDiagnostics
using ClimaUtilities.ClimaArtifacts
using ClimaUtilities
import ClimaCore
import ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
const FT = Float64

include(joinpath(pkgdir(ClimaLand), "experiments/integrated/performance/land_spmd.jl"))

@info "Starting Run 2: Spatially-Varying Canopy Height"

# Setup (same as Run 1)
device_suffix = ClimaComms.device() isa ClimaComms.CUDADevice ? "gpu" : "cpu"
context = ClimaComms.context(device_suffix)
start_date = DateTime("2008-03-01")
stop_date = DateTime("2009-03-01")
Δt = Float64(450)
earth_param_set = LP.LandParameters(FT)
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
domain = create_domain_function(FT, context, device_suffix)

# Run 2: Spatially-varying height
root_path_varying = "comparison_varying_height_$(device_suffix)"
diagnostics_outdir_varying = joinpath(root_path_varying, "global_diagnostics")
outdir_varying = ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir_varying)

@info "Using spatially-varying canopy height from CLM data"

# This will use spatially-varying height
function setup_model_varying(FT, start_date, stop_date, Δt, domain, toml_dict)
    # Read CLM canopy height data
    z_top_field = ClimaLand.Canopy.clm_canopy_height(domain.space.surface, FT)
    
    # Cap heights to respect atmospheric forcing constraint
    z_atm = FT(10.0)  # ERA5 atmospheric reference height  
    height_field = ClimaLand.Canopy.effective_canopy_height(z_top_field, z_atm; buffer = FT(2.0))
    
    # Create model using snowy_land_pmodel logic with custom height
    model = setup_model_full(
        FT,
        start_date,
        stop_date,
        Δt,
        domain,
        toml_dict;
        prescribed_biomass_model_args = (; height = height_field)
    )
    
    return model
end

model_varying = setup_model_varying(FT, start_date, stop_date, Δt, domain, toml_dict)
simulation_varying = LandSimulation(start_date, stop_date, Δt, model_varying; outdir = outdir_varying)

CP.log_parameter_information(toml_dict, joinpath(root_path_varying, "parameters.toml"))

@info "Starting simulation with spatially-varying canopy heights..."
@time ClimaLand.Simulations.solve!(simulation_varying)

@info "Run 2 completed! Output saved to: $root_path_varying"

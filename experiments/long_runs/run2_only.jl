# Run 2 ONLY: Spatially-varying canopy height
# Assumes Run 1 is already complete in comparison_default_height_gpu/
# This uses the exact same setup as compare_canopy_heights.jl but only runs Run 2

# Just run the comparison script's Run 2 portion by including it
# and checking if Run 1 output exists to skip it

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

@info "="^50
@info "RUN 2 ONLY: Spatially-Varying Canopy Height"
@info "="^50

# Setup
device_suffix = ClimaComms.device() isa ClimaComms.CUDADevice ? "gpu" : "cpu"
context = ClimaComms.context(device_suffix)
start_date = DateTime("2008-03-01")
stop_date = DateTime("2009-03-01")
Δt = Float64(450)
earth_param_set = LP.LandParameters(FT)
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
domain = create_domain_function(FT, context, device_suffix)

# Output directory for Run 2
root_path = "comparison_varying_height_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir = ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

@info "Using spatially-varying canopy height from CLM data with capping"

# Read CLM canopy height
@info "Reading CLM canopy height data..."
z_top_field = ClimaLand.Canopy.clm_canopy_height(domain.space.surface, FT)

# Cap heights to respect atmospheric forcing constraint
z_atm = FT(10.0)  # ERA5 atmospheric reference height
@info "Applying height capping (z_atm = $z_atm m, buffer = 2 m)..."
height_field = ClimaLand.Canopy.effective_canopy_height(z_top_field, z_atm; buffer = FT(2.0))

@info "Setting up model with spatially-varying canopy height..."
# Use the full model setup with custom height
model = setup_model_full(
    FT,
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict;
    prescribed_biomass_model_args = (; height = height_field)
)

simulation = LandSimulation(start_date, stop_date, Δt, model; outdir = outdir)

CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))

@info "Starting 1-year simulation (2008-03-01 to 2009-03-01)..."
@time ClimaLand.Simulations.solve!(simulation)

@info "="^50
@info "RUN 2 COMPLETED!"
@info "Output saved to: $root_path"
@info "="^50

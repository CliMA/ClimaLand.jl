using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaLand
import ClimaParams
import ClimaLand.Parameters as LP
import ClimaTimeSteppers
using Dates
import ClimaUtilities
import ClimaCore

# First, run a simulation for 2 steps, saving a checkpoint every step.
# Then, restart the simulation after step 1 and complete the simulation up to end of step 2.
# Compare the final states, they should be identical.
FT = Float32
Δt = 900.0
start_date = Dates.DateTime(2008)
t0 = ClimaUtilities.TimeManager.ITime(0, Second(1), start_date)
stop_date = start_date + Second(2*Δt)
root_path = "land_restart"
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(root_path)

domain = ClimaLand.Domains.global_box_domain(FT);
toml_dict = LP.create_toml_dict(FT);
forcing = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = true,
);
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
);
checkpoint_period = Second(Δt);
checkpoint_cb =
    ClimaLand.CheckpointCallback(checkpoint_period, output_dir, t0; model);
simulation = ClimaLand.Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt,
    land,
    diagnostics = nothing,
    user_callbacks = (checkpoint_cb,)
);
ClimaLand.Simulations.solve!(simulation);

# Now, let's restart from the checkpoint (we have the pass the root_path, not the output_dir)
restart_file = ClimaLand.find_restart(root_path)
#closure of restart_file...
set_ic! = make_set_initial_state_from_checkpoint(restart_file, model)

t_restart = ClimaLand.initial_time_from_checkpoint(restart_file; model)
simulation_restart = ;

for p in propertynames(Y.bucket)
    @test getproperty(Y.bucket, p) == getproperty(Y_restart.bucket, p)
end

rm(output_dir; recursive = true)

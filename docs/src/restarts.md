# Restarting Simulations

`ClimaLand` provides functionality to save and load simulation checkpoints,
allowing you to restart simulations from a previous state. This is particularly
useful for long-running simulations or if you want to experiment with different
configurations starting from a specific point in the simulation.


## Saving Checkpoints

To save a simulation checkpoint, you can use the `ClimaLand.save_checkpoint`
function. This function takes the current state `Y`, the simulation time `t`,
and the output directory as arguments. Optionally you can provide the
`ClimaLand` model object model. This will store the hash of the model in the
checkpoint file. You can use this information to ensure that you are restarting
the simulation with the same model that was used to generate the checkpoint.

```julia
ClimaLand.save_checkpoint(Y, t, output_dir; model)
```

Most typically, this function is not called directly. Instead, it is called as a
callback.

In ClimaLand, you can automate the process of saving checkpoints using the
`CheckpointCallback`. This callback allows you to specify the frequency at which
checkpoints are saved and handles the saving process during the simulation.

To use the `CheckpointCallback`, you need to create an instance of it and pass
it to the solve function along with your other callbacks.

Example:

```julia

# ... your ClimaLand simulation setup ...

# Create a CheckpointCallback to save checkpoints every 6 hours
checkpoint_cb = CheckpointCallback(Dates.Hour(6), output_dir, start_date; model, dt)

# Add the callback to the callback set
cb = ClimaTimeSteppers.CallbackSet(checkpoint_cb, other_callbacks...)

# Run the simulation with the callbacks
sol = ClimaTimeSteppers.solve(prob, ode_algo; dt = Δt, callback = cb)

# ... your ClimaLand simulation analysis ...
```

In this example, the `CheckpointCallback` will save a checkpoint every 6 hours
during the simulation. You can customize the checkpoint_frequency to control how
often checkpoints are saved. You can also pass the `ClimaLand` model object model
to store its hash in the checkpoint file. This information can be used later to
ensure that you are restarting the simulation with the same model that was used
to generate the checkpoint.

If `dt` is passed, `CheckpointCallback` will also check that it is consistent
with the checkpoint frequency.


## Restarting from a Checkpoint

To restart a simulation from a checkpoint, you can use the
`ClimaLand.find_restart` function to locate the most recent checkpoint file in
the output directory. Then, you can use the `ClimaLand.read_checkpoint` function
to load the state vector and simulation time from the checkpoint file.

```julia
restart_file = ClimaLand.find_restart(output_dir)
Y, t = ClimaLand.read_checkpoint(restart_file; model)
```

The time is returned as an `ITime` if the checkpoint was saved with one carrying
a start date, and as a number of seconds otherwise. Simulations built from a
`start_date` and a `stop_date` use `ITime`s. Such a simulation is restarted by
continuing from the time of the checkpoint while keeping the epoch of the
original run, which is the convention `ClimaCoupler` uses: the restarted
simulation counts its time from the original start date rather than from the
checkpoint. The final time and the time step have to be built with that epoch
too:

```julia
import ClimaUtilities.TimeManager: ITime, epoch

restart_file = ClimaLand.find_restart(output_dir)
t_restart = ClimaLand.initial_time_from_checkpoint(restart_file; model)
start_date = epoch(t_restart)
t0, tf, Δt = promote(
    t_restart,
    ITime(Dates.value(Second(stop_date - start_date)); epoch = start_date),
    ITime(Δt_seconds; epoch = start_date),
)
simulation = LandSimulation(
    t0,
    tf,
    Δt,
    model;
    set_ic! = (Y, p, t, model) ->
        ClimaLand.set_initial_conditions_from_checkpoint!(Y, restart_file; model),
)
```

## Output Structure

`ClimaLand` utilizes the `OutputPathGenerator` from `ClimaUtilities` to manage
the output directory structure. By default, it uses the `ActiveLinkStyle`, which
creates a sequence of numbered subfolders within the base output directory.

For example, if your base output directory is output, the following structure
will be created:
```
output/
├── output_0000/
│   └── ... checkpoint files ...
├── output_0001/
│   └── ... checkpoint files ...
├── output_0002/
│   └── ... checkpoint files ...
└── output_active -> output_0002/
```

The output_active symbolic link always points to the most recent output
subfolder, making it easy to access the latest simulation results.

### Checkpoint File Structure

When using the `CheckpointCallback`, the checkpoints are saved as HDF5 files
within the numbered output subfolders. The files are named using the following
convention:
```
day<day_number>.<seconds_since_midnight>.hdf5
```
For example, a checkpoint saved at day 10, 3600 seconds after midnight would be
named `day10.3600.hdf5`.

### Checkpoints, drivers, and accumulated diagnostics

At the moment, `ClimaLand` does not support working with accumulated diagnostics
across restarts. The present limitations are best illustrated with an example.

Suppose you are saving 30-day averages and stop the simulation at day 45. If you
do so, you'll find output for day 30 and the checkpoint at day 45. Then, if you
restart the simulation, you'll see that the next diagnostic output will be at
day 75, and not day 60. In other words, the counter starts from 0 with every
restart.

!!! note

    If you care about accurate accumulated diagnostics, make sure to line up your
    checkpoint and diagnostic frequencies.

Similarly, `ClimaLand` does not guarantee drivers values will be set consistently
across restarts. That is, the values of the forcing variables from the last step
of the original run will not coincide with the values of the forcing variables
upon the initial step of the restart run, unless the frequencies of
checkpointing and updating the drivers are compatible.

!!! note

    The reason why `ClimaLand` does not support these features (at the moment), is that
    updating drivers and diagnostics are implemented as callbacks. Callbacks have some
    internal memory that is not saved in the restart files.

### Checkpoints and the cache

Only the state `Y` is checkpointed; the cache is rebuilt from it and from the forcing
when the simulation restarts.

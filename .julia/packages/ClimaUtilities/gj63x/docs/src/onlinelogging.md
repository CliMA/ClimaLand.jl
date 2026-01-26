# Online Simulation Progress Reporting

The `OnlineLogging` module provides tools to monitor and report the simulation
progress and other information of simulations.

Currently, the only feature implemented is timing report, to print information
about current step, simulation time, and average performance. With
`report_walltime`, you will see progress information like the following
```
┌ Info: Progress
│   simulation_time = "2 seconds"
│   n_steps_completed = 20
│   wall_time_per_step = "10 milliseconds, 100 microseconds"
│   wall_time_total = "1 second, 10 milliseconds"
│   wall_time_remaining = "808 milliseconds, 35 microseconds"
│   wall_time_spent = "202 milliseconds, 8 microseconds"
│   percent_complete = "20.0%"
│   estimated_sypd = "0.027"
│   date_now = 2024-12-03T16:01:13.660
└   estimated_finish_date = 2024-12-03T16:01:14.660
```

## `WallTimeInfo` and `report_walltime`

`WallTimeInfo` is a struct that holds information about wall time (the time you
see on your watch, not the simulation time) and that can be used to report
timing information with [`report_walltime`](@ref ClimaUtilities.OnlineLogging.report_walltime).

`WallTimeInfo` keeps track and accumulates how much time has elapsed since the
last time it was updated. In this, `WallTimeInfo` tries to automatically remove
the compilation overhead that your simulation might run into in the first step
(this is accomplished by ignoring the first step and doubling the cost of the
second step to compensate).

The simplest way to use `WallTimeInfo` is to make `report_walltime` a callback.
Here is an example using a `SciML` integrator
```julia
import SciMLBase
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime
# Create the WallTimeInfo
wt = WallTimeInfo()

# Define a schedule that defines how often to report. Here, we follow the signature
# required by SciML. This function is triggered every 10 steps
function every10steps(u, t, integrator)
    return mod(integrator.step, 10) == 0
end

# Next, define the callback
report_callback = SciMLBase.DiscreteCallback(every10steps,
    let wt = wt
        integrator -> report_walltime(wt, integrator)
    end)

# The let wt = wt is not strickly required, but it can improve type-stability and performance

# Now that we have the callback, we can pass it to the SciML constructor for the integrator
```

!!! todo

    Describe schedules when we add them to `ClimaUtilities`

The `report_walltime` function prints various timing statistics:

* **`simulation_time`:** The elapsed time within the simulation.
* **`n_steps_completed`:** The number of simulation steps completed.
* **`wall_time_per_step`:** The average wall clock time per step.
* **`wall_time_total`:** Estimated total wall clock time for the simulation.
* **`wall_time_remaining`:** Estimated remaining wall clock time.
* **`wall_time_spent`:** Total wall clock time already spent.
* **`percent_complete`:** Percentage of the simulation completed.
* **`estimated_sypd`:** Estimated simulated years per day (or days per day if
  the rate is slow).
* **`date_now`:** The current date and time.
* **`estimated_finish_date`:** The estimated date and time of completion.

!!! note

    The estimated values (like `wall_time_remaining` and `estimated_sypd`) are
    based on the average wall time per step and can fluctuate, especially early
    in the simulation.  They become more reliable as the simulation progresses.

## API

```@docs
ClimaUtilities.OnlineLogging.WallTimeInfo
ClimaUtilities.OnlineLogging.report_walltime
```

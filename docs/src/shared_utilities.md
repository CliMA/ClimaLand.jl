# ClimaLand Shared Utilities

## State NaN Counter - `count_nans_state`
We have implemented a function `count_nans_state` which recursively goes
through the entire state and displays the number of NaNs found for each state
variable. This function is intended to be used to debug simulations to
determine quantitatively if a simulation is stable.

If NaNs are found for a particular variable, this will be displayed via
a warning printed to the console. The `verbose` argument toggles whether
the function prints output when no NaNs are found.

If a ClimaCore Field is provided as `mask`, the function will only count NaNs
in the state variables where the mask is 1. This is intended to be used with
the land/sea mask, to avoid counting NaNs over the ocean. Note this assumes
the mask is 1 over land and 0 over ocean.

This function does not distinguish between surface or subsurface
variables, so a variable defined on the subsurface will display more NaNs
than one defined on the surface, even if they are NaN at the same
spatial locations in the horizontal.

### Usage examples
This function can be used to inspect the state after a simulation finishes
running by calling `count_nans_state(sol.u[end])`.

This function can be used throughout the duration of a simulation by
triggering it via a callback. The `NaNCheckCallback` is designed for this
purpose, and can be set up as follows:
```julia
nancheck_interval = Dates.Month(1)
nancheck_cb = ClimaLand.NaNCheckCallback(nancheck_interval, start_date; dt = Δt)
```
and then included along with any other callbacks in a `SciMLBase.CallbackSet`.

Please see our longrun experiments to see examples of this callback in action!

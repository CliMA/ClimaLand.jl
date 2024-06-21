"""This file contains the site-generic time variables for running ClimaLand on
fluxtower sites. These work in tandem with the site-specific timing parameters
found in the {site-ID}_simulation.jl files in each site directory."""

N_spinup_days = 30
N_days = N_spinup_days + 30
tf = Float64(t0 + 3600 * 24 * N_days)
t_spinup = Float64(t0 + N_spinup_days * 3600 * 24)
dt_save = 3600.0
saveat = Array(t_spinup:dt_save:tf)

# Set up timestepper
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

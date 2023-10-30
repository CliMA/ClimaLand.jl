"""This file contains the site-generic time variables for running CliMALSM on 
fluxtower sites. These work in tandem with the site-specific timing parameters 
found in the {site-ID}_simulation.jl files in each site directory."""

N_spinup_days = 30
N_days = N_spinup_days + 30
tf = Float64(t0 + 3600 * 24 * N_days)
t_spinup = Float64(t0 + N_spinup_days * 3600 * 24)
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

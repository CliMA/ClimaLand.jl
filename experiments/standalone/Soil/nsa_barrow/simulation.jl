t0 = FT(150 * 3600 * 24)
N_days_spinup = 0
N_days = N_days_spinup + 60
dt = FT(5)
tf = t0 + FT(3600 * 24 * N_days)
# Number of timesteps between saving output
n = 60
t_spinup = t0 + FT(N_days_spinup * 3600 * 24)
saveat = Array(t0:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

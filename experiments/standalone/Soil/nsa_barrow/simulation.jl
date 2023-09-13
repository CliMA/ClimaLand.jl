t0 = FT(170 * 3600 * 24)
N_days_spinup = 45
N_days = N_days_spinup + 45
dt = FT(60)
tf = t0 + FT(3600 * 24 * N_days)
# Number of timesteps between saving output
n = 60
t_spinup = t0 + FT(N_days_spinup * 3600 * 24)
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

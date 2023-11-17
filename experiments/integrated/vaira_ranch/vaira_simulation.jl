t0 = 21.0 * 3600 * 24
N_days_spinup = 30
N_days = N_days_spinup + 60
dt = 40.0
tf = t0 + 3600 * 24 * N_days

# Number of timesteps between saving output
n = 45
t_spinup = t0 + N_days_spinup * 3600 * 24
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

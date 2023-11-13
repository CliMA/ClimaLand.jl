t0 = Float64(120 * 3600 * 24) # start mid year
N_spinup_days = 10
N_days = N_spinup_days + 30
dt = Float64(120)
tf = Float64(t0 + 3600 * 24 * N_days)

# Number of timesteps between saving output
n = 15
t_spinup = Float64(t0 + N_spinup_days * 3600 * 24)
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

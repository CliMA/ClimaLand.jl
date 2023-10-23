t0 = FT(120 * 3600 * 24)# start mid year
N_days_spinup = 10
N_days = 30 + N_days_spinup
dt = FT(120)
tf = t0 + FT(3600 * 24 * N_days)
# Number of timesteps between saving output
n = 30
t_spinup = t0 + FT(N_days_spinup * 3600 * 24)
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

t0 = FT(120 * 3600 * 24)# start mid year
N_days = 120
dt = FT(120)
tf = t0 + FT(3600 * 24 * N_days)
# Number of timesteps between saving output
n = 30
t_spinup = t0 + FT(10 * 3600 * 24)
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

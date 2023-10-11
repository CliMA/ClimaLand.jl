t0 = FT(0 * 3600 * 24)# start mid year
N_spinup_days = 0
N_days = N_spinup_days + 180
dt = FT(12)
tf = t0 + FT(3600 * 24 * N_days)
# Number of timesteps between saving output
n = 150
t_spinup = t0 + FT(N_spinup_days * 3600 * 24)
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

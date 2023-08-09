t0 = FT(120 * 3600 * 24)# start mid year
N_days = 120
tf = t0 + FT(3600 * 24 * N_days)
dt = FT(60)
n = 60
saveat = Array(t0:(n * dt):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)

N_spinup_days = 30
N_days = N_spinup_days + 30
tf = Float64(t0 + 3600 * 24 * N_days)
t_spinup = Float64(t0 + N_spinup_days * 3600 * 24)
saveat = Array(t_spinup:(n * dt):tf)
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

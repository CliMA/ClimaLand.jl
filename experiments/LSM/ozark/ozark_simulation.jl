t0 = FT(0)
N_days = 365
tf = t0 + FT(3600 * 24 * N_days)
dt = FT(30);
hourly = Array(t0:(3600):(t0 + N_days * 3600 * 24))
timestepper = Euler()

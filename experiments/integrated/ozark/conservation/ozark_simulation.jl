t0 = FT(120 * 3600 * 24)# start mid year
N_days = 10
tf = t0 + FT(3600 * 24 * N_days)
dt = FT(150)
n = 1
saveat = Array(t0:(n * dt):tf)
timestepper = CTS.ARS111()
norm_condition = CTS.MaximumError(FT(1e-8))
conv_checker = CTS.ConvergenceChecker(; norm_condition)
max_iterations = 20
# Set up timestepper
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = max_iterations,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

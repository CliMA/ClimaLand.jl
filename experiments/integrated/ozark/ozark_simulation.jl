
# Column dimensions - separation of layers at the top and bottom of the column:
dz_bottom = FT(1.5)
dz_top = FT(0.025)

# Stem and leaf compartments and their heights:
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9) # m
h_leaf = FT(9.5) # m

# TIME STEPPING:
t0 = Float64(120 * 3600 * 24)# start mid year
N_days = 10
tf = t0 + Float64(3600 * 24 * N_days)
dt = Float64(150)
n = 1
saveat = Array(t0:(n * dt):tf)

timestepper = CTS.ARS111()
# Select conv. condition based on float type due to different precision
err = (FT == Float64) ? 1e-8 : 1e-4
norm_condition = CTS.MaximumError(err)
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

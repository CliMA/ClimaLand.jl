"""
Run the NEON calibration pipeline sequentially with clear step logging.

Steps:
1) set_Station.jl
2) generate_observations.jl
3) run_calibration.jl
"""
using Pkg
Pkg.activate(".buildkite")
Pkg.update()

function run_step(step_name, script_path)
	println("\n=== START: " * step_name * " ===")
	t_start = time()
	try
		include(script_path)
		elapsed = round(time() - t_start; digits = 2)
		println("=== DONE: " * step_name * " (" * string(elapsed) * " s) ===")
	catch err
		elapsed = round(time() - t_start; digits = 2)
		println("=== FAILED: " * step_name * " after " * string(elapsed) * " s ===")
		showerror(stdout, err, catch_backtrace())
		println()
		rethrow()
	end
end

println("Starting NEON calibration pipeline...")

run_step(
	"Set station",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/set_Station.jl",
)
run_step(
	"Generate observations",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl",
)
run_step(
	"Run calibration",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_calibration.jl",
)
run_step(
	"Plot EKI diagnostics",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/plot_eki_diagnostics.jl",
)

ENV["CALL_FIGURES_DIR"] = joinpath(ENV["CALL_OUTPUT_1"], "prior_mean_optimized")
run_step(
	"Run prior mean with optimized parameters",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_prior_mean.jl",
)

println("\nNEON calibration pipeline completed successfully.")
"""
Run the NEON calibration pipeline sequentially with clear step logging.

Steps:
1) set_Station.jl
2) generate_observations.jl
3) run_calibration.jl
"""

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

ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2021-01-01";  ENV["NEON_STOP_DATE"] = "2021-12-31"

run_step(
	"Set station",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/set_Station_mult.jl",
)
run_step(
	"Generate observations",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl",
)
run_step(
	"Run calibration",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_calibration.jl",
)

ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2022-01-01";  ENV["NEON_STOP_DATE"] = "2022-12-31"

run_step(
	"Set station",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/set_Station_mult.jl",
)
run_step(
	"Generate observations",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl",
)
run_step(
	"Run calibration",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_calibration.jl",
)

ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2023-01-01";  ENV["NEON_STOP_DATE"] = "2023-12-31"

run_step(
	"Set station",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/set_Station_mult.jl",
)
run_step(
	"Generate observations",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl",
)
run_step(
	"Run calibration",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_calibration.jl",
)

ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2024-01-01";  ENV["NEON_STOP_DATE"] = "2024-12-31"

run_step(
	"Set station",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/set_Station_mult.jl",
)
run_step(
	"Generate observations",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl",
)
run_step(
	"Run calibration",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_calibration.jl",
)

ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2025-01-01";  ENV["NEON_STOP_DATE"] = "2025-12-31"

run_step(
	"Set station",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/set_Station_mult.jl",
)
run_step(
	"Generate observations",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl",
)
run_step(
	"Run calibration",
	"/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_calibration.jl",
)


println("\nNEON calibration pipeline completed successfully.")
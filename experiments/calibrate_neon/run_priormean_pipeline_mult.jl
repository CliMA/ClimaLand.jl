"""
Run the NEON prior-mean forward model sequentially over multiple date ranges
with clear step logging.

Steps (per date range):
1) run_prior_mean.jl
"""

using Dates


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

println("Starting NEON prior-mean pipeline...")

site_id = "NEON-sjer"
settingsdesc = "NeonExpTo005_main20260602" #name of folder to save results in, e.g. "SPINUPdays_calDepth"

run_prior_script = "/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_prior_mean.jl"
generate_observations_script = "/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl"

date_ranges = [
	("2019-01-02", "2019-12-31"),
	("2020-01-01", "2020-12-31"),
	("2021-01-01", "2021-12-31"),
	("2022-01-01", "2022-12-31"),
	("2023-01-01", "2023-12-31"),
	("2024-01-01", "2024-12-31"),
	("2025-01-01", "2025-12-31"),
	]

ENV["NEON_SPINUP_DAYS"] = "20"
ENV["NEON_N_ITERATIONS"] = "10"
ENV["CALL_DEPTH"] = "0.02" #in X.XXm
ENV["outfolder"] = "output"
ENV["NEON_SITE_ID"] = site_id

# IMPORTANT: do NOT declare `SITE_ID`, `N_ITERATIONS`, `SPINUP_DAYS`, etc. as
# globals or consts here. The included script (`run_prior_mean.jl`) already
# declares them as `const`. If we declare them here as plain globals first,
# Julia raises:
#   "cannot declare Main.X constant; it was already declared global"
# when the included script tries to make them const. So: pass everything via
# ENV only, and let the included script read ENV itself.

for (start_date_str, stop_date_str) in date_ranges
	ENV["NEON_SITE_ID"] = site_id
	ENV["NEON_START_DATE"] = start_date_str
	ENV["NEON_STOP_DATE"] = stop_date_str

	start_date = DateTime(start_date_str)
	stop_date  = DateTime(stop_date_str)
	#replace . with _ in depth for folder naming
	spinup_days_str = get(ENV, "NEON_SPINUP_DAYS", "20")
	n_iter_str      = get(ENV, "NEON_N_ITERATIONS", "10")
	Caldepthnumstr  = replace(string(get(ENV, "CALL_DEPTH", "0.00")), "." => "_")
	Caldepth        = "$(Caldepthnumstr)M"
	ENV["CALL_OUTPUT_PATH"] = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_calibration/$(site_id)/$(site_id)_$(Date(start_date))_$(Date(stop_date))/SpinUP-$(spinup_days_str)d/CalDepth-$(Caldepth)/$(n_iter_str)-It/$(settingsdesc)/"
	println("NEON_SITE_ID: ", ENV["NEON_SITE_ID"])
	println("Outputpath:   ", ENV["CALL_OUTPUT_PATH"])

	run_step("Generate observations", generate_observations_script)
	run_step("Run prior mean", run_prior_script)
end


println("\nNEON prior-mean pipeline completed successfully.")

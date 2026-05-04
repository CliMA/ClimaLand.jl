"""
Run the NEON calibration pipeline sequentially with clear step logging.

Steps:
1) set_Station.jl
2) generate_observations.jl
3) run_calibration.jl
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

println("Starting NEON calibration pipeline...")

site_id = "NEON-moab"
settingsdesc = "NEONextrapReal_newSO2_dt450s" #name of folder to save results in, e.g. "SPINUPdays_calDepth"

generate_observations_script = "/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/generate_observations.jl"
run_calibration_script = "/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/run_calibration.jl"

date_ranges = [
	("2024-01-01", "2024-12-31"),
	("2025-01-01", "2025-12-31"),
	]

ENV["NEON_SPINUP_DAYS"] = "20"
ENV["NEON_N_ITERATIONS"] = "10"
ENV["CALL_DEPTH"] = "0.06" #in X.XXm

const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const N_ITERATIONS = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
const SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))

for (start_date, stop_date) in date_ranges
	ENV["NEON_SITE_ID"] = site_id
	ENV["NEON_START_DATE"] = start_date
	ENV["NEON_STOP_DATE"] = stop_date

	start_date = DateTime(get(ENV, "NEON_START_DATE", string(Date(2009,1,1))))
	stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(2009,12,31))))
	#replace . with _ in depth for folder naming
	Caldepthnumstr = replace(string(get(ENV, "CALL_DEPTH", "0.00")), "." => "_")
	Caldepth = "$(Caldepthnumstr)M"
	ENV["CALL_OUTPUT_PATH"] = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_calibration/$(SITE_ID)/$(SITE_ID)_$(Date(start_date))_$(Date(stop_date))/SpinUP-$(SPINUP_DAYS)d/CalDepth-$(Caldepth)/$(N_ITERATIONS)-It/$(settingsdesc)/"
	println("NEON_SITE_ID: ", ENV["NEON_SITE_ID"])
	
	run_step("Generate observations", generate_observations_script)
	run_step("Run calibration", run_calibration_script)
end


println("\nNEON calibration pipeline completed successfully.")
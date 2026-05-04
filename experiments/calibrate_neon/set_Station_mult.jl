#import Pkg
#Pkg.add("CSV")
using Dates


# Callibration settings
#station and years get set by /home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/start_calibration_pipeline_mult.jl

ENV["NEON_SPINUP_DAYS"] = "20"
ENV["NEON_N_ITERATIONS"] = "10"
ENV["CALL_DEPTH"] = "0.06" #in X.XXm
settingsdesc = "NEONextrapReal_newSO2_dt450s" #name of folder to save results in, e.g. "SPINUPdays_calDepth"
# TimeVaryNoise_
# Output path for both observations and calibration results
const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const N_ITERATIONS = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
const SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
start_date = DateTime(get(ENV, "NEON_START_DATE", string(Date(2009,1,1))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(2009,12,31))))
#replace . with _ in depth for folder naming
Caldepthnumstr = replace(string(get(ENV, "CALL_DEPTH", "0.00")), "." => "_")
Caldepth = "$(Caldepthnumstr)M"
ENV["CALL_OUTPUT_PATH"] = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_calibration/$(SITE_ID)/$(SITE_ID)_$(Date(start_date))_$(Date(stop_date))/SpinUP-$(SPINUP_DAYS)d/CalDepth-$(Caldepth)/$(N_ITERATIONS)-It/$(settingsdesc)/"
#ENV["CALL_OUTPUT_PATH"] = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_calibration/NEON-cper/NEON-cper_2021-01-01_2021-12-31/SpinUP-20d/CalDepth-6cm/10-It/Oldexp_newDAMM4Param_dt180s"
println("NEON_SITE_ID: ", ENV["NEON_SITE_ID"])
#println("NEON_START_DATE: ", ENV["NEON_START_DATE"])
#println("NEON_STOP_DATE: ", ENV["NEON_STOP_DATE"])
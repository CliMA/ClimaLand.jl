#import Pkg
#Pkg.add("CSV")
using Dates


# Callibration settings
#cper 2019
# CPER	2017, 2019
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2017-01-02";  ENV["NEON_STOP_DATE"] = "2017-12-31" 
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2021-01-01";  ENV["NEON_STOP_DATE"] = "2021-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2024-01-01";  ENV["NEON_STOP_DATE"] = "2024-12-31"
# JORN	2019
#ENV["NEON_SITE_ID"] = "NEON-jorn"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
# SJER	2019
#ENV["NEON_SITE_ID"] = "NEON-sjer"; ENV["NEON_START_DATE"] = "2019-01-02";  ENV["NEON_STOP_DATE"] = "2019-12-31"
# SRER	2019
# ENV["NEON_SITE_ID"] = "NEON-srer"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
# STER	2018, 2019
#ENV["NEON_SITE_ID"] = "NEON-ster"; ENV["NEON_START_DATE"] = "2018-01-02";  ENV["NEON_STOP_DATE"] = "2018-12-31"
#ENV["NEON_SITE_ID"] = "NEON-ster"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
# MOAB	2018
#ENV["NEON_SITE_ID"] = "NEON-moab"; ENV["NEON_START_DATE"] = "2018-01-02";  ENV["NEON_STOP_DATE"] = "2018-12-31"

ENV["NEON_SPINUP_DAYS"] = "20"
ENV["NEON_N_ITERATIONS"] = "10"
settingsdesc = "SpinUP-20d_CalDepth-2cm_10-It_NEONextrap_RHupdate" #name of folder to save results in, e.g. "SPINUPdays_calDepth"

# Output path for both observations and calibration results
const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const N_ITERATIONS = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
const SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
start_date = DateTime(get(ENV, "NEON_START_DATE", string(Date(2009,1,1))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(2009,12,31))))

ENV["CALL_OUTPUT_PATH"] = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/$(SITE_ID)_$(Date(start_date))_$(Date(stop_date))/$(settingsdesc)/"

println("NEON_SITE_ID: ", ENV["NEON_SITE_ID"])
#println("NEON_START_DATE: ", ENV["NEON_START_DATE"])
#println("NEON_STOP_DATE: ", ENV["NEON_STOP_DATE"])
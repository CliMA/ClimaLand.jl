#import Pkg
#Pkg.add("CSV")
using Dates


# Callibration settings
#more years
#cper 2019
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2025-12-31" 
#ENV["NEON_SITE_ID"] = "NEON-jorn"; ENV["NEON_START_DATE"] = "2019-01-02";  ENV["NEON_STOP_DATE"] = "2021-12-31"
#ENV["NEON_SITE_ID"] = "NEON-jorn"; ENV["NEON_START_DATE"] = "2020-01-01";  ENV["NEON_STOP_DATE"] = "2021-12-31"
#ENV["NEON_SITE_ID"] = "NEON-jorn"; ENV["NEON_START_DATE"] = "2018-01-02";  ENV["NEON_STOP_DATE"] = "2021-12-31"



# all years
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2017-01-02";  ENV["NEON_STOP_DATE"] = "2025-12-31"
#ENV["NEON_SITE_ID"] = "NEON-moab"; ENV["NEON_START_DATE"] = "2018-01-02";  ENV["NEON_STOP_DATE"] = "2025-12-31"
#ENV["NEON_SITE_ID"] = "NEON-ster"; ENV["NEON_START_DATE"] = "2019-01-02";  ENV["NEON_STOP_DATE"] = "2025-12-31"
#ENV["NEON_SITE_ID"] = "NEON-sjer"; ENV["NEON_START_DATE"] = "2019-01-02";  ENV["NEON_STOP_DATE"] = "2025-12-31"
#ENV["NEON_SITE_ID"] = "NEON-jorn"; ENV["NEON_START_DATE"] = "2018-01-02";  ENV["NEON_STOP_DATE"] = "2025-12-31"
#ENV["NEON_SITE_ID"] = "NEON-srer"; ENV["NEON_START_DATE"] = "2018-01-01";  ENV["NEON_STOP_DATE"] = "2025-12-31"


#single years
# CPER
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2017-01-02";  ENV["NEON_STOP_DATE"] = "2017-12-31" 
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2020-01-01";  ENV["NEON_STOP_DATE"] = "2020-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2021-01-01";  ENV["NEON_STOP_DATE"] = "2021-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2022-01-01";  ENV["NEON_STOP_DATE"] = "2022-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2023-01-01";  ENV["NEON_STOP_DATE"] = "2023-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2024-01-01";  ENV["NEON_STOP_DATE"] = "2024-12-31"
#ENV["NEON_SITE_ID"] = "NEON-cper"; ENV["NEON_START_DATE"] = "2025-01-01";  ENV["NEON_STOP_DATE"] = "2025-12-31"
# JORN	2019
#ENV["NEON_SITE_ID"] = "NEON-jorn"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
ENV["NEON_SITE_ID"] = "NEON-jorn"; ENV["NEON_START_DATE"] = "2021-01-01";  ENV["NEON_STOP_DATE"] = "2021-12-31"

# SJER	2019
#ENV["NEON_SITE_ID"] = "NEON-sjer"; ENV["NEON_START_DATE"] = "2019-01-02";  ENV["NEON_STOP_DATE"] = "2019-12-31"
# SRER	2019
# ENV["NEON_SITE_ID"] = "NEON-srer"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
# STER	2018, 2019
#ENV["NEON_SITE_ID"] = "NEON-ster"; ENV["NEON_START_DATE"] = "2018-01-02";  ENV["NEON_STOP_DATE"] = "2018-12-31"
#ENV["NEON_SITE_ID"] = "NEON-ster"; ENV["NEON_START_DATE"] = "2019-01-01";  ENV["NEON_STOP_DATE"] = "2019-12-31"
# MOAB	2018
#ENV["NEON_SITE_ID"] = "NEON-moab"; ENV["NEON_START_DATE"] = "2018-01-02";  ENV["NEON_STOP_DATE"] = "2018-12-31"
#ENV["NEON_SITE_ID"] = "NEON-moab"; ENV["NEON_START_DATE"] = "2023-01-01";  ENV["NEON_STOP_DATE"] = "2023-12-31"

ENV["NEON_SPINUP_DAYS"] = "20"
ENV["NEON_N_ITERATIONS"] = "5"
ENV["CALL_DEPTH"] = "0.06" #in X.XXm
settingsdesc = "NeonExpTo005_10layers_newSWCInit_main20260504_dt180s"#"NEONextrapReal_newSO2_dt180s" #name of folder to save results in, e.g. "SPINUPdays_calDepth"
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
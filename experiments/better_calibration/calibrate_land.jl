const START_DATE = Dates.DateTime(2008, 12, 01)
const NELEMENTS = (101, 15)
const SHORT_NAMES = ["lhf", "shf", "swu", "lwu"]
const SAMPLE_DATE_RANGES = [(Dates.DateTime(2008, 12, 1), Dates.DateTime(2009, 9, 1))]
const MINIBATCH_SIZE = 1
const SPINUP = Dates.Month(3)

const OUTPUT_DIR = "experiments/better_calibration/land_model/output"
const N_ITERATIONS = 20

# For now, we will use `data_sources.jl` for the leaderboard, since it is
# the easiest option, but it would be better to make your own `data_source.jl`
# and preprocess the observational data to match the simulation data as opposed
# to updating both the simulation and observational data (e.g. with ILAMB data).
include("../long_runs/leaderboard/data_sources.jl")

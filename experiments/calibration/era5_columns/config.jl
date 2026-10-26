# Settings shared by the column selection, the observations, the forward model,
# and the driver of the ERA5 column-ensemble calibration.

import Dates
import EnsembleKalmanProcesses as EKP
import ClimaLand

const FT = Float64

"""
    GlobalConfig

The configuration `experiments/calibration/configs/energy_fluxes.jl` of the
global calibration, whose calibrated variables, noise, prior, and number of
iterations this calibration uses. Its dates, resolution, and minibatches only
apply to the global calibration.
"""
module GlobalConfig
import Dates
import ClimaLand
import EnsembleKalmanProcesses as EKP
const CALIBRATION_DIR =
    joinpath(pkgdir(ClimaLand), "experiments", "calibration")
include(joinpath(CALIBRATION_DIR, "api.jl"))
# Output directory of the global calibration, which the config requires
const OUTPUT_DIR = ""
include(joinpath(CALIBRATION_DIR, "configs", "energy_fluxes.jl"))
end

"""Calibrated diagnostics, in the order they are stacked in the observation
vector."""
const SHORT_NAMES = GlobalConfig.CALIBRATE_CONFIG.short_names

"""Start of the simulations. The first year is spinup and is not compared to
observations."""
const SPINUP_START = Dates.DateTime(2007, 3, 1)

"""Start of the seasonal means compared to observations (the start date of the
default long run)."""
const CALIBRATION_START = Dates.DateTime(2008, 3, 1)

"""End of the simulations and of the seasonal means compared to
observations."""
const STOP_DATE = Dates.DateTime(2010, 3, 1)

const Δt = 900.0

const N_ITERATIONS = GlobalConfig.CALIBRATE_CONFIG.n_iterations

const RNG_SEED = 42

const N_COLUMNS = 200

"""Variables whose observed seasonal climatology the columns are selected by:
the calibrated variables, with GOSIF GPP (`sif_gpp`) in place of FLUXCOM GPP,
which is missing in most deserts and ends in 2013."""
const SELECTION_NAMES = replace(SHORT_NAMES, "gpp" => "sif_gpp")

"""Years of the observed seasonal climatology the columns are selected by: the
years of GOSIF GPP in the `inversion_nee` artifact."""
const CLIMATOLOGY_YEARS = 2002:2020

"""Standard deviation of the observational and structural error of a seasonal
mean at one column, in observation units: the square root of the variances
`NOISE_SCALARS` of the global config. A column with weight `w` (normalized to a
mean of 1) gets the variance `NOISE_STD[name]^2 / w`, so the loss approximates
an area-weighted global misfit."""
const NOISE_STD = Dict(
    name => sqrt(variance) for (name, variance) in GlobalConfig.NOISE_SCALARS
)

"""
    get_calibration_prior()

Return the prior of the calibrated parameters, from the global config. With
`TransformUnscented`, the ensemble has `2p + 1` members for `p` parameters.
"""
get_calibration_prior() = GlobalConfig.get_calibration_prior()

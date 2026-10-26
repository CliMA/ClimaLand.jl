# Settings shared by the column selection, the observations, the forward model,
# and the driver of the ERA5 column-ensemble calibration.

import Dates
import EnsembleKalmanProcesses as EKP

const FT = Float64

"""Calibrated diagnostics, in the order they are stacked in the observation
vector."""
const SHORT_NAMES = ["lhf", "shf", "swu", "lwu", "gpp"]

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

const N_ITERATIONS = 5

const RNG_SEED = 42

const N_COLUMNS = 200

"""Years of the observed seasonal climatology the columns are selected by: the
years of GOSIF GPP in the `inversion_nee` artifact."""
const CLIMATOLOGY_YEARS = 2002:2020

"""Standard deviation of the observational and structural error of a seasonal
mean at one column, in observation units: W m^-2 for the energy fluxes and
g C m^-2 day^-1 for GPP. A column with weight `w` (normalized to a mean of 1)
gets the variance `NOISE_STD[name]^2 / w`, so the loss approximates an
area-weighted global misfit."""
const NOISE_STD = Dict(
    "lhf" => 15.0,
    "shf" => 15.0,
    "swu" => 10.0,
    "lwu" => 10.0,
    "gpp" => 1.5,
)

"""Calibrated factor that multiplies the default PAR and NIR reflectances of
stems, so that the two change together."""
const STEM_REFLECTANCE_SCALE = "stem_reflectance_scale"

"""
    get_calibration_prior()

Return the prior of the calibrated parameters. Except for `pmodel_α` (an
acclimation timescale of two weeks) and `canopy_subcanopy_min_reference_height`
(whose default, 0.5 m, is its lower bound), the priors are centered on the
values in `toml/default_parameters.toml`.

With `TransformUnscented`, the ensemble has `2 * 16 + 1 = 33` members.
"""
function get_calibration_prior()
    priors = [
        # Photosynthesis and its soil moisture limitation: GPP and LHF
        EKP.constrained_gaussian("pmodel_cstar", 0.366, 0.1, 0.05, 1.0),
        EKP.constrained_gaussian("pmodel_β_c3", 24.9, 15.0, 1.0, 300.0),
        EKP.constrained_gaussian("pmodel_β_c4", 8.89, 5.0, 0.5, 100.0),
        EKP.constrained_gaussian("pmodel_α", 1 / 14, 0.03, 0.01, 0.5),
        EKP.constrained_gaussian("moisture_stress_c", 0.591, 0.2, 0.05, 2.0),
        # Stem and root area indices, as multiples of the maximum LAI
        EKP.constrained_gaussian("SAI_coeff", 0.16, 0.08, 0.01, 0.6),
        EKP.constrained_gaussian("RAI_coeff", 1.0, 0.5, 0.1, 5.0),
        # Canopy turbulent exchange and heat storage: LHF and SHF
        EKP.constrained_gaussian("leaf_Cd", 0.07, 0.03, 0.005, 0.5),
        EKP.constrained_gaussian("canopy_z_0m_coeff", 0.0823, 0.04, 0.001, 0.5),
        EKP.constrained_gaussian("canopy_z_0b_coeff", 0.073, 0.04, 0.0001, 0.5),
        EKP.constrained_gaussian("canopy_d_coeff", 0.67, 0.1, 0.3, 0.9),
        EKP.constrained_gaussian("ac_canopy", 2.5e3, 1.0e3, 500.0, 1.0e4),
        # Turbulent exchange of the ground below the canopy
        EKP.constrained_gaussian(
            "canopy_subcanopy_wind_extinction_coefficient",
            0.5,
            0.25,
            0.05,
            2.0,
        ),
        EKP.constrained_gaussian(
            "canopy_subcanopy_min_reference_height",
            1.0,
            0.4,
            0.5,
            3.0,
        ),
        # Radiation: SWU and LWU
        EKP.constrained_gaussian(STEM_REFLECTANCE_SCALE, 1.0, 0.3, 0.2, 2.5),
        EKP.constrained_gaussian("canopy_K_lw", 1.16, 0.3, 0.1, 3.0),
    ]
    return EKP.combine_distributions(priors)
end

"""
    model_parameters(calibrated, default_toml)

Return the model parameters set by the `calibrated` values, a dictionary from
the names of the prior to values. `STEM_REFLECTANCE_SCALE` scales the values of
`stem_PAR_reflectance` and `stem_NIR_reflectance` in `default_toml`; any other
name is a model parameter.
"""
function model_parameters(calibrated, default_toml)
    parameters = Dict{String, Float64}()
    for (name, value) in calibrated
        if name == STEM_REFLECTANCE_SCALE
            for stem_name in ("stem_PAR_reflectance", "stem_NIR_reflectance")
                parameters[stem_name] = value * default_toml[stem_name]
            end
        else
            parameters[name] = value
        end
    end
    return parameters
end

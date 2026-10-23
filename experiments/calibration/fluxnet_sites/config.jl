# Settings shared by the site data, the forward model, and the driver of the
# FLUXNET site calibration.

import Dates
import EnsembleKalmanProcesses as EKP

const FT = Float64

"""Calibrated diagnostics, in the order they are stacked in the observation
vector."""
const SHORT_NAMES = ["lhf", "shf", "swu", "lwu", "gpp"]

"""Simulated years compared to observations. Each site's own years are mapped
onto them (see `source_years`), so all sites share one calendar."""
const CALIBRATION_YEARS = 2008:2010

"""Start of the simulations. The first year is spinup and is not compared to
observations."""
const SPINUP_START = Dates.DateTime(first(CALIBRATION_YEARS) - 1)

"""Start of the monthly means compared to observations."""
const CALIBRATION_START = Dates.DateTime(first(CALIBRATION_YEARS))

"""End of the simulations and of the monthly means compared to observations."""
const STOP_DATE = Dates.DateTime(last(CALIBRATION_YEARS) + 1)

const Δt = 900.0

const N_ITERATIONS = 5

const RNG_SEED = 42

"""Standard deviation of the observational and structural error of a monthly
mean at one site, in observation units: W m^-2 for the energy fluxes and
g C m^-2 day^-1 for GPP. It includes the mismatch between the footprint of a
tower and the column it is compared to."""
const NOISE_STD = Dict(
    "lhf" => 20.0,
    "shf" => 20.0,
    "swu" => 15.0,
    "lwu" => 10.0,
    "gpp" => 2.0,
)

"""Minimum fraction of the half-hours of a calendar year (local time) with
every forcing variable but CO2 observed or gap-filled, for the year to be
simulated."""
const MIN_FORCING_COVERAGE = 0.95

"""Minimum fraction of the half-hours of a month with a measured or good
quality gap-filled flux (quality flag 0 or 1), for its monthly mean to be
compared."""
const MIN_GOOD_FRACTION = 0.5

"""Minimum fraction of the half-hours of a month with a value, for its monthly
mean to be compared."""
const MIN_COVERAGE = 0.9

"""Height (m) of the forcing above the canopy of a column whose tower is not
above it. Where a site has no canopy height, that of the column comes from the
CLM map, which can exceed the height of a tower over short vegetation in a
forested grid cell."""
const MIN_CANOPY_CLEARANCE = 2.0

"""Calibrated factor that multiplies the default PAR and NIR reflectances of
stems, so that the two change together."""
const STEM_REFLECTANCE_SCALE = "stem_reflectance_scale"

"""
    get_calibration_prior()

Return the prior of the calibrated parameters, centered on the values in
`toml/default_parameters.toml`. The canopy roughness coefficients, whose
defaults are orders of magnitude below the usual values, have log-normal
priors.

With `TransformUnscented`, the ensemble has `2 * 18 + 1 = 37` members.
"""
function get_calibration_prior()
    priors = [
        # Photosynthesis and its soil moisture limitation: GPP and LHF
        EKP.constrained_gaussian("pmodel_cstar", 0.463, 0.1, 0.05, 1.0),
        EKP.constrained_gaussian("pmodel_β_c3", 112.0, 50.0, 1.0, 300.0),
        EKP.constrained_gaussian("pmodel_β_c4", 1.76, 1.5, 0.1, 100.0),
        EKP.constrained_gaussian("pmodel_α", 0.0113, 0.008, 0.002, 0.5),
        EKP.constrained_gaussian("moisture_stress_c", 0.759, 0.25, 0.05, 3.0),
        # Stem and root area indices, as multiples of the maximum LAI
        EKP.constrained_gaussian("SAI_coeff", 0.158, 0.08, 0.01, 0.6),
        EKP.constrained_gaussian("RAI_coeff", 0.22, 0.15, 0.02, 5.0),
        # Canopy turbulent exchange and heat storage: LHF and SHF
        EKP.constrained_gaussian("leaf_Cd", 1.17, 0.5, 0.01, 3.0),
        EKP.constrained_gaussian("canopy_z_0m_coeff", 3.26e-6, 1e-4, 0, Inf),
        EKP.constrained_gaussian("canopy_z_0b_coeff", 0.00157, 0.016, 0, Inf),
        EKP.constrained_gaussian("canopy_d_coeff", 0.67, 0.1, 0.3, 0.9),
        EKP.constrained_gaussian("ac_canopy", 2.5e3, 1.0e3, 500.0, 1.0e4),
        # Turbulent exchange of the ground below the canopy
        EKP.constrained_gaussian(
            "canopy_subcanopy_wind_extinction_coefficient",
            0.0556,
            0.05,
            0.005,
            2.0,
        ),
        EKP.constrained_gaussian(
            "canopy_subcanopy_min_reference_height",
            2.86,
            1.0,
            0.5,
            10.0,
        ),
        # Bare soil evaporation and heat exchange
        EKP.constrained_gaussian(
            "maximum_dry_soil_layer_depth",
            0.00296,
            0.003,
            0.0,
            Inf,
        ),
        EKP.constrained_gaussian(
            "soil_scalar_roughness_length",
            0.00302,
            0.003,
            0.0,
            Inf,
        ),
        # Radiation: SWU and LWU
        EKP.constrained_gaussian(STEM_REFLECTANCE_SCALE, 1.0, 0.2, 0.2, 1.25),
        EKP.constrained_gaussian("canopy_K_lw", 0.999, 0.3, 0.1, 3.0),
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

import Dates

import ClimaLand
import ClimaAnalysis
import ClimaAnalysis: OutputVar, NCCatalog
import ClimaUtilities.ClimaArtifacts: @clima_artifact

"""
    preprocess_sim_var(var::OutputVar)

Preprocess the simulation variable `var` based on its short name.

Return the preprocessed `OutputVar`.
"""
function preprocess_sim_var(var::OutputVar)
    return _preprocess_sim_var(var, Val(Symbol(ClimaAnalysis.short_name(var))))
end

function _preprocess_sim_var(var, ::Val{:et})
    (ClimaAnalysis.units(var) == "kg m^-2 s^-1") && (
        var = ClimaAnalysis.convert_units(
            var,
            "mm / day",
            conversion_function = units -> units * 86400.0,
        )
    )
    return var
end

function _preprocess_sim_var(var, ::Val{:gpp})
    # converting from `mol CO2 m^-2 s^-1` in sim to `g C m-2 day-1` in obs
    (ClimaAnalysis.units(var) == "mol CO2 m^-2 s^-1") && (
        var = ClimaAnalysis.convert_units(
            var,
            "g m-2 day-1",
            conversion_function = units -> units * 86400.0 * 12.011,
        )
    )
    return var
end

function _preprocess_sim_var(var, ::Val{:er})
    (ClimaAnalysis.units(var) == "mol CO2 m^-2 s^-1") && (
        var = ClimaAnalysis.convert_units(
            var,
            "g m-2 day-1",
            conversion_function = units -> units * 86400.0 * 12.011,
        )
    )
    return var
end

function _preprocess_sim_var(var, ::Val{:nee})
    (ClimaAnalysis.units(var) == "mol CO2 m^-2 s^-1") && (
        var = ClimaAnalysis.convert_units(
            var,
            "g m-2 day-1",
            conversion_function = units -> units * 86400.0 * 12.011,
        )
    )
    return var
end

function _preprocess_sim_var(var, ::Val{:hr})
    # heterotrophic respiration: convert `mol CO2 m^-2 s^-1` in sim to
    # `g C m-2 day-1` in obs (used for the inv_hr inversion target).
    (ClimaAnalysis.units(var) == "mol CO2 m^-2 s^-1") && (
        var = ClimaAnalysis.convert_units(
            var,
            "g m-2 day-1",
            conversion_function = units -> units * 86400.0 * 12.011,
        )
    )
    return var
end

_preprocess_sim_var(var, ::Val{:lwu}) = var
_preprocess_sim_var(var, ::Val{:lhf}) = var
_preprocess_sim_var(var, ::Val{:shf}) = var
_preprocess_sim_var(var, ::Val{:swu}) = var
# LAI is already dimensionless (m^2 m^-2), matching the MODIS obs; no conversion.
_preprocess_sim_var(var, ::Val{:lai}) = var
_preprocess_sim_var(var, name::Val) =
    error("Preprocessing var with short name ($name) is not defined")

"""
    AbstractDataLoader

Supertype for all data loaders.

`AbstractDataLoader`s have to provide two methods: `available_vars` and
`Base.get`.

`AbstractDataLoader`s load preprocessed data as `OutputVar`s that can be broadly
used for any kind of calibration and making leaderboard plots.
"""
abstract type AbstractDataLoader end

"""
    available_vars(data_loader::AbstractDataLoader)

Return the available preprocessed variables in `data_loader`.
"""
available_vars(data_loader::AbstractDataLoader) = data_loader.available_vars

function Base.show(io::IO, data_loader::AbstractDataLoader)
    vars = sort(collect(data_loader.available_vars))
    printstyled(io, nameof(typeof(data_loader)), bold = true, color = :green)
    print(io, ": ")
    print(io, join(vars, ", "))
end

"""
    preprocess(data_loader::AbstractDataLoader, var, ::Val{varname}) where {varname}

Preprocess `var` with short name `varname` using `data_loader`.

Subtypes of `AbstractDataLoader` must implement this method for each variable
they support. Throws an error if no preprocessing method is defined for
`varname`.
"""
function preprocess(
    data_loader::AbstractDataLoader,
    _,
    ::Val{varname},
) where {varname}
    error(
        "No preprocessing function is found for $varname. Add a method for preprocess(::$(typeof(data_loader)), ::Val{:$varname}) that preprocess this variable as a OutputVar",
    )
end

"""
    ERA5DataLoader

A struct for loading preprocessed ERA5 data as `OutputVar`s.
"""
struct ERA5DataLoader <: AbstractDataLoader
    """A catalog built from NetCDF files, designed to initialize `OutputVar`s.
    See ClimaAnalysis documentation for more information about NCCatalog."""
    catalog::NCCatalog

    """A list of available variables to load."""
    available_vars::Set{String}
end

const ERA5_TO_CLIMA_NAMES = [
    "mslhf" => "lhf",
    "msshf" => "shf",
    "msuwlwrf" => "lwu",
    "msuwswrf" => "swu",
]
const STANDARD_UNITS = Dict("W m**-2" => "W m^-2", "W m-2" => "W m^-2")

"""
    ERA5DataLoader(; era5_to_clima_names = ERA5_TO_CLIMA_NAMES)

Construct a data loader that can be used to load preprocessed ERA5 monthly
time-averaged data in `OutputVar`, where
- the short name, sign of the data, and units match CliMA conventions
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th),
- units match the variables in the output of the CliMA diagnostics.

The ERA5 data comes from the
`era5_monthly_averages_surface_single_level_1979_2024` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/era5_monthly_averages_single_level_1979_2024)
for more information about this artifact.

The keyword argument `era5_to_clima_names` is a vector of pairs mapping
ERA5 name to CliMA name.
"""
function ERA5DataLoader(; era5_to_clima_names = ERA5_TO_CLIMA_NAMES)
    artifact_dir =
        @clima_artifact"era5_monthly_averages_surface_single_level_1979_2024"
    flux_file = joinpath(
        artifact_dir,
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, flux_file, era5_to_clima_names...)
    return ERA5DataLoader(catalog, Set(last.(era5_to_clima_names)))
end

"""
    get(loader::ERA5DataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the ERA5
dataset.
"""
function Base.get(loader::ERA5DataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, add it to ERA5_TO_CLIMA_NAMES as a pair mapping ERA5 name to CliMA name and create a new ERA5DataLoader",
    )
    var = get(catalog, short_name;)
    ClimaAnalysis.transform_dates!(var, Dates.firstdayofmonth)
    return preprocess(loader, var, Val(Symbol(short_name)))
end

"""
    preprocess(::ERA5DataLoader, var, ::Val{varname}) where {varname}

Preprocess `var` with short name `varname`.
"""
preprocess(::ERA5DataLoader, var, ::Val{:shf}) =
    _preprocess_var(var; flip_sign = true)
preprocess(::ERA5DataLoader, var, ::Val{:lhf}) =
    _preprocess_var(var; flip_sign = true)
preprocess(::ERA5DataLoader, var, ::Val{:lwu}) = _preprocess_var(var)
preprocess(::ERA5DataLoader, var, ::Val{:swu}) = _preprocess_var(var)

function _preprocess_var(var; flip_sign = false)
    short_name = ClimaAnalysis.short_name(var)

    flip_sign && replace!(val -> -val, var)
    issorted(ClimaAnalysis.latitudes(var)) ||
        ClimaAnalysis.reverse_dim!(var, ClimaAnalysis.latitude_name(var))

    # Longitudes in other datasets is 0 to 360 degrees, but longitudes in CliMA
    # diagnostics range from -180 to 180 degrees
    if ClimaAnalysis.has_longitude(var)
        var = ClimaAnalysis.shift_longitude(var, -180.0, 180.0)
    end

    var_units = ClimaAnalysis.units(var)
    var_units in keys(STANDARD_UNITS) &&
        (var = ClimaAnalysis.set_units(var, STANDARD_UNITS[var_units]))

    # Functions in ClimaAnalysis can mutate the short name. Calibration uses
    # short names to check the observational and simulation data, so we keep the
    # same short name as before.
    ClimaAnalysis.set_short_name!(var, short_name)
    return var
end

"""
    ILAMBDataLoader

A struct for loading preprocessed ILAMB data as `OutputVar`s.
"""
struct ILAMBDataLoader <: AbstractDataLoader
    """A catalog built from NetCDF files, designed to initialize `OutputVar`s.
    See ClimaAnalysis documentation for more information about NCCatalog."""
    catalog::NCCatalog

    """A list of available variables to load."""
    available_vars::Set{String}
end

"""
    ILAMBDataLoader()

Construct a data loader which you can use to load preprocessed ILAMB monthly
time-averaged data as `OutputVar`s, where
- the short name and units match CliMA conventions,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th),
- units match the variables in the output of the CliMA diagnostics.

The supported variables are `et` (MODIS evapotranspiration), `gpp` (FLUXCOM
gross primary production), `lwu` (CERESed4.2 upwelling longwave radiation), `er`
(FLUXCOM ecosystem respiration), and `nee` (FLUXCOM net ecosystem exchange).
"""
function ILAMBDataLoader()
    et_filepath =
        ClimaLand.Artifacts.ilamb_dataset_path("evspsbl_MODIS_et_0.5x0.5.nc")
    gpp_filepath = ClimaLand.Artifacts.ilamb_dataset_path("gpp_FLUXCOM_gpp.nc")
    lwu_filepath =
        ClimaLand.Artifacts.ilamb_dataset_path("rlus_CERESed4.2_rlus.nc")
    er_filepath = ClimaLand.Artifacts.ilamb_respiration_nee_dataset_path(
        "reco_FLUXCOM_reco.nc",
    )
    nee_filepath = ClimaLand.Artifacts.ilamb_respiration_nee_dataset_path(
        "nee_FLUXCOM_nee.nc",
    )

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, et_filepath, "et" => "et")
    ClimaAnalysis.add_file!(catalog, gpp_filepath, "gpp" => "gpp")
    ClimaAnalysis.add_file!(catalog, lwu_filepath, "rlus" => "lwu")
    ClimaAnalysis.add_file!(catalog, er_filepath, "reco" => "er")
    ClimaAnalysis.add_file!(catalog, nee_filepath, "nee" => "nee")

    return ILAMBDataLoader(catalog, Set(["et", "gpp", "lwu", "er", "nee"]))
end

"""
    get(loader::ILAMBDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the ILAMB
dataset.
"""
function Base.get(loader::ILAMBDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars ||
        error("$short_name is not available to load")
    var = get(catalog, short_name;)
    ClimaAnalysis.transform_dates!(var, Dates.firstdayofmonth)
    return preprocess(loader, var, Val(Symbol(short_name)))
end

"""
    preprocess(::ILAMBDataLoader, var, ::Val{:et})

Preprocess `var` for evapotranspiration from ILAMB MODIS data.
"""
function preprocess(::ILAMBDataLoader, var, ::Val{:et})
    if ClimaAnalysis.units(var) == "kg/m2/s"
        var = ClimaAnalysis.convert_units(
            var,
            "mm / day",
            conversion_function = units -> units * 86400.0,
        )
    end
    replace!(var, missing => NaN)
    return _preprocess_var(var)
end

"""
    preprocess(::ILAMBDataLoader, var, ::Val{:gpp})

Preprocess `var` for gross primary production from ILAMB FLUXCOM data.
"""
function preprocess(::ILAMBDataLoader, var, ::Val{:gpp})
    ClimaAnalysis.dim_units(var, "lon") == "degree" &&
        ClimaAnalysis.set_dim_units!(var, "lon", "degrees_east")
    ClimaAnalysis.dim_units(var, "lat") == "degree" &&
        ClimaAnalysis.set_dim_units!(var, "lat", "degrees_north")
    replace!(var, missing => NaN)
    return _preprocess_var(var)
end

"""
    preprocess(::ILAMBDataLoader, var, ::Val{:lwu})

Preprocess `var` for upwelling longwave radiation from ILAMB CERESed4.2 data.
"""
preprocess(::ILAMBDataLoader, var, ::Val{:lwu}) = _preprocess_var(var)

"""
    preprocess(::ILAMBDataLoader, var, ::Val{:er})

Preprocess `var` for ecosystem respiration from ILAMB FLUXCOM data.
"""
function preprocess(::ILAMBDataLoader, var, ::Val{:er})
    ClimaAnalysis.dim_units(var, "lon") == "degree" &&
        ClimaAnalysis.set_dim_units!(var, "lon", "degrees_east")
    ClimaAnalysis.dim_units(var, "lat") == "degree" &&
        ClimaAnalysis.set_dim_units!(var, "lat", "degrees_north")
    replace!(var, missing => NaN)
    return _preprocess_var(var)
end

"""
    preprocess(::ILAMBDataLoader, var, ::Val{:nee})

Preprocess `var` for net ecosystem exchange from the NOAA CarbonTracker CT2022
product.
"""
function preprocess(::ILAMBDataLoader, var, ::Val{:nee})
    ClimaAnalysis.dim_units(var, "lon") == "degree" &&
        set_dim_units!(var, "lon", "degrees_east")
    ClimaAnalysis.dim_units(var, "lat") == "degree" &&
        set_dim_units!(var, "lat", "degrees_north")
    replace!(var, missing => NaN)
    # nee_FLUXCOM_nee.nc has no _FillValue attribute, so the netCDF
    # default fill (~9.97e36) leaks through as real data — set it to NaN.
    replace!(x -> (!ismissing(x) && abs(x) > 1e20) ? NaN : x, var)
    return _preprocess_var(var)
end

"""
    get_mask_dict(data_loader::ERA5DataLoader)

Return a dictionary mapping short names to a function which takes in `sim_var`,
an `OutputVar` containing simulation data, and `obs_var`, an `OutputVar`
containing observational data, and returns a masking function. The masking
function is used to correctly normalize the global bias and global RMSE.
"""
function get_mask_dict(data_loader::ERA5DataLoader)
    # Dict for loading in masks
    mask_dict = Dict{String, Any}()

    make_mask_fn =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.apply_oceanmask
        end

    mask_dict["shf"] = make_mask_fn
    mask_dict["lhf"] = make_mask_fn
    mask_dict["swu"] = make_mask_fn
    mask_dict["lwu"] = make_mask_fn

    @assert keys(mask_dict) == available_vars(data_loader)
    return mask_dict
end

"""
    get_mask_dict(data_loader::ILAMBDataLoader)

Return a dictionary mapping short names to a function which takes in `sim_var`,
a `OutputVar` containing simulation data, and `obs_var`, a `OutputVar`
containing observational data, and returns a masking function. The masking
function is used to correctly normalize the global bias and global RMSE.
"""
function get_mask_dict(data_loader::ILAMBDataLoader)
    # Dict for loading in masks
    mask_dict = Dict{String, Any}()

    mask_dict["lwu"] =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.apply_oceanmask
        end

    make_mask_fn =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.make_lonlat_mask(
                ClimaAnalysis.slice(
                    obs_var,
                    time = ClimaAnalysis.times(obs_var) |> first,
                );
                set_to_val = isnan,
            )
        end

    mask_dict["et"] = make_mask_fn
    mask_dict["gpp"] = make_mask_fn
    mask_dict["er"] = make_mask_fn
    mask_dict["nee"] = make_mask_fn

    @assert keys(mask_dict) == available_vars(data_loader)
    return mask_dict
end

"""
    get_compare_vars_biases_plot_extrema(; annual = false)

Return a dictionary mapping short names to ranges for the bias plots.

If `annual = true`, the limits are tightened since the errors in annual means
are smaller.

To add a variable to the leaderboard, add a key-value pair to the dictionary
`compare_vars_biases_plot_extrema` whose key is the short name of the variable
and whose value is a tuple `(lower, upper)` setting the color scale bounds for
the bias plots.
"""
function get_compare_vars_biases_plot_extrema(; annual = false)
    factor = annual ? 1 / sqrt(2) : 1.0 # treats each season as ~ independent
    compare_vars_biases_plot_extrema = Dict(
        "et" => (-2.0, 2.0) .* factor,
        "gpp" => (-6.0, 6.0) .* factor,
        "er" => (-6.0, 6.0) .* factor,
        "nee" => (-4.0, 4.0) .* factor,
        "hr" => (-4.0, 4.0) .* factor,
        "lwu" => (-40.0, 40.0) .* factor,
        "shf" => (-50.0, 50.0) .* factor,
        "lhf" => (-40.0, 40.0) .* factor,
        "swu" => (-50.0, 50.0) .* factor,
        "lai" => (-3.0, 3.0) .* factor,
    )
    return compare_vars_biases_plot_extrema
end

"""
    _monthly_total_to_daily_rate!(obs_var)

In-place: convert `obs_var.data` from `g C m⁻² month⁻¹` to a daily rate by
dividing each time slice by its real days-in-month (28–31), rather than the
constant 365.25/12 ≈ 30.4375. Avoids the ~5% source of error the constant
days-per-month assumption introduces. Sets the units string to `"g m-2 day-1"`.
"""
function _monthly_total_to_daily_rate!(obs_var)
    native = ClimaAnalysis.units(obs_var)
    native == "g C m-2 month-1" || error(
        "expected a monthly total in \"g C m-2 month-1\", got \"$native\"",
    )
    for (i, date) in enumerate(ClimaAnalysis.dates(obs_var))
        slice = ClimaAnalysis.view_select(
            obs_var;
            by = ClimaAnalysis.Index(),
            time = i,
        )
        slice.data ./= Dates.daysinmonth(date)
    end
    ClimaAnalysis.set_units!(obs_var, "g m-2 day-1")
    return obs_var
end

"""
    get_inversion_obs_var_dict()

Inversion-derived NEE/GPP/ER/Rh observations from the `inversion_nee` artifact
(`derived_nee_gpp_er_rh_2002_2020.nc`, monthly 1°×1°, 2002–2020), returned as
full-timeseries `OutputVar`s (the framework sets the per-sample reference date
and windows by season in `generate_observations.jl`).

Variables (all positive = source except gpp = uptake):
  - `nee` from CarbonTracker CT2022, g C m⁻² month⁻¹
  - `gpp` from GOSIF-GPP v2,        g C m⁻² month⁻¹
  - `er`  residual = nee + gpp,     g C m⁻² month⁻¹
  - `rh`  from Hashimoto 2015,      **g C m⁻² day⁻¹** (native units)

All four are returned in g C m⁻² day⁻¹ (nee/gpp/er divided by actual
days-in-month; rh already daily), with the units *string* set to `"g m-2 day-1"`
to match the model diagnostics — ClimaCalibrate's `UnitsChecker` does a raw
string compare. Short_names are retagged to `inv_nee`/`sif_gpp`/`res_er`/`inv_hr`
so they don't collide with the FLUXCOM `gpp`/`er`/`nee` keys.
"""
function get_inversion_obs_var_dict()
    inversion_path = ClimaLand.Artifacts.inversion_nee_dataset_path()

    _load(varname) = begin
        obs_var = ClimaAnalysis.OutputVar(inversion_path, varname)
        ClimaAnalysis.dim_units(obs_var, "lon") == "degree" &&
        ClimaAnalysis.set_dim_units!(obs_var, "lon", "degrees_east")
        ClimaAnalysis.dim_units(obs_var, "lat") == "degree" &&
        ClimaAnalysis.set_dim_units!(obs_var, "lat", "degrees_north")
        ClimaAnalysis.transform_dates!(obs_var, Dates.firstdayofmonth)
        return ClimaAnalysis.replace(obs_var, missing => NaN)
    end

    obs_var_dict = Dict{String, Any}()
    for (alias, varname, monthly) in (
        ("inv_nee", "nee", true),
        ("sif_gpp", "gpp", true),
        ("res_er", "er", true),
        ("inv_hr", "rh", false), # rh already a daily rate; relabeled below
    )
        obs_var = _load(varname)
        if monthly
            # Convert to a daily rate and set units to "g m-2 day-1".
            _monthly_total_to_daily_rate!(obs_var)
        else
            # rh is already a daily rate; relabel "g C m-2 day-1" → "g m-2 day-1"
            # so the UnitsChecker's raw string compare matches the model hr
            # diagnostic (a mismatch NaN-fills the member's whole G column).
            native = ClimaAnalysis.units(obs_var)
            native == "g C m-2 day-1" || error(
                "inversion `$varname` expected units \"g C m-2 day-1\", got \"$native\"",
            )
            ClimaAnalysis.set_units!(obs_var, "g m-2 day-1")
        end
        # Like ILAMB obs: sort lat ascending and shift lon to [-180, 180] so
        # preprocess_single_obs_var aligns these with the model grid.
        # `_preprocess_var` leaves "g m-2 day-1" and the short_name untouched.
        obs_var = _preprocess_var(obs_var)
        obs_var.attributes["short_name"] = alias
        obs_var_dict[alias] = obs_var
    end
    return obs_var_dict
end

# Map the calibration aliases (artifact-derived short names) to the model
# diagnostic short names used by the simulation output and the leaderboard.
# `inv_hr` (Hashimoto Rh) is intentionally excluded: the leaderboard shows the
# MODIS `lai` target in its place (see FlagshipCarbonMetricsDataLoader).
const _INVERSION_ALIAS_TO_MODEL_NAME =
    Dict("inv_nee" => "nee", "sif_gpp" => "gpp", "res_er" => "er")

"""
    FlagshipCarbonMetricsDataLoader

Loads our flagship carbon-cycle observations for the leaderboard: CT2022 NEE,
GOSIF GPP, and residual ER from the `inversion_nee` artifact, plus the MODIS
`lai` target. Like `ILAMBDataLoader`, it serves monthly obs keyed by the model
short names `nee`/`gpp`/`er`/`lai`, but sourced from these products instead of
FLUXCOM.
"""
struct FlagshipCarbonMetricsDataLoader <: AbstractDataLoader
    """Preprocessed inversion `OutputVar`s, keyed by model short name."""
    obs_var_dict::Dict{String, Any}

    """A list of available variables to load."""
    available_vars::Set{String}
end

"""
    FlagshipCarbonMetricsDataLoader()

Construct a data loader for the inversion-derived carbon targets and MODIS LAI.
The carbon variables are already preprocessed (monthly total → daily rate,
latitude sorted ascending, longitude shifted to [-180, 180], units set to
`g m-2 day-1`) by `get_inversion_obs_var_dict`; here they are re-keyed and
re-tagged from `inv_nee`/`sif_gpp`/`res_er` to `nee`/`gpp`/`er`. The `lai` target
(m^2 m^-2) is loaded from `get_modis_lai_obs_var`.
"""
function FlagshipCarbonMetricsDataLoader()
    inversion_dict = get_inversion_obs_var_dict()
    obs_var_dict = Dict{String, Any}()
    for (alias, model_name) in _INVERSION_ALIAS_TO_MODEL_NAME
        obs_var = inversion_dict[alias]
        obs_var.attributes["short_name"] = model_name
        obs_var_dict[model_name] = obs_var
    end
    lai_obs = get_modis_lai_obs_var()
    lai_obs.attributes["short_name"] = "lai"
    obs_var_dict["lai"] = lai_obs
    return FlagshipCarbonMetricsDataLoader(
        obs_var_dict,
        Set(keys(obs_var_dict)),
    )
end

"""
    get(loader::FlagshipCarbonMetricsDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the model short name `short_name` (one of
`nee`, `gpp`, `er`, `lai`).
"""
function Base.get(loader::FlagshipCarbonMetricsDataLoader, short_name::String)
    short_name in loader.available_vars ||
        error("$short_name is not available to load")
    return loader.obs_var_dict[short_name]
end

"""
    get_mask_dict(data_loader::FlagshipCarbonMetricsDataLoader)

Return a dictionary mapping model short names to a masking function used to
normalize the global bias and RMSE. The inversion carbon variables
(`nee`/`gpp`/`er`) are masked where the observation is missing (NaN), mirroring
the ILAMB carbon-variable masks. `lai` uses `ClimaAnalysis.apply_oceanmask`
instead: MODIS LAI is finite (≈0), not NaN, over ocean, so the `isnan` mask
would mask nothing.
"""
function get_mask_dict(data_loader::FlagshipCarbonMetricsDataLoader)
    mask_dict = Dict{String, Any}()

    make_mask_fn =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.make_lonlat_mask(
                ClimaAnalysis.slice(
                    obs_var,
                    time = ClimaAnalysis.times(obs_var) |> first,
                );
                set_to_val = isnan,
            )
        end

    for short_name in available_vars(data_loader)
        mask_dict[short_name] =
            short_name == "lai" ?
            ((sim_var, obs_var) -> ClimaAnalysis.apply_oceanmask) : make_mask_fn
    end

    @assert keys(mask_dict) == available_vars(data_loader)
    return mask_dict
end

"""
    get_modis_lai_obs_var(; years = 2000:2020)

MODIS LAI (Yuan et al., monthly 1°×1°, per-year files) as a single `OutputVar`
spanning `years`, keyed `lai`, units `m^2 m^-2` (native, matching the model
`lai` diagnostic). Latitude is sorted ascending and longitude shifted to
[-180, 180] via `_preprocess_var`.

The native samples are ~30.4-day spaced and calendar-unaligned, so they are
linearly interpolated in time onto calendar month-starts to match the
month-start-dated model `lai` diagnostic. Only interior month-starts are kept
(endpoints would need extrapolation), so `years` should bracket the window of
interest by at least one month on each side.
"""
function get_modis_lai_obs_var(; years = 2000:2020)
    paths = [
        ClimaLand.Artifacts.modis_lai_single_year_path(; year) for year in years
    ]
    obs_var = ClimaAnalysis.OutputVar(paths, "lai")
    obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

    # Resample onto calendar month-starts (interior only) by linear time
    # interpolation, so seasonal averaging produces the same canonical
    # month-start season dates as the model's monthly `lai` diagnostic.
    start_date = Dates.DateTime(obs_var.attributes["start_date"])
    data_dates = ClimaAnalysis.dates(obs_var)
    month_starts = filter(
        d -> first(data_dates) <= d <= last(data_dates),
        Dates.firstdayofmonth(first(data_dates)):Dates.Month(1):last(
            data_dates,
        ),
    )
    target_seconds = [
        Float64(Dates.value(Dates.Second(d - start_date))) for d in month_starts
    ]
    obs_var = ClimaAnalysis.resampled_as(obs_var; time = target_seconds)

    obs_var = _preprocess_var(obs_var)
    obs_var.attributes["short_name"] = "lai"
    obs_var.attributes["units"] = "m^2 m^-2"
    return obs_var
end

"""
    get_calibration_obs_var_dict(; short_names = nothing)

Return a dictionary mapping short names to `OutputVar` containing preprocessed
observational data for calibration. This combines ERA5 energy flux variables
(`lhf`, `shf`, `lwu`, `swu`), ILAMB variables (`gpp`, `er`, `nee`), the
inversion-derived carbon targets (`inv_nee`, `sif_gpp`, `res_er`, `inv_hr`)
from the `inversion_nee` artifact, and the MODIS `lai` target.

If `short_names` is provided, only the requested variables are returned (and
the MODIS LAI file load is skipped unless `lai` is requested).
"""
function get_calibration_obs_var_dict(; short_names = nothing)
    obs_var_dict = Dict{String, Any}()
    # ERA5 energy fluxes
    era5_data_loader = ERA5DataLoader()
    # ILAMB GPP, ecosystem respiration (FLUXCOM reco), and NEE
    ilamb_data_loader = ILAMBDataLoader()

    for short_name in available_vars(era5_data_loader)
        obs_var_dict[short_name] = get(era5_data_loader, short_name)
    end

    # Add more variables to obs_var_dict
    for short_name in ("gpp", "er", "nee")
        obs_var_dict[short_name] = get(ilamb_data_loader, short_name)
    end

    # Inversion-derived carbon targets (CT2022 NEE + GOSIF GPP + residual ER +
    # Hashimoto Rh), keyed by inv_nee/sif_gpp/res_er/inv_hr.
    merge!(obs_var_dict, get_inversion_obs_var_dict())

    # MODIS LAI target for optimal-LAI calibration. Loading the 21 per-year
    # files is comparatively expensive, so only build it when requested.
    if isnothing(short_names) || "lai" in short_names
        obs_var_dict["lai"] = get_modis_lai_obs_var()
    end

    if !isnothing(short_names)
        filter!(p -> p.first in short_names, obs_var_dict)
    end

    return obs_var_dict
end

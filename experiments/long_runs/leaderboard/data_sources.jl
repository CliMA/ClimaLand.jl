import ClimaAnalysis
import ClimaLand

"""
    get_sim_var_dict(diagnostics_folder_path)

Return a dictionary mapping short names to `OutputVar` containing preprocessed
simulation data. This is used by the function `compute_leaderboard`.

To add a variable for the leaderboard, add a key-value pair to the dictionary
`sim_var_dict` whose key is the short name of the variable and the value is an
anonymous function that returns a `OutputVar`. For each variable, any
preprocessing should be done in the corresponding anonymous function which
includes unit conversion and shifting the dates.

The variable should have only three dimensions: time, longitude, and latitude.
"""
function get_sim_var_dict(diagnostics_folder_path)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()

    sim_var_dict["lwu"] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "lwu",
            )
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end

    sim_var_dict["et"] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "et",
            )
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            (ClimaAnalysis.units(sim_var) == "kg m^-2 s^-1") && (
                sim_var = ClimaAnalysis.convert_units(
                    sim_var,
                    "mm / day",
                    conversion_function = units -> units * 86400.0,
                )
            )
            return sim_var
        end


    sim_var_dict["gpp"] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "gpp",
            )
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            # converting from to `mol CO2 m^-2 s^-1` in sim to `g C m-2 day-1` in obs
            (ClimaAnalysis.units(sim_var) == "mol CO2 m^-2 s^-1") && (
                sim_var = ClimaAnalysis.convert_units(
                    sim_var,
                    "g m-2 day-1",
                    conversion_function = units -> units * 86400.0 * 12.011,
                )
            )
            return sim_var
        end

    sim_var_dict["lhf"] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "lhf",
            )
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end

    sim_var_dict["shf"] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "shf",
            )
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end

    sim_var_dict["swu"] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "swu",
            )
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end
    return sim_var_dict
end

"""
    get_obs_var_dict(data_source)

Return a dictionary mapping short names to `OutputVar` containing preprocessed
observational data. This is used by the function `compute_leaderboard`.The
argument `data_source` can either be "ERA5" or "ILAMB".

To add a variable for the leaderboard, add a key-value pair to the dictionary
`obs_var_dict` whose key is the short name of the variable and the value is an
anonymous function that returns a `OutputVar`. The function must take in a
start date which is used to align the times in the observational data to match
the simulation data. The short name must be the same as in `sim_var_dict` in the
function `sim_var_dict`. For each variable, any preprocessing is done in the
corresponding anonymous function which includes unit conversion and shifting the
dates.

The variable should have only three dimensions: latitude, longitude, and time.
"""
function get_obs_var_dict(data_source)
    if uppercase(data_source) == "ILAMB"
        return get_ilamb_obs_var_dict()
    elseif uppercase(data_source) == "ERA5"
        return get_era5_obs_var_dict()
    else
        error("Did not expect $data_source for the data source")
    end
end

"""
    get_ilamb_obs_var_dict()

Get ILAMB observational data. See `get_obs_var_dict` for more information for
how to add new variables.
"""
function get_ilamb_obs_var_dict()
    # Dict for loading in observational data
    obs_var_dict = Dict{String, Any}()
    obs_var_dict["et"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                ClimaLand.Artifacts.ilamb_dataset_path(
                    "evspsbl_MODIS_et_0.5x0.5.nc",
                ),
                "et",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            (ClimaAnalysis.units(obs_var) == "kg/m2/s") && (
                obs_var = ClimaAnalysis.convert_units(
                    obs_var,
                    "mm / day",
                    conversion_function = units -> units * 86400.0,
                )
            )
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)
            return obs_var
        end

    obs_var_dict["gpp"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                ClimaLand.Artifacts.ilamb_dataset_path("gpp_FLUXCOM_gpp.nc"),
                "gpp",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            ClimaAnalysis.dim_units(obs_var, "lon") == "degree" &&
                (obs_var.dim_attributes["lon"]["units"] = "degrees_east")
            ClimaAnalysis.dim_units(obs_var, "lat") == "degree" &&
                (obs_var.dim_attributes["lat"]["units"] = "degrees_north")
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)
            return obs_var
        end

    obs_var_dict["lwu"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                ClimaLand.Artifacts.ilamb_dataset_path(
                    "rlus_CERESed4.2_rlus.nc",
                ),
                "rlus",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            ClimaAnalysis.units(obs_var) == "W m-2" &&
                (obs_var = ClimaAnalysis.set_units(obs_var, "W m^-2"))
            return obs_var
        end
    return obs_var_dict
end

"""
    get_era5_obs_var_dict()

Get ERA5 observational data. See `get_obs_var_dict` for more information for
how to add new variables.
"""
function get_era5_obs_var_dict()
    obs_var_dict = Dict{String, Any}()
    era5_data_path = joinpath(
        ClimaLand.Artifacts.era5_monthly_averages_single_level_path(),
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )

    obs_var_dict["lhf"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "mslhf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            (ClimaAnalysis.units(obs_var) == "W m**-2") && (
                obs_var = ClimaAnalysis.convert_units(
                    obs_var,
                    "W m^-2",
                    conversion_function = units -> units * -1.0,
                )
            )
            return obs_var
        end

    obs_var_dict["shf"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "msshf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )

            (ClimaAnalysis.units(obs_var) == "W m**-2") && (
                obs_var = ClimaAnalysis.convert_units(
                    obs_var,
                    "W m^-2",
                    conversion_function = units -> units * -1.0,
                )
            )

            return obs_var
        end

    obs_var_dict["lwu"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "msuwlwrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            (ClimaAnalysis.units(obs_var) == "W m**-2") &&
                (obs_var = ClimaAnalysis.set_units(obs_var, "W m^-2"))
            return obs_var
        end

    obs_var_dict["swu"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "msuwswrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )

            (ClimaAnalysis.units(obs_var) == "W m**-2") &&
                (obs_var = ClimaAnalysis.set_units(obs_var, "W m^-2"))
            return obs_var
        end
    return obs_var_dict
end

"""
    get_mask_dict()

Return a dictionary mapping short names to a function which takes in `sim_var`,
a `OutputVar` containing simulation data, and `obs_var`, a `OutputVar`
containing observational data, and return a masking function. The argument
`data_source` is either `ILAMB` or `ERA5`.

To add a variable to the leaderboard, add a key-value pair to the dictionary
`mask_dict` whose key is the same short name in `sim_var_dict` and the value is
a function that takes in a `OutputVar` representing simulation data and a
`OutputVar` representing observational data and returns a masking function or
`nothing` if a masking function is not needed. The masking function is used to
correctly normalize the global bias and global RMSE.
"""
function get_mask_dict(data_source)
    # Dict for loading in masks
    mask_dict = Dict{String, Any}()

    # LWU is found in both ERA5 and ILAMB
    mask_dict["lwu"] =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.apply_oceanmask
        end

    if uppercase(data_source) == "ILAMB"
        mask_dict["et"] =
            (sim_var, obs_var) -> begin
                return ClimaAnalysis.make_lonlat_mask(
                    ClimaAnalysis.slice(
                        obs_var,
                        time = ClimaAnalysis.times(obs_var) |> first,
                    );
                    set_to_val = isnan,
                )
            end

        mask_dict["gpp"] =
            (sim_var, obs_var) -> begin
                return ClimaAnalysis.make_lonlat_mask(
                    ClimaAnalysis.slice(
                        obs_var,
                        time = ClimaAnalysis.times(obs_var) |> first,
                    );
                    set_to_val = isnan,
                )
            end
    elseif uppercase(data_source) == "ERA5"
        mask_dict["shf"] =
            (sim_var, obs_var) -> begin
                return ClimaAnalysis.apply_oceanmask
            end

        mask_dict["lhf"] =
            (sim_var, obs_var) -> begin
                return ClimaAnalysis.apply_oceanmask
            end

        mask_dict["swu"] =
            (sim_var, obs_var) -> begin
                return ClimaAnalysis.apply_oceanmask
            end
    else
        error("Did not expect $data_source for the data source")
    end
    return mask_dict
end

"""
    get_compare_vars_biases_plot_extrema()

Return a dictionary mapping short names to ranges for the bias plots.

To add a variable to the leaderboard, add a key-value pair to the dictionary
`compare_vars_biases_plot_extrema` whose key is a short name key is the same
short name in `sim_var_pfull_dict` in the function `get_sim_var_pfull_dict` and
the value is a tuple, where the first element is the lower bound and the last
element is the upper bound for the bias plots.
"""
function get_compare_vars_biases_plot_extrema()
    compare_vars_biases_plot_extrema = Dict(
        "et" => (-2.0, 2.0),
        "gpp" => (-6.0, 6.0),
        "lwu" => (-40.0, 40.0),
        "shf" => (-50.0, 50.0),
        "lhf" => (-40.0, 40.0),
        "swu" => (-50.0, 50.0),
    )
    return compare_vars_biases_plot_extrema
end

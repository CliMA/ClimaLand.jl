"""
    get_sim_var_dict()

Retrieve the following variables from the simulation output directory:
- lhf
- shf
- lwu
- swu
- lwd
- swd
- lwn (computed as lwd - lwu)
- swn (computed as swd - swu)
"""
function get_sim_var_dict(simdir)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()

    # Get LHF by converting from ET
    toml_dict = LP.create_toml_dict(Float64)
    earth_param_set = LP.LandParameters(toml_dict)
    # _LH_v0 = LP.LH_v0(earth_param_set) # J/kg
    sim_var_dict["lhf"] =
        () -> begin
            sim_var_lhf = get(simdir, short_name = "lhf") # units (W/m²)
            sim_var_lhf.attributes["long_name"] = "Latent heat flux"
            sim_var_lhf.attributes["units"] = "W m⁻²"

            sim_var_lhf =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var_lhf)

            return sim_var_lhf
        end

    # Read in SHF
    sim_var_dict["shf"] =
        () -> begin
            sim_var_shf = get(simdir, short_name = "shf") # units (W/m²)
            sim_var_shf.attributes["long_name"] = "Sensible heat flux"
            sim_var_shf.attributes["units"] = "W m⁻²"

            sim_var_shf =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var_shf)
            return sim_var_shf
        end

    # Read in LWU
    sim_var_dict["lwu"] =
        () -> begin
            sim_var_lwu = get(simdir, short_name = "lwu") # units (W/m²)
            sim_var_lwu.attributes["long_name"] = "Upward longwave radiation"
            sim_var_lwu.attributes["units"] = "W m⁻²"

            sim_var_lwu =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var_lwu)
            return sim_var_lwu
        end

    # Read in SWU
    sim_var_dict["swu"] =
        () -> begin
            sim_var_swu = get(simdir, short_name = "swu") # units (W/m²)
            sim_var_swu.attributes["long_name"] = "Upward shortwave radiation"
            sim_var_swu.attributes["units"] = "W m⁻²"

            sim_var_swu =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var_swu)
            return sim_var_swu
        end

    # Read in LWD
    sim_var_dict["lwd"] =
        () -> begin
            sim_var_lwd = get(simdir, short_name = "lwd") # units (W/m²)
            sim_var_lwd.attributes["long_name"] = "Downward longwave radiation"
            sim_var_lwd.attributes["units"] = "W m⁻²"

            sim_var_lwd =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var_lwd)
            return sim_var_lwd
        end

    # Read in SWD
    sim_var_dict["swd"] =
        () -> begin
            sim_var_swd = get(simdir, short_name = "swd") # units (W/m²)
            sim_var_swd.attributes["long_name"] = "Downward shortwave radiation"
            sim_var_swd.attributes["units"] = "W m⁻²"

            sim_var_swd =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var_swd)
            return sim_var_swd
        end

    # Compute LWN = LWD - LWU
    sim_var_dict["lwn"] =
        () -> begin
            lwd = sim_var_dict["lwd"]()
            lwu = sim_var_dict["lwu"]()
            lwn = lwd - lwu
            lwn.attributes["long_name"] = "Net longwave radiation"
            lwn.attributes["units"] = "W m⁻²"
            lwn.attributes["start_date"] = lwu.attributes["start_date"]
            return lwn
        end

    # Compute SWN = SWD - SWU
    sim_var_dict["swn"] =
        () -> begin
            swd = sim_var_dict["swd"]()
            swu = sim_var_dict["swu"]()
            swn = swd - swu
            swn.attributes["long_name"] = "Net shortwave radiation"
            swn.attributes["units"] = "W m⁻²"
            swn.attributes["start_date"] = swu.attributes["start_date"]
            return swn
        end

    return sim_var_dict
end

"""
    get_obs_var_dict()

Retrieve the following variables from ERA5 observational data:
- lhf
- shf
- lwu
- swu
- lwd
- swd
- lwn (computed as lwd - lwu)
- swn (computed as swd - swu)

Note: `get_obs_var_dict` in ext/land_sim_vis/leaderboard/data_sources.jl uses
the artifact `era5_monthly_averages_single_level_path`, which does not include lwd and swd.

This function does the same, but also gets lwd and swd from `era5_land_forcing_data2008_folder_path`.
"""
function get_obs_var_dict(
    comparison_start_date,
    comparison_end_date;
    is_local = false,
)
    # contains monthly mslhf, msshf, msuwlwrf, msuwswrf
    era5_data_path = joinpath(
        ClimaLand.Artifacts.era5_monthly_averages_single_level_path(),
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    ) #lwu
    # contains hourly msdwlwrf, msdwswrf
    if is_local
        era5_land_forcing_data_path =
            ClimaLand.Artifacts.era5_land_forcing_data2008_lowres_path()
    else
        era5_land_forcing_data_path = ClimaLand.Artifacts.find_era5_year_paths(
            comparison_start_date, # start_date
            comparison_end_date, # stop_date
        ) #lwd
    end

    # Dict for loading in observational data
    obs_var_dict = Dict{String, Any}()
    obs_var_dict["lhf"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "mslhf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m² = J/m²s)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            (ClimaAnalysis.units(obs_var) == "W m**-2") && (
                obs_var = ClimaAnalysis.convert_units(
                    obs_var,
                    "W m⁻²",
                    conversion_function = units -> units * -1.0,
                )
            )

            obs_var.attributes["short_name"] = "lhf"
            obs_var.attributes["start_date"] = start_date
            obs_var.attributes["long_name"] = "Latent heat flux"

            return obs_var
        end

    # Read in SHF
    obs_var_dict["shf"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "msshf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m² = J/m²s)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)


            (ClimaAnalysis.units(obs_var) == "W m**-2") && (
                obs_var = ClimaAnalysis.convert_units(
                    obs_var,
                    "W m⁻²",
                    conversion_function = units -> units * -1.0,
                )
            )

            obs_var.attributes["short_name"] = "shf"
            obs_var.attributes["start_date"] = start_date
            obs_var.attributes["long_name"] = "Sensible heat flux"

            return obs_var
        end

    # Read in LWU
    obs_var_dict["lwu"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "msuwlwrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m²)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            obs_var.attributes["short_name"] = "lwu"
            obs_var.attributes["start_date"] = start_date
            obs_var.attributes["long_name"] = "Upward longwave radiation"
            obs_var.attributes["units"] = "W m⁻²"

            return obs_var
        end

    # Read in SWU
    obs_var_dict["swu"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_data_path,
                "msuwswrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m²)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            obs_var.attributes["short_name"] = "swu"
            obs_var.attributes["start_date"] = start_date
            obs_var.attributes["long_name"] = "Upward shortwave radiation"
            obs_var.attributes["units"] = "W m⁻²"

            return obs_var
        end

    # Read in LWD
    obs_var_dict["lwd"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_land_forcing_data_path,
                "msdwlwrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m²)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            obs_var.attributes["short_name"] = "lwd"
            obs_var.attributes["start_date"] = start_date
            obs_var.attributes["long_name"] = "Downward longwave radiation"
            obs_var.attributes["units"] = "W m⁻²"

            return obs_var
        end

    # Read in SWD
    obs_var_dict["swd"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                era5_land_forcing_data_path,
                "msdwswrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m²)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            obs_var.attributes["short_name"] = "swd"
            obs_var.attributes["start_date"] = start_date
            obs_var.attributes["long_name"] = "Downward shortwave radiation"
            obs_var.attributes["units"] = "W m⁻²"

            return obs_var
        end


    # Compute LWN = LWD - LWU
    obs_var_dict["lwn"] =
        (start_date) -> begin
            # LWD is very large because it's hourly, so it takes a while to load
            lwd = obs_var_dict["lwd"](start_date)
            lwu = obs_var_dict["lwu"](start_date)

            monthly_times = lwu.dims["time"]
            lwd_monthly_data = similar(lwu.data)
            # average hourly LWD data to monthly so we can compare with monthly LWU
            for (idx, time) in enumerate(monthly_times)
                # collect all hourly time indices in this month
                # Note both datasets are saved to the first day of the month by ClimaAnalysis
                time_idxs = findall(x -> (x == time), lwd.dims["time"])

                # average hourly data to monthly
                lwd_monthly_data[:, :, idx] =
                    mean(lwd.data[:, :, time_idxs], dims = 3)
            end

            # Compute LWN = LWD - LWU
            lwn_data = lwd_monthly_data - lwu.data

            attribs = Dict(
                "short_name" => "lwn",
                "long_name" => "Net longwave radiation",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )

            lwn = ClimaAnalysis.OutputVar(
                attribs,
                lwu.dims,
                lwu.dim_attributes,
                lwn_data,
            )
            return lwn
        end

    # Compute SWN = SWD - SWU
    obs_var_dict["swn"] =
        (start_date) -> begin
            # SWD is very large because it's hourly, so it takes a while to load
            swd = obs_var_dict["swd"](start_date)
            swu = obs_var_dict["swu"](start_date)

            monthly_times = swu.dims["time"]
            swd_monthly_data = similar(swu.data)
            # average hourly SWD data to monthly so we can compare with monthly SWU
            for (idx, time) in enumerate(monthly_times)
                # collect all hourly time indices in this month
                time_idxs = findall(x -> (x == time), swd.dims["time"])

                # average hourly data to monthly
                swd_monthly_data[:, :, idx] =
                    mean(swd.data[:, :, time_idxs], dims = 3)
            end

            # Compute SWN = SWD - SWU
            swn_data = swd_monthly_data - swu.data

            attribs = Dict(
                "short_name" => "swn",
                "long_name" => "Net shortwave radiation",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )

            swn = ClimaAnalysis.OutputVar(
                attribs,
                swu.dims,
                swu.dim_attributes,
                swn_data,
            )
            return swn
        end

    return obs_var_dict
end

using ClimaAnalysis

"""
    get_sim_var_dict(simdir)

Returns a dictionary of functions that load in simulation data from `simdir`.

These include the following variables:
- lhf: latent heat flux (W/m²)
- shf: sensible heat flux (W/m²)
- lwu: upward longwave radiation (W/m²)
- swu: upward shortwave radiation (W/m²)
- lwd: downward longwave radiation (W/m²)
- swd: downward shortwave radiation (W/m²)
- lwn: net longwave radiation (W/m²)
- swn: net shortwave radiation (W/m²)
"""
function get_sim_var_dict(simdir)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()

    # Read in LHF
    sim_var_dict["lhf"] =
        () -> begin
            sim_var_lhf = get(simdir, short_name = "lhf") # units (W/m²)
            sim_var_lhf.attributes["long_name"] = "Latent heat flux"
            sim_var_lhf.attributes["units"] = "W m⁻²"
            return sim_var_lhf
        end

    # Read in SHF
    sim_var_dict["shf"] =
        () -> begin
            sim_var_shf = get(simdir, short_name = "shf") # units (W/m²)
            sim_var_shf.attributes["long_name"] = "Sensible heat flux"
            sim_var_shf.attributes["units"] = "W m⁻²"
            return sim_var_shf
        end

    # Read in LWU
    sim_var_dict["lwu"] =
        () -> begin
            sim_var_lwu = get(simdir, short_name = "lwu") # units (W/m²)
            sim_var_lwu.attributes["long_name"] = "Upward longwave radiation"
            sim_var_lwu.attributes["units"] = "W m⁻²"
            return sim_var_lwu
        end

    # Read in SWU
    sim_var_dict["swu"] =
        () -> begin
            sim_var_swu = get(simdir, short_name = "swu") # units (W/m²)
            sim_var_swu.attributes["long_name"] = "Upward shortwave radiation"
            sim_var_swu.attributes["units"] = "W m⁻²"
            return sim_var_swu
        end

    # Read in LWD
    sim_var_dict["lwd"] =
        () -> begin
            sim_var_lwd = get(simdir, short_name = "lwd") # units (W/m²)
            sim_var_lwd.attributes["long_name"] = "Downward longwave radiation"
            sim_var_lwd.attributes["units"] = "W m⁻²"
            return sim_var_lwd
        end

    # Read in SWD
    sim_var_dict["swd"] =
        () -> begin
            sim_var_swd = get(simdir, short_name = "swd") # units (W/m²)
            sim_var_swd.attributes["long_name"] = "Upward shortwave radiation"
            sim_var_swd.attributes["units"] = "W m⁻²"
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

Returns a dictionary of functions that load in observational data.
All data is in hourly intervals.

The following variables come from the monthly ERA5 dataset (`era5_monthly_averages_2008`):
- lhf: latent heat flux (W/m²)
- shf: sensible heat flux (W/m²)
- lwu: upward longwave radiation (W/m²)
- swu: upward shortwave radiation (W/m²)

The following variables come from the hourly ERA5 dataset (`era5_land_forcing_data2008`),
and are averaged to monthly intervals:
- lwd: downward longwave radiation (W/m²)
- swd: downward shortwave radiation (W/m²)

From these, we compute:
- lwn: net longwave radiation (W/m²)
- swn: net shortwave radiation (W/m²)
"""
function get_obs_var_dict()
    # contains monthly mslhf, msshf, msuwlwrf, msuwswrf
    era5_data_path = joinpath(
        ClimaLand.Artifacts.era5_monthly_averages_2008_folder_path(),
        "era5_monthly_surface_fluxes_200801-200812.nc",
    )
    # contains hourly msdwlwrf, msdwswrf (these take a long time to load)
    era5_land_forcing_data_path = joinpath(
        ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(
            lowres = false,
        ),
        "era5_2008_1.0x1.0.nc",
    )

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

            # Add attributes and make values positive to account for convention difference
            attribs = Dict(
                "short_name" => "lhf",
                "long_name" => "Latent heat flux",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )

            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims[ClimaAnalysis.longitude_name(obs_var)], 360)
            new_data =
                -1 * cat(obs_var.data, obs_var.data[[1], :, :], dims = 1)
            obs_var_pos = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                new_data,
            )
            return obs_var_pos
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

            # Add attributes and make values positive to account for convention difference
            attribs = Dict(
                "short_name" => "shf",
                "long_name" => "Sensible heat flux",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )
            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims[ClimaAnalysis.longitude_name(obs_var)], 360)
            new_data =
                -1 * cat(obs_var.data, obs_var.data[[1], :, :], dims = 1)
            obs_var_pos = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                new_data,
            )

            return obs_var_pos
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

            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims[ClimaAnalysis.longitude_name(obs_var)], 360)
            new_data = cat(obs_var.data, obs_var.data[[1], :, :], dims = 1)

            # Manually set attributes
            attribs = Dict(
                "short_name" => "lwu",
                "long_name" => "Upward longwave radiation",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )

            obs_var_shifted = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                new_data,
            )
            return obs_var_shifted
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

            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims[ClimaAnalysis.longitude_name(obs_var)], 360)
            new_data = cat(obs_var.data, obs_var.data[[1], :, :], dims = 1)
            # Manually set attributes
            attribs = Dict(
                "short_name" => "swu",
                "long_name" => "Upward shortwave radiation",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )

            obs_var_shifted = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                new_data,
            )
            return obs_var_shifted
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

            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims[ClimaAnalysis.longitude_name(obs_var)], 360)
            new_data = cat(obs_var.data, obs_var.data[[1], :, :], dims = 1)
            # Manually set attributes
            attribs = Dict(
                "short_name" => "lwd",
                "long_name" => "Downward longwave radiation",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )

            obs_var_shifted = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                new_data,
            )
            return obs_var_shifted
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

            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims[ClimaAnalysis.longitude_name(obs_var)], 360)
            new_data = cat(obs_var.data, obs_var.data[[1], :, :], dims = 1)
            # Manually set attributes
            attribs = Dict(
                "short_name" => "swd",
                "long_name" => "Downward shortwave radiation",
                "units" => "W m⁻²",
                "start_date" => start_date,
            )

            obs_var_shifted = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                new_data,
            )
            return obs_var_shifted
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

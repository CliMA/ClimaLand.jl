
function get_sim_var_dict(simdir)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()

    # Get LHF by converting from ET
    earth_param_set = LP.LandParameters(Float64)
    # _LH_v0 = LP.LH_v0(earth_param_set) # J/kg
    sim_var_dict["lhf"] =
        () -> begin
            sim_var_lhf = get(simdir, short_name = "lhf") # units (W/m²)
            sim_var_lhf.attributes["long_name"] = "Latent heat flux"
            sim_var_lhf.attributes["units"] = "W m⁻²"

            # sim_var_et = get(simdir, short_name = "et") # units (kg/m²s)

            # attribs = Dict(
            #     "short_name" => "lhf",
            #     "long_name" => "Latent heat flux",
            #     "units" => "W m⁻²",
            #     "start_date" => sim_var_et.attributes["start_date"],
            # )
            # sim_var_lhf = ClimaAnalysis.OutputVar(
            #     attribs,
            #     sim_var_et.dims,
            #     sim_var_et.dim_attributes,
            #     sim_var_et.data * _LH_v0, # W/m²
            # )
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

    # # Compute LWN = LWD - LWU
    # sim_var_dict["lwn"] =
    #     () -> begin
    #         lwd = sim_var_dict["lwd"]()
    #         lwu = sim_var_dict["lwu"]()
    #         lwn = lwd - lwu
    #         lwn.attributes["long_name"] = "Net longwave radiation"
    #         lwn.attributes["units"] = "W m⁻²"
    #         return lwn
    #     end

    # # Compute SWN = SWD - SWU
    # sim_var_dict["swn"] =
    #     () -> begin
    #         swd = sim_var_dict["swd"]()
    #         swu = sim_var_dict["swu"]()
    #         swn = swd - swu
    #         swn.attributes["long_name"] = "Net shortwave radiation"
    #         swn.attributes["units"] = "W m⁻²"
    #         return swn
    #     end

    return sim_var_dict
end

function get_obs_var_dict()
    era5_data_path = joinpath(
        ClimaLand.Artifacts.era5_monthly_averages_2008_folder_path(),
        "era5_monthly_surface_fluxes_200801-200812.nc",
    )
    # era5_land_forcing_data_path = joinpath(
    #     ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(
    #         lowres = false,
    #     ),
    #     "era5_2008_0.9x1.25.nc",
    # )

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
            push!(obs_var.dims["longitude"], 360)
            new_data =
                -1 * cat(obs_var.data, obs_var.data[[1], :, :], dims = 1)
            # # relabel longitude values from [0, 359] to [-180, 179]
            # obs_var.dims["longitude"] = obs_var.dims["longitude"] .- 180
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
            push!(obs_var.dims["longitude"], 360)
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
            push!(obs_var.dims["longitude"], 360)
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
            push!(obs_var.dims["longitude"], 360)
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
                era5_data_path,
                "msdwlwrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m²)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims["longitude"], 360)
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
                era5_data_path,
                "msdwswrf",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            ) # units (W/m²)
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)

            # Copy data at lon = 0 to lon = 360 to avoid white lines
            push!(obs_var.dims["longitude"], 360)
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


    # # Compute LWN = LWD - LWU
    # obs_var_dict["lwn"] =
    #     (start_date) -> begin
    #         lwd = obs_var_dict["lwd"](start_date)
    #         lwu = obs_var_dict["lwu"](start_date)
    #         lwn = lwd - lwu
    #         lwn.attributes["long_name"] = "Net longwave radiation"
    #         lwn.attributes["units"] = "W m⁻²"
    #         return lwu
    #     end

    # # Compute SWN = SWD - SWU
    # obs_var_dict["swn"] =
    #     (start_date) -> begin
    #         swd = obs_var_dict["swd"](start_date)
    #         swu = obs_var_dict["swu"](start_date)
    #         swn = swd - swu
    #         swn.attributes["long_name"] = "Net shortwave radiation"
    #         swn.attributes["units"] = "W m⁻²"
    #         return swn
    #     end

    return obs_var_dict
end

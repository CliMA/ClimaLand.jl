
function get_sim_var_dict(simdir)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()

    # Get LHF by converting from ET
    earth_param_set = LP.LandParameters(Float64)
    _LH_v0 = LP.LH_v0(earth_param_set) # J/kg
    sim_var_dict["lhf"] =
        () -> begin
            sim_var_et = get(simdir, short_name = "et") # units (kg/m²s)

            attribs = Dict(
                "short_name" => "lhf",
                "long_name" => "Latent heat flux",
                "units" => "W/m²",
                "start_date" => sim_var_et.attributes["start_date"],
            )
            sim_var_lhf = ClimaAnalysis.OutputVar(
                attribs,
                sim_var_et.dims,
                sim_var_et.dim_attributes,
                sim_var_et.data * _LH_v0, # W/m²
            )
            return sim_var_lhf
        end

    # Read in SHF
    sim_var_dict["shf"] =
        () -> begin
            sim_var_shf = get(simdir, short_name = "shf") # units (W/m²)
            sim_var_shf.attributes["long_name"] = "Sensible heat flux"
            sim_var_shf.attributes["units"] = "W/m²"
            return sim_var_shf
        end

    # Read in LWU
    sim_var_dict["lwu"] =
        () -> begin
            sim_var_lwu = get(simdir, short_name = "lwu") # units (W/m²)
            sim_var_lwu.attributes["long_name"] = "Upward longwave radiation"
            sim_var_lwu.attributes["units"] = "W/m²"
            return sim_var_lwu
        end

    # Read in SWU
    sim_var_dict["swu"] =
        () -> begin
            sim_var_swu = get(simdir, short_name = "swu") # units (W/m²)
            sim_var_swu.attributes["long_name"] = "Upward shortwave radiation"
            sim_var_swu.attributes["units"] = "W/m²"
            return sim_var_swu
        end

    # Read in LWD (not to plot, but for energy balance)
    sim_var_dict["lwd"] =
        () -> begin
            sim_var_lwd = get(simdir, short_name = "lwd") # units (W/m²)
            sim_var_lwd.attributes["long_name"] = "Downward longwave radiation"
            sim_var_lwd.attributes["units"] = "W/m²"
            return sim_var_lwd
        end

    # Read in SWD
    sim_var_dict["swd"] =
        () -> begin
            sim_var_swd = get(simdir, short_name = "swd") # units (W/m²)
            sim_var_swd.attributes["long_name"] = "Upward shortwave radiation"
            sim_var_swd.attributes["units"] = "W/m²"
            return sim_var_swd
        end

    return sim_var_dict
end

function get_obs_var_dict()
    era5_data_path = joinpath(
        ClimaLand.Artifacts.era5_monthly_averages_2008_folder_path(),
        "era5_monthly_surface_fluxes_200801-200812.nc",
    )
    # era5_monthly_surface_fluxes_200801-200812.nc
    # era5_data_path = joinpath(ClimaLand.Artifacts.era5_monthly_averages_2008_folder_path(lowres = true), "era5_2008_1.0x1.0_lowres.nc")
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
                "units" => "W/m²",
                "start_date" => start_date,
            )
            obs_var_pos = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                -1 * obs_var.data,
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
                "units" => "W/m²",
                "start_date" => start_date,
            )
            obs_var_pos = ClimaAnalysis.OutputVar(
                attribs,
                obs_var.dims,
                obs_var.dim_attributes,
                -1 * obs_var.data,
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

            # Manually set attributes
            obs_var.attributes["short_name"] = "lwu"
            obs_var.attributes["long_name"] = "Upward longwave radiation"
            obs_var.attributes["units"] = "W/m²"
            obs_var.attributes["start_date"] = start_date
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

            # Manually set attributes
            obs_var.attributes["short_name"] = "swu"
            obs_var.attributes["long_name"] = "Upward shortwave radiation"
            obs_var.attributes["units"] = "W/m²"
            obs_var.attributes["start_date"] = start_date
            return obs_var
        end

    return obs_var_dict
end

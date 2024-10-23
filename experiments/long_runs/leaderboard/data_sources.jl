import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import GeoMakie
import CairoMakie

# For plotting bias plots
compare_vars_biases_plot_extrema =
    Dict("et" => (-3.0, 3.0), "gpp" => (-8.0, 8.0), "lwu" => (-40.0, 40.0))

# Dict for loading in simulation data
sim_var_dict = Dict{String,Any}()

for short_name in ["et", "lwu"]
    sim_var_dict[short_name] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = short_name,
            )
            # Remove the line below later when start_date is added to the diagnostics
            haskey(sim_var.attributes, "start_date") || (sim_var.attributes["start_date"] = "2012-01-01T00:00:00")
            sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            sim_var = ClimaAnalysis.shift_by_years(sim_var, 9) # TODO: REMOVE LATER
            return sim_var
        end
end

for short_name in ["gpp"]
    sim_var_dict[short_name] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = short_name,
            )
            # Remove the line below later when start_date is added to the diagnostics
            haskey(sim_var.attributes, "start_date") || (sim_var.attributes["start_date"] = "2012-01-01T00:00:00")
            sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            sim_var = ClimaAnalysis.shift_by_years(sim_var, 9) # TODO: REMOVE LATER
            sim_var = ClimaAnalysis.convert_units(
                sim_var,
                "g m-2 day-1",
                conversion_function = units -> units * 86400.0 * 12.011,
            )
            return sim_var
        end
end

# Dict for loading in observational data
obs_var_dict = Dict{String,Any}()
obs_var_dict["et"] =
    (start_date) -> begin
        obs_var = ClimaAnalysis.OutputVar(
            "obs_data/evspsbl/GLEAMv3.3a/et.nc",
            "et",
            new_start_date = start_date,
            shift_by = Dates.firstdayofmonth,
        )
        ClimaAnalysis.units(obs_var) == "kg/m2/s" &&
            (obs_var = ClimaAnalysis.set_units(obs_var, "kg m^-2 s^-1"))
        obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)
        return obs_var
    end

obs_var_dict["gpp"] =
    (start_date) -> begin
        obs_var = ClimaAnalysis.OutputVar(
            "obs_data/gpp/FLUXCOMM/gpp_reverse_dims_flux_comm.nc",
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
            "obs_data/rlus/CERESed4.2/rlus.nc",
            "rlus",
            new_start_date = start_date,
            shift_by = Dates.firstdayofmonth,
        )
        ClimaAnalysis.units(obs_var) == "W m-2" &&
            (obs_var = ClimaAnalysis.set_units(obs_var, "W m^-2"))
        return obs_var
    end

# Dict for loading in masks
mask_dict = Dict{String,Any}()

mask_dict["et"] =
    (sim_var, obs_var) -> begin
        return ClimaAnalysis.Visualize.oceanmask(),
        ClimaAnalysis.make_lonlat_mask(
            ClimaAnalysis.slice(obs_var, time = ClimaAnalysis.times(obs_var) |> first);
            set_to_zero = isnan,
        )
    end

mask_dict["gpp"] =
    (sim_var, obs_var) -> begin
        return ClimaAnalysis.Visualize.oceanmask(),
        ClimaAnalysis.make_lonlat_mask(
            ClimaAnalysis.slice(obs_var, time = ClimaAnalysis.times(obs_var) |> first);
            set_to_zero = isnan,
        )
    end

mask_dict["lwu"] = (sim_var, obs_var) -> begin
    return nothing, nothing
end

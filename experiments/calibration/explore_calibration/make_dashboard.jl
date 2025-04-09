# TODO:
# Currently, on load_button.value is a bit bugged (sometime works, sometimes not).
# works most of the time though. especially if all menu to 1st value.

import EnsembleKalmanProcesses as EKP
import GeoMakie as GM
using WGLMakie
using Bonito
using JLD2
using Statistics
using Printf

# These will be in EKP get_metadata in the future
include("make_training_locations.jl")
nelements = (50, 10)
locations = make_training_locations(nelements)

########### Some useful functions ###########################

function RMSE(x, z)
    return sqrt(mean((x.-z).^2))
end

function error_abs(x, z)
    return mean(abs.(x.-z))
end

"""
    slicevar(v, n, n_vars)

Return the slice for variable `n` (1-based) out of `n_vars` total variables,
from the flattened vector `v`.
"""
function slicevar(v, n, n_vars)
    # Calculate the stride (how many positions to move for the next site/location)
    stride = n_vars * 4  # 4 seasons per location per variable

    # Calculate indices for variable n
    start_offset = (n - 1) * 4  # Starting offset for this variable

    # Create indices list for this variable (handling all seasons for all locations)
    idx = vcat([(i + start_offset):(i + start_offset + 3) for i in 1:stride:length(v)]...)

    return v[idx]
end


function add_seasonal_access(data_dict, vars)
    result = Dict{Any, Any}()

    for (site_id, site_data) in data_dict
        result[site_id] = Dict{String, Any}()

        for var_key in vars
            # Create a new dictionary for each variable with seasonal slices
            if haskey(site_data, var_key)
                values = site_data[var_key]
                result[site_id][var_key] = Dict{String, Vector{Float64}}(
                                                                         "DJF" => values[1:4:end],
                                                                         "MAM" => values[2:4:end],
                                                                         "JJA" => values[3:4:end],
                                                                         "SON" => values[4:4:end]
                                                                        )
            end
        end
    end

    return result
end

function add_seasonal_access_g(data_dict, vars)
    result = Dict{Any, Any}()

    for (outer_key, outer_val) in data_dict
        result[outer_key] = Dict{Any, Any}()

        for (middle_key, middle_val) in outer_val
            result[outer_key][middle_key] = Dict{String, Any}()

            for var_key in vars
                if haskey(middle_val, var_key)
                    values = middle_val[var_key]
                    # Create seasonal slices for each variable
                    result[outer_key][middle_key][var_key] = Dict{String, Vector{Float64}}(
                                                                                           "DJF" => values[1:4:end],
                                                                                           "MAM" => values[2:4:end],
                                                                                           "JJA" => values[3:4:end],
                                                                                           "SON" => values[4:4:end]
                                                                                          )
                end
            end
        end
    end

    return result
end

# Calculate mean lat weight averaged
function cosine_weighted_global_mean(values::Vector{Float64}, lats::Vector{Float64})
    @assert length(values) == length(lats) "Length mismatch between values and latitudes"

    weights = cosd.(lats)  # cosd for degrees
    numerator = sum(values .* weights)
    denominator = sum(weights)

    return numerator / denominator
end

# Function to load data and process it
function load_and_process_data(eki_file, prior_file, variable_file, locations_file="all_locations.jld2")
    println("Loading: EKI=$(eki_file), Prior=$(prior_file)")

    # Load eki object, locations, prior
    eki = JLD2.load_object(eki_file)
    prior = include(prior_file)
    variable_list_file = include(variable_file)
#    @load locations_file locations


    # Get basic information
    n_ensembles = EKP.get_N_ens(eki)
    n_iterations = EKP.get_N_iterations(eki)
    errors = eki.error
    normalized_errors = errors ./ errors[1] .* 100

    # Get all g
    g_all = [EKP.get_g(eki, i) for i in 1:n_iterations]

    # Get all y
    obs_series = EKP.get_observation_series(eki)
    y_obs = obs_series.observations
    y_all = [EKP.get_obs(y_obs[i]) for i in 1:n_iterations]

    # Get all constrained parameters
    params = EKP.get_ϕ(prior, eki)
    param_dict = Dict(i => [params[i][:, j] for j in 1:size(params[i], 2)] for i in eachindex(params))
    params_name = prior.name

    # Get variable list - extract from the prior or config if available
    # This is a bit of an assumption - adjust based on how your prior stores variable names
    n_vars = length(variable_list)

    # Process y_data
    y_data = Dict()
    for iteration_n in 1:length(y_all)
        y_data[iteration_n] = Dict()
        for (var_idx, var_name) in enumerate(variable_list)
            y_data[iteration_n][var_name] = slicevar(y_all[iteration_n], var_idx, n_vars)
        end
    end

    # Process g_data
    g_data = Dict()
    for iteration_n in 1:length(g_all)
        g = g_all[iteration_n]
        g_data[iteration_n] = Dict()

        # Create indices for each variable
        idxs = Dict()
        for (var_idx, var_name) in enumerate(variable_list)
            idxs[var_name] = slicevar(1:size(g, 1), var_idx, n_vars)
        end

        # Fill in g_data
        for ensemble in 1:size(g, 2)
            g_data[iteration_n][ensemble] = Dict(
                                                 var => g[idxs[var], ensemble] for var in variable_list
                                                )
        end
    end

    # Create structures with seasonal access
    seasonal_y_data = add_seasonal_access(y_data, variable_list)
    seasonal_g_data = add_seasonal_access_g(g_data, variable_list)

    # Extract locations
    lons = map(x -> x[1], locations)
    lats = map(x -> x[2], locations)

    # Define RMSE benchmarks (can be customized based on your variables)
    rmse_benchmarks = Dict()
    for var in variable_list
        # Default benchmarks, can be customized
        if var == "lhf"
            rmse_benchmarks[var] = 20
        elseif var == "shf"
            rmse_benchmarks[var] = 15
        elseif var == "swu"
            rmse_benchmarks[var] = 25
        elseif var == "lwu"
            rmse_benchmarks[var] = 10
        else
            # Default for other variables
            rmse_benchmarks[var] = 20
        end
    end

    return Dict(
                "eki" => eki,
                "prior" => prior,
                "n_ensembles" => n_ensembles,
                "n_iterations" => n_iterations,
                "errors" => errors,
                "normalized_errors" => normalized_errors,
                "g_all" => g_all,
                "y_all" => y_all,
                "params" => params,
                "param_dict" => param_dict,
                "params_name" => params_name,
                "variable_list" => variable_list,
                "y_data" => y_data,
                "g_data" => g_data,
                "seasonal_y_data" => seasonal_y_data,
                "seasonal_g_data" => seasonal_g_data,
                "lons" => lons,
                "lats" => lats,
                "rmse_benchmarks" => rmse_benchmarks
               )
end

function update_fig(load_button, menu_calibration, menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, ax_anomalies, ax_sm, seasonal_g_data, seasonal_y_data, lons, lats)

        # Use @lift to properly handle observables
    m_v = menu_var.value
    m_i = menu_iter.value
    m_m = menu_m.value
    m_s = menu_season.value

    g = @lift(seasonal_g_data[$m_i][$m_m][$m_v][$m_s])
    y = @lift(seasonal_y_data[$m_i][$m_v][$m_s])
    anomalies = @lift($g .- $y)
    rmse_y_g = @lift(string("RMSE = ", round(RMSE($g, $y), digits=1), " W m⁻²"))

    min_p = @lift(minimum(vcat($g, $y)))
    max_p = @lift(maximum(vcat($g, $y)))
    limits_p = @lift(($min_p, $max_p))

    min_ano = @lift(minimum($anomalies))
    max_ano = @lift(maximum($anomalies))
    limits_ano = (-30, 30)

    p_g = heatmap!(ax_g, lons, lats, g, colorrange = limits_p)
    p_y = heatmap!(ax_y, lons, lats, y, colorrange = limits_p)
    p_ano = heatmap!(ax_anomalies, lons, lats, anomalies, colorrange = limits_ano, colormap = cgrad(:bluesreds, categorical = false), highclip = :red, lowclip = :blue)

    cl = @lift($m_v * " (W m⁻²)")
    cb = Colorbar(fig[1, 3], colorrange = limits_p, label = cl, height = 300, tellheight = false)

    cl_ano = @lift($m_v * " (W m⁻²)")
    cb_ano = Colorbar(fig[2, 3], colorrange = limits_ano, label = cl_ano, height = 300, tellheight = false, colormap = cgrad(:bluesreds, categorical = false), highclip = :red, lowclip = :blue)

    y_seasonal_means = @lift([cosine_weighted_global_mean(seasonal_y_data[$m_i][$m_v][season], lats) for season in ["DJF", "MAM", "JJA", "SON"]])
    y_seasonal_means_1 = @lift([cosine_weighted_global_mean(seasonal_y_data[1][$m_v][season], lats) for season in ["DJF", "MAM", "JJA", "SON"]])

    g_seasonal_means = @lift([cosine_weighted_global_mean(seasonal_g_data[$m_i][$m_m][$m_v][season], lats) for season in ["DJF", "MAM", "JJA", "SON"]])
    g_seasonal_means_1 = @lift([cosine_weighted_global_mean(seasonal_g_data[1][1][$m_v][season], lats) for season in ["DJF", "MAM", "JJA", "SON"]])

    min_sm = 0 # @lift(minimum(vcat($y_seasonal_means, $g_seasonal_means)))
    max_sm = @lift(maximum(vcat($y_seasonal_means, $g_seasonal_means)) + 10)
    limits_sm = @lift(($min_sm, $max_sm))

    line_y_1 = lines!(ax_sm, 1:4, y_seasonal_means_1, color= (:green, 0.3), linestyle = :dash)
    lines_g_1 = lines!(ax_sm, 1:4, g_seasonal_means_1, color= (:black, 0.3), linestyle = :dash)
    lines_g = lines!(ax_sm, 1:4, g_seasonal_means, color= :black)
    lines_y = lines!(ax_sm, 1:4, y_seasonal_means, color= :green)
    text!(ax_sm, 0.1, 0.1, text = rmse_y_g, align = (:left, :top), space = :relative)

    seasons = ["DJF", "MAM", "JJA", "SON"]
    current_s = @lift(findfirst(==($m_s), seasons))
    current_s_array = @lift([$current_s, $current_s])
    lines!(ax_sm, current_s_array, [0, 1000], color= :red, linewidth = 3)
    @lift(ylims!(ax_sm, $min_sm, $max_sm))
    axislegend(ax_sm, [lines_y, lines_g], ["era5", "ClimaLand"])

    fig
    return fig
end

function get_calibration_sets(directory="ekifiles")
    # List all files in the directory
    all_files = readdir(directory)

    # Find all EKI files
    eki_files = filter(f -> startswith(f, "eki_file_") && endswith(f, ".jld2"), all_files)

    # Create a mapping of base names to file pairs
    calibration_sets = Dict{String, Dict{String, String}}()

    for eki_file in eki_files
        # Extract the base name (remove prefix and extension)
        base_name = replace(eki_file, "eki_file_" => "", ".jld2" => "")

        # Check if corresponding prior file exists
        prior_file = "priors_$(base_name).jl"

        # Check if corresponding variable list file exists
        var_file = "variables_$(base_name).jl"

        # Only add to calibration sets if both required files exist
        if prior_file in all_files && var_file in all_files
            # Store the paths to all three files
            calibration_sets[base_name] = Dict(
                "eki" => joinpath(directory, eki_file),
                "prior" => joinpath(directory, prior_file),
                "variable_list" => joinpath(directory, var_file)
            )
        end
    end

    return calibration_sets
end

app = App(title="CliCal v0.2.0") do
    # Get available calibration sets
    calibration_sets = get_calibration_sets()
    calibration_names = sort(collect(keys(calibration_sets)))

    # Create a single menu for selecting calibration sets
    menu_calibration = Dropdown(calibration_names; index = 1) # only swu file. works!

    # Button to load data
    load_button = Button("Load Data")

    # Initialize figure
    fig = Figure(size = (1800, 1100), fontsize = 22)

    # Setup default axes
    ax_y = GM.GeoAxis(
                      fig[1, 1];
                      dest = "+proj=wintri",
                      title = "Era5 data (y)",
                     )
    lines!(ax_y, GM.coastlines())

    ax_g = GM.GeoAxis(
                      fig[1, 2];
                      dest = "+proj=wintri",
                      title = "ClimaLand (g)",
                     )
    lines!(ax_g, GM.coastlines())

    ax_anomalies = GM.GeoAxis(
                              fig[2, 2];
                              dest = "+proj=wintri",
                              title = "Anomalies: ClimaLand (g) - Era5 (y)",
                             )
    lines!(ax_anomalies, GM.coastlines())

    ax_sm = Axis(fig[2, 1],
                 title = "Seasonal means",
                 limits = (0.99, 4.01, 0, 400),
                 ylabel = "Value (W m⁻²)",
                 xticks = (1:4, ["DJF", "MAM", "JJA", "SON"]),
                 xlabel = "Season",
                )

    cal_set = menu_calibration.value

    # Get file paths for the selected calibration set using @lift
    selected_set = @lift(calibration_sets[$(cal_set)])
    eki_file = @lift($(selected_set)["eki"])
    prior_file = @lift($(selected_set)["prior"])
    variable_file = @lift($(selected_set)["variable_list"])

    # Load and process data
    loaded_data = Observable(load_and_process_data(eki_file[], prior_file[], variable_file[]))

    variable_list_vec = @lift($loaded_data["variable_list"])  # Default
    menu_var = Dropdown(variable_list_vec[])
    menu_iter = Dropdown(1:loaded_data[]["n_iterations"])
    menu_m = Dropdown(1:loaded_data[]["n_ensembles"])
    menu_season = Dropdown(["DJF", "MAM", "JJA", "SON"])

    year_x = @lift(2008+$(menu_iter.value))
    title_fig = @lift("$($(menu_season.value)) $($(menu_var.value)), iteration $($(menu_iter.value)), ensemble $($(menu_m.value)), year $($(year_x))")
    Label(fig[0, :], title_fig, fontsize=30, tellwidth = false)

    # Handle load button click with proper Observable handling
    on(load_button.value) do _
            # Get file paths for the selected calibration set without using @lift
            # to test - dev
            # menu_calibration.value[] = "swu_asnow_zenith"

            selected_set[] = calibration_sets[menu_calibration.value[]]
            eki_file[] = selected_set[]["eki"]
            prior_file[] = selected_set[]["prior"]
            variable_file[] = selected_set[]["variable_list"]

            # Load and process data
            loaded_data.val = load_and_process_data(eki_file[], prior_file[], variable_file[])

            # Update menus with new options
            variable_list_vec.val = loaded_data[]["variable_list"]
            menu_var.options.val = variable_list_vec[]
            menu_iter.options.val = 1:loaded_data[]["n_iterations"]
            menu_m.options.val = 1:loaded_data[]["n_ensembles"]

            # Important: Reset the values to the first item in each list to avoid KeyError
            menu_var.value[] = variable_list_vec[][1]    # Select first variable
            menu_iter.value[] = 1                  # Select first iteration
            menu_m.value[] = 1                     # Select first ensemble

            # now we can update safely
            variable_list_vec[] = loaded_data[]["variable_list"]
            menu_var.options[] = variable_list_vec[]
            menu_iter.options[] = 1:loaded_data[]["n_iterations"]
            menu_m.options[] = 1:loaded_data[]["n_ensembles"]
    end

    # Update display
    maps = update_fig(load_button, menu_calibration, menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, ax_anomalies, ax_sm, loaded_data[]["seasonal_g_data"], loaded_data[]["seasonal_y_data"], loaded_data[]["lons"], loaded_data[]["lats"])

    # Return the main layout
    return DOM.div(
                   style="display: grid; grid-template-columns: 300px 380px 1fr; gap: 20px; width: 100%;",

                   # Left column
                   DOM.div(
                           style="display: flex; flex-direction: column; gap: 20px;",
                           Card(
                                title="Calibration Set Selection",
                                DOM.div(
                                        style="display: flex; flex-direction: column; gap: 15px; padding: 10px;",
                                        DOM.div("Select Set", menu_calibration),
                                        DOM.div(load_button)
                                       )
                               ),
                           Card(
                                title="Parameters",
                                DOM.div(
                                        style="display: flex; flex-direction: column; gap: 15px; padding: 10px;",
                                        DOM.div("Variable", menu_var),
                                        DOM.div("Iteration", menu_iter),
                                        DOM.div("Ensemble", menu_m),
                                        DOM.div("Season", menu_season)
                                       )
                               ),
                           Card(
                                title="Error Metrics",
                                DOM.div(
                                        style="display: flex; flex-direction: column; gap: 15px; padding: 10px;",
                                        @lift(begin
                                                  if !haskey($loaded_data, "errors")
                                                      return DOM.div("Load data to see error metrics")
                                                  end

                                                  errors = $loaded_data["errors"]
                                                  normalized_errors = $loaded_data["normalized_errors"]
                                                  current_iter = $(menu_iter.value)

                                                  if current_iter == 0
                                                      return DOM.div("No iteration selected")
                                                  end

                                                  DOM.div(
                                                          style="display: flex; flex-direction: column; gap: 5px;",
                                                          DOM.div(
                                                                  style="display: grid; grid-template-columns: 0.8fr 1.1fr 1.1fr; border-bottom: 1px solid #999; padding-bottom: 5px; font-weight: bold;",
                                                                  DOM.span("Iteration"),
                                                                  DOM.span("Absolute"),
                                                                  DOM.span("Normalized (%)")
                                                                 ),
                                                          [DOM.div(
                                                                   style=string(
                                                                                "display: grid; grid-template-columns: 0.8fr 1.1fr 1.1fr; ",
                                                                                "padding: 5px; border-bottom: 1px solid #eee; ",
                                                                                i == current_iter ? "background-color: #f0f8ff; font-weight: bold;" : ""
                                                                               ),
                                                                   DOM.span("$(i)"),
                                                                   DOM.span(string(round(errors[i]/1e6, digits=2), " × 10⁶")),
                                                                   DOM.div(
                                                                           style=string(
                                                                                        "display: flex; align-items: center; ",
                                                                                        "color: ", i > 1 && normalized_errors[i] < normalized_errors[i-1] ? "green" : "inherit"
                                                                                       ),
                                                                           DOM.span(round(normalized_errors[i], digits=1)),
                                                                           i > 1 ? DOM.span(
                                                                                            style=string(
                                                                                                         "margin-left: 5px; font-size: 0.8em; ",
                                                                                                         "color: ", normalized_errors[i] < normalized_errors[i-1] ? "green" : "red"
                                                                                                        ),
                                                                                            normalized_errors[i] < normalized_errors[i-1] ? "↓" : "↑"
                                                                                           ) : DOM.span("")
                                                                          )
                                                                  ) for i in 1:length(errors)]
                                                         )
                                              end)
                                       )
                               )
                          ),

                   # Middle column
                   DOM.div(
                           style="display: flex; flex-direction: column; gap: 20px;",
                           Card(
                                title="Parameter Values",
                                DOM.div(
                                        style="display: flex; flex-direction: column; gap: 10px; padding: 10px;",
                                        @lift(begin
                                                  if !haskey($loaded_data, "param_dict")
                                                      return DOM.div("Load data to see parameter values")
                                                  end

                                                  param_dict = $loaded_data["param_dict"]
                                                  params_name = $loaded_data["params_name"]

                                                  current_iter = $(menu_iter.value)
                                                  current_m = $(menu_m.value)

                                                  if current_iter == 0 || current_m == 0
                                                      return DOM.div("No iteration or ensemble selected")
                                                  end

                                                  if !haskey(param_dict, current_iter) || length(param_dict[current_iter]) < current_m
                                                      return DOM.div("Parameter data not available for this iteration/ensemble")
                                                  end

                                                  param_values = param_dict[current_iter][current_m]
                                                  # For relative values, compare with initial parameters
                                                  params_initial = param_dict[1][1]
                                                  relative_values = param_values ./ params_initial

                                                  DOM.div(
                                                          style="display: flex; flex-direction: column; gap: 5px;",
                                                          DOM.div(
                                                                  style="display: grid; grid-template-columns: 1.2fr 1fr 1fr; border-bottom: 1px solid #999; padding-bottom: 5px; font-weight: bold;",
                                                                  DOM.span("Parameter"),
                                                                  DOM.span("Value"),
                                                                  DOM.span("Relative")
                                                                 ),
                                                          [DOM.div(
                                                                   style="display: grid; grid-template-columns: 1.2fr 1fr 1fr; padding: 5px; border-bottom: 1px solid #eee;",
                                                                   DOM.span(style="font-weight: bold;", params_name[i]),
                                                                   DOM.span(
                                                                            # Format based on magnitude
                                                                            let val = param_values[i]
                                                                                if abs(val) < 0.001 && val != 0
                                                                                    @sprintf("%.3e", val)  # Scientific for very small numbers
                                                                                elseif abs(val) > 10000
                                                                                    @sprintf("%.3e", val)  # Scientific for very large numbers
                                                                                else
                                                                                    @sprintf("%.4g", val)  # General format with 4 significant digits
                                                                                end
                                                                            end
                                                                           ),
                                                                   DOM.div(
                                                                           style=string(
                                                                                        "display: flex; align-items: center; ",
                                                                                        "color: ", relative_values[i] > 1 ? "green" : (relative_values[i] < 1 ? "red" : "black")
                                                                                       ),
                                                                           DOM.span(@sprintf("%.2f", relative_values[i])),  # Always 2 decimal places
                                                                           DOM.span(
                                                                                    style="margin-left: 5px; font-size: 0.8em;",
                                                                                    relative_values[i] > 1 ? "↑" : (relative_values[i] < 1 ? "↓" : "")
                                                                                   )
                                                                          )
                                                                  ) for i in 1:length(params_name)]
                                                         )
                                              end)
                                       )
                               ),
                           Card(
                                title="RMSE Metrics",
                                DOM.div(
                                        style="display: flex; flex-direction: column; gap: 15px; padding: 10px;",
                                        @lift(begin
                                                  if !haskey($loaded_data, "g_all")
                                                      return DOM.div("Load data to see RMSE metrics")
                                                  end

                                                  g_all = $loaded_data["g_all"]
                                                  y_all = $loaded_data["y_all"]
                                                  g_data = $loaded_data["g_data"]
                                                  y_data = $loaded_data["y_data"]
                                                  rmse_benchmarks = $loaded_data["rmse_benchmarks"]

                                                  current_iter = $(menu_iter.value)
                                                  current_var = $(menu_var.value)
                                                  current_m = $(menu_m.value)

                                                  if current_iter == 0 || current_var == "" || current_m == 0
                                                      return DOM.div("No iteration, variable, or ensemble selected")
                                                  end

                                                  rmse_clm_value = rmse_benchmarks[current_var]

                                                  DOM.div(
                                                          style="display: flex; flex-direction: column; gap: 15px;",
                                                          # Section headers
                                                          DOM.div(
                                                                  style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px;",
                                                                  DOM.h4(style="margin: 0; color: #333; grid-column: 1;", "Overall RMSE"),
                                                                  DOM.h4(style="margin: 0; color: #333; grid-column: 2;", "$(current_var) RMSE (LSMs: $(rmse_clm_value) W m⁻²)")
                                                                 ),
                                                          # Table headers
                                                          DOM.div(
                                                                  style="display: grid; grid-template-columns: 0.4fr 0.6fr 0.4fr 0.6fr; gap: 15px; border-bottom: 1px solid #999; padding-bottom: 5px; font-weight: bold;",
                                                                  DOM.span("Iter."),
                                                                  DOM.span("RMSE (W m⁻²)"),
                                                                  DOM.span("Iter."),
                                                                  DOM.span("RMSE (W m⁻²)")
                                                                 ),
                                                          # Table rows
                                                          [DOM.div(
                                                                   style="display: grid; grid-template-columns: 0.4fr 0.6fr 0.4fr 0.6fr; gap: 15px; padding: 5px; border-bottom: 1px solid #eee;",
                                                                   # Overall RMSE column
                                                                   DOM.span(
                                                                            style=i == current_iter ? "font-weight: bold;" : "",
                                                                            "$(i)"
                                                                           ),
                                                                   DOM.div(
                                                                           style=string(
                                                                                        "display: flex; align-items: center; ",
                                                                                        i == current_iter ? "font-weight: bold;" : "",
                                                                                        "color: ", i > 1 && RMSE(g_all[i][:,1], y_all[i]) < RMSE(g_all[i-1][:,1], y_all[i-1]) ? "green" : "inherit"
                                                                                       ),
                                                                           DOM.span(@sprintf("%.3f", RMSE(g_all[i][:,1], y_all[i]))),
                                                                           i > 1 ? DOM.span(
                                                                                            style=string(
                                                                                                         "margin-left: 5px; font-size: 0.8em; ",
                                                                                                         "color: ", RMSE(g_all[i][:,1], y_all[i]) < RMSE(g_all[i-1][:,1], y_all[i-1]) ? "green" : "red"
                                                                                                        ),
                                                                                            RMSE(g_all[i][:,1], y_all[i]) < RMSE(g_all[i-1][:,1], y_all[i-1]) ? "↓" : "↑"
                                                                                           ) : DOM.span("")
                                                                          ),
                                                                   # Variable RMSE column
                                                                   DOM.span(
                                                                            style=i == current_iter ? "font-weight: bold;" : "",
                                                                            "$(i)"
                                                                           ),
                                                                   DOM.div(
                                                                           style=string(
                                                                                        "display: flex; align-items: center; ",
                                                                                        i == current_iter ? "font-weight: bold;" : "",
                                                                                        "color: ", RMSE(g_data[i][current_m][current_var], y_data[i][current_var]) < rmse_clm_value ? "green" : "red"
                                                                                       ),
                                                                           DOM.span(@sprintf("%.3f", RMSE(g_data[i][current_m][current_var], y_data[i][current_var]))),
                                                                           DOM.span(
                                                                                    style="margin-left: 5px; font-size: 0.8em;",
                                                                                    i > 1 && RMSE(g_data[i][current_m][current_var], y_data[i][current_var]) <
                                                                                    RMSE(g_data[i-1][current_m][current_var], y_data[i-1][current_var]) ? "↓" :
                                                                                    (i > 1 ? "↑" : "")
                                                                                   )
                                                                          )
                                                                  ) for i in 1:length(g_all)]
                                                         )
                                              end)
                                       )
                               )
                          ),

                   # Right column (maps)
                   maps
                  )
end
# http://localhost:9384/browser-display

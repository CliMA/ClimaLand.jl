# Calibrated parameter file, error metrics, figures, and summary of the CRUJRA
# column-ensemble calibration. Included after `model_interface.jl`.

import CairoMakie
import GeoMakie
import Printf

const UNITS = Dict(
    "lhf" => "W m⁻²",
    "shf" => "W m⁻²",
    "swu" => "W m⁻²",
    "lwu" => "W m⁻²",
    "gpp" => "g C m⁻² day⁻¹",
)

"""
    write_parameter_file(path, prior, ekp)

Write the mean of the final parameter ensemble of `ekp`, in TOML, to `path`.
"""
function write_parameter_file(path, prior, ekp)
    names = EKP.get_name(prior)
    values = EKP.get_ϕ_mean_final(prior, ekp)
    open(path, "w") do io
        println(
            io,
            """
            # Calibrated by experiments/calibration/crujra_columns: $N_ITERATIONS iterations of
            # TransformUnscented EKI on $(length(read_columns().weight)) columns of the default land model (P-model,
            # prescribed MODIS LAI) forced by CRUJRA, against seasonal means of FLUXCOM LHF, SHF, GPP
            # and CERES SWU, LWU from $(Dates.Date(CALIBRATION_START)) to $(Dates.Date(STOP_DATE)).""",
        )
        haskey(ENV, "BUILDKITE_BUILD_URL") &&
            println(io, "# Buildkite build: $(ENV["BUILDKITE_BUILD_URL"])")
        haskey(ENV, "BUILDKITE_COMMIT") &&
            println(io, "# Commit: $(ENV["BUILDKITE_COMMIT"])")
        for (name, value) in zip(names, values)
            println(io, "\n[\"$name\"]")
            println(io, "value = $(round(value, sigdigits = 6))")
            println(io, "type = \"float\"")
        end
    end
    return nothing
end

"""
    error_metrics(simulated, observed, masks, weight)

Return, for each of `SHORT_NAMES`, the area-weighted RMSE and bias of the
`simulated` seasonal means over the calibrated entries (`masks`), and the loss:
the mean over all calibrated entries of the squared error normalized by the
observation noise variance (see `make_observation`).
"""
function error_metrics(simulated, observed, masks, weight)
    n_seasons = size(observed[first(SHORT_NAMES)], 1)
    relative_weight = entry_weights(weight, n_seasons)
    rmse, bias = Dict{String, Float64}(), Dict{String, Float64}()
    normalized_squared_errors = Float64[]
    for name in SHORT_NAMES
        mask = masks[name]
        w = relative_weight[mask]
        error = (vec(simulated[name]) .- vec(observed[name]))[mask]
        rmse[name] = sqrt(sum(w .* error .^ 2) / sum(w))
        bias[name] = sum(w .* error) / sum(w)
        append!(normalized_squared_errors, w .* error .^ 2 ./ NOISE_STD[name]^2)
    end
    loss = sum(normalized_squared_errors) / length(normalized_squared_errors)
    return (; rmse, bias, loss)
end

"""
    column_loss_change(default, calibrated, observed, masks)

Return, for each column, the change in its noise-normalized squared error
summed over variables and seasons, from `default` to `calibrated`.
"""
function column_loss_change(default, calibrated, observed, masks)
    n_seasons, n_columns = size(observed[first(SHORT_NAMES)])
    change = zeros(n_columns)
    for name in SHORT_NAMES
        valid = reshape(masks[name], n_seasons, n_columns)
        squared_error(simulated) =
            ifelse.(valid, (simulated .- observed[name]) .^ 2, 0.0)
        change .+=
            vec(
                sum(
                    squared_error(calibrated[name]) .-
                    squared_error(default[name]),
                    dims = 1,
                ),
            ) ./ NOISE_STD[name]^2
    end
    return change
end

"""
    plot_parameters(path, ekp, prior)

Plot the constrained parameters and the loss over the iterations of `ekp`.
"""
function plot_parameters(path, ekp, prior)
    n_parameters = sum(length.(EKP.batch(prior)))
    n_cols = 4
    n_rows = cld(n_parameters + 1, n_cols)
    fig = CairoMakie.Figure(size = (400 * n_cols, 320 * n_rows))
    for i in 1:n_parameters
        EKP.Visualize.plot_ϕ_over_iters(
            fig[fldmod1(i, n_cols)...],
            ekp,
            prior,
            i,
        )
    end
    EKP.Visualize.plot_error_over_iters(
        fig[fldmod1(n_parameters + 1, n_cols)...],
        ekp,
        error_metric = "loss",
    )
    CairoMakie.save(path, fig)
    return nothing
end

"""
    plot_rmse(path, metrics_default, metrics_calibrated)

Plot the area-weighted RMSE of each variable with the default and the
calibrated parameters.
"""
function plot_rmse(path, metrics_default, metrics_calibrated)
    fig = CairoMakie.Figure(size = (250 * length(SHORT_NAMES), 380))
    for (i, name) in enumerate(SHORT_NAMES)
        before, after =
            metrics_default.rmse[name], metrics_calibrated.rmse[name]
        ax = CairoMakie.Axis(
            fig[1, i],
            title = uppercase(name),
            ylabel = "RMSE [$(UNITS[name])]",
            xticks = (1:2, ["default", "calibrated"]),
        )
        CairoMakie.barplot!(
            ax,
            1:2,
            [before, after],
            color = [:gray, :firebrick],
        )
        CairoMakie.text!(
            ax,
            2,
            after;
            text = Printf.@sprintf("%+.0f%%", 100 * (after - before) / before),
            align = (:center, :bottom),
        )
        CairoMakie.ylims!(ax, 0, 1.15 * max(before, after))
    end
    CairoMakie.Label(
        fig[0, :],
        "Area-weighted RMSE of seasonal means over the calibration columns",
        font = :bold,
    )
    CairoMakie.save(path, fig)
    return nothing
end

"""
    plot_columns(path, longlat, weight, loss_change)

Map the columns, sized by area weight and colored by the change of their
noise-normalized squared error (negative where the calibration improved them).
"""
function plot_columns(path, longlat, weight, loss_change)
    fig = CairoMakie.Figure(size = (1200, 650))
    ax = GeoMakie.GeoAxis(
        fig[1, 1];
        dest = "+proj=robin",
        title = "Calibration columns",
    )
    CairoMakie.lines!(
        ax,
        GeoMakie.coastlines();
        color = :gray50,
        linewidth = 0.5,
    )
    limit = max(maximum(abs, loss_change), eps())
    plot = CairoMakie.scatter!(
        ax,
        first.(longlat),
        last.(longlat);
        color = loss_change,
        colormap = CairoMakie.Reverse(:RdBu),
        colorrange = (-limit, limit),
        markersize = 6 .+ 14 .* sqrt.(weight ./ maximum(weight)),
        strokecolor = :black,
        strokewidth = 0.4,
    )
    CairoMakie.Colorbar(
        fig[1, 2],
        plot;
        label = "Change in normalized squared error (calibrated − default)",
    )
    CairoMakie.save(path, fig)
    return nothing
end

function ClimaCalibrate.analyze_iteration(
    ::CRUJRAColumnsInterface,
    ekp,
    g_ensemble,
    prior,
    output_dir,
    iteration,
)
    @info "Iteration $iteration" loss = last(EKP.get_error(ekp)) parameters =
        Dict(zip(EKP.get_name(prior), EKP.get_ϕ_mean_final(prior, ekp)))
    results_dir = joinpath(output_dir, "results")
    mkpath(results_dir)
    plot_parameters(joinpath(results_dir, "parameters.png"), ekp, prior)
    return nothing
end

"""
    write_summary(path, prior, parameter_file, metrics_default, metrics_calibrated, artifact_dir)

Write a Markdown summary of the calibration to `path`, embedding the figures
uploaded as Buildkite artifacts under `artifact_dir`.
"""
function write_summary(
    path,
    prior,
    parameter_file,
    metrics_default,
    metrics_calibrated,
    artifact_dir,
)
    default_toml = LP.create_toml_dict(FT)
    calibrated_toml = TOML.parsefile(parameter_file)
    open(path, "w") do io
        println(
            io,
            """
            ### CRUJRA column calibration

            $N_ITERATIONS iterations of TransformUnscented EKI with $(2 * length(EKP.get_name(prior)) + 1) members on \
            $(length(read_columns().weight)) columns, forced by CRUJRA from $(Dates.Date(SPINUP_START)) to \
            $(Dates.Date(STOP_DATE)) (first year is spinup), against seasonal means of FLUXCOM LHF, SHF, GPP and \
            CERES SWU, LWU.

            | Parameter | Default | Calibrated |
            |---|---|---|""",
        )
        for name in EKP.get_name(prior)
            println(
                io,
                "| `$name` | $(round(default_toml[name], sigdigits = 4)) | $(round(calibrated_toml[name]["value"], sigdigits = 4)) |",
            )
        end
        println(
            io,
            """

            | Variable | RMSE default | RMSE calibrated | Bias default | Bias calibrated |
            |---|---|---|---|---|""",
        )
        for name in SHORT_NAMES
            values = (
                metrics_default.rmse[name],
                metrics_calibrated.rmse[name],
                metrics_default.bias[name],
                metrics_calibrated.bias[name],
            )
            println(
                io,
                "| $(uppercase(name)) [$(UNITS[name])] | ",
                join([Printf.@sprintf("%.3g", v) for v in values], " | "),
                " |",
            )
        end
        println(
            io,
            """

            Loss (noise-normalized squared error per observation): \
            $(Printf.@sprintf("%.3g", metrics_default.loss)) with the default parameters, \
            $(Printf.@sprintf("%.3g", metrics_calibrated.loss)) calibrated.

            <img src="artifact://$artifact_dir/parameters.png" alt="Parameters and loss over iterations">
            <img src="artifact://$artifact_dir/rmse.png" alt="RMSE before and after calibration">
            <img src="artifact://$artifact_dir/columns.png" alt="Calibration columns">

            <details><summary><code>toml/crujra_parameters.toml</code></summary>

            ```toml
            $(read(parameter_file, String))
            ```
            </details>""",
        )
    end
    return nothing
end

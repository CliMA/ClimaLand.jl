# Calibrated parameter file, error metrics, figures, and summary of the FLUXNET
# site calibration. Included after `model_interface.jl`.

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
    description(sites, ekp)

Return a description of the calibration that produced `ekp` on `sites`.
"""
description(sites, ekp) = """
    $(EKP.get_N_iterations(ekp)) iterations of TransformUnscented EKI with \
    $(EKP.get_N_ens(ekp)) members on $(length(sites.site_IDs)) FLUXNET sites, \
    forced by the site meteorology and compared to monthly means of the observed \
    LHF, SHF, SWU, LWU and GPP over up to $(length(CALIBRATION_YEARS)) years per site"""

"""
    write_parameter_file(path, sites, prior, ekp)

Write the model parameters set by the mean of the final parameter ensemble of
`ekp`, in TOML, to `path`.
"""
function write_parameter_file(path, sites, prior, ekp)
    calibrated =
        Dict(zip(EKP.get_name(prior), EKP.get_ϕ_mean_final(prior, ekp)))
    parameters = Dict(
        name => round(value, sigdigits = 6) for (name, value) in
        model_parameters(calibrated, LP.create_toml_dict(FT))
    )
    header = """
        # Calibrated by experiments/calibration/fluxnet_sites: $(description(sites, ekp)).
        """
    haskey(ENV, "BUILDKITE_BUILD_URL") &&
        (header *= "# Buildkite build: $(ENV["BUILDKITE_BUILD_URL"])\n")
    haskey(ENV, "BUILDKITE_COMMIT") &&
        (header *= "# Commit: $(ENV["BUILDKITE_COMMIT"])\n")
    write_parameter_toml(path, parameters; header)
    return nothing
end

"""
    error_metrics(simulated, observed, masks)

Return, for each of `SHORT_NAMES`, the RMSE and bias of the `simulated` monthly
means over the calibrated entries (`masks`), and the loss: the mean over all
calibrated entries of the squared error normalized by the observation noise
variance (see `make_observation`).
"""
function error_metrics(simulated, observed, masks)
    rmse, bias = Dict{String, Float64}(), Dict{String, Float64}()
    normalized_squared_errors = Float64[]
    for name in SHORT_NAMES
        error = (vec(simulated[name]) .- vec(observed[name]))[masks[name]]
        rmse[name] = isempty(error) ? NaN : sqrt(Statistics.mean(error .^ 2))
        bias[name] = isempty(error) ? NaN : Statistics.mean(error)
        append!(normalized_squared_errors, error .^ 2 ./ NOISE_STD[name]^2)
    end
    loss = Statistics.mean(normalized_squared_errors)
    return (; rmse, bias, loss)
end

"""
    site_loss_change(default, calibrated, observed, masks)

Return, for each site, the change in its noise-normalized squared error summed
over variables and months, from `default` to `calibrated`.
"""
function site_loss_change(default, calibrated, observed, masks)
    change = zeros(size(observed[first(SHORT_NAMES)], 2))
    for name in SHORT_NAMES
        valid = reshape(masks[name], size(observed[name]))
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

Plot the RMSE of each variable with the default and the calibrated parameters.
"""
function plot_rmse(path, metrics_default, metrics_calibrated)
    fig = CairoMakie.Figure(size = (250 * length(SHORT_NAMES), 380))
    for (i, name) in enumerate(SHORT_NAMES)
        before, after =
            metrics_default.rmse[name], metrics_calibrated.rmse[name]
        isfinite(before) || continue
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
        "RMSE of monthly means over the calibration sites",
        font = :bold,
    )
    CairoMakie.save(path, fig)
    return nothing
end

"""
    plot_sites(path, longlat, loss_change)

Map the sites, colored by the change of their noise-normalized squared error
(negative where the calibration improved them).
"""
function plot_sites(path, longlat, loss_change)
    fig = CairoMakie.Figure(size = (1200, 650))
    ax = GeoMakie.GeoAxis(
        fig[1, 1];
        dest = "+proj=robin",
        title = "Calibration sites",
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
        markersize = 10,
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

"""
    plot_seasonal_cycle(path, sites, default, calibrated, observed, masks)

Plot, for each variable and hemisphere, the mean by calendar month over the
calibrated entries of the observed, default, and calibrated monthly means.
"""
function plot_seasonal_cycle(path, sites, default, calibrated, observed, masks)
    calendar_month = Dates.month.(calibration_months())
    hemispheres = filter(
        h -> any(h.in_hemisphere, last.(sites.longlat)),
        [
            (; name = "Northern", in_hemisphere = >=(0)),
            (; name = "Southern", in_hemisphere = <(0)),
        ],
    )
    fig = CairoMakie.Figure(
        size = (300 * length(SHORT_NAMES), 300 * length(hemispheres)),
    )
    for (row, hemisphere) in enumerate(hemispheres),
        (col, name) in enumerate(SHORT_NAMES)

        valid = reshape(masks[name], size(observed[name]))
        valid[:, .!hemisphere.in_hemisphere.(last.(sites.longlat))] .= false
        ax = CairoMakie.Axis(
            fig[row, col],
            title = "$(uppercase(name)), $(hemisphere.name) hemisphere",
            xlabel = "Month",
            ylabel = UNITS[name],
            xticks = 1:12,
        )
        for (values, label, color) in (
            (observed[name], "observed", :black),
            (default[name], "default", :gray),
            (calibrated[name], "calibrated", :firebrick),
        )
            means = map(1:12) do month
                entries = valid[calendar_month .== month, :]
                any(entries) || return NaN
                Statistics.mean(values[calendar_month .== month, :][entries])
            end
            CairoMakie.lines!(ax, 1:12, means; label, color)
        end
        row == 1 && col == 1 && CairoMakie.axislegend(ax; position = :lt)
    end
    CairoMakie.save(path, fig)
    return nothing
end

function ClimaCalibrate.analyze_iteration(
    ::FluxnetSitesInterface,
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
    write_summary(path, sites, prior, ekp, parameter_file, metrics_default, metrics_calibrated, artifact_dir)

Write a Markdown summary of the calibration to `path`, embedding the figures
uploaded as Buildkite artifacts under `artifact_dir`.
"""
function write_summary(
    path,
    sites,
    prior,
    ekp,
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
            ### FLUXNET site calibration

            $(description(sites, ekp)).

            | Parameter | Default | Calibrated |
            |---|---|---|""",
        )
        for name in sort!(collect(keys(calibrated_toml)))
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
            <img src="artifact://$artifact_dir/seasonal_cycle.png" alt="Mean seasonal cycle over the sites">
            <img src="artifact://$artifact_dir/sites.png" alt="Calibration sites">

            <details><summary><code>$(basename(parameter_file))</code></summary>

            ```toml
            $(read(parameter_file, String))
            ```
            </details>""",
        )
    end
    return nothing
end

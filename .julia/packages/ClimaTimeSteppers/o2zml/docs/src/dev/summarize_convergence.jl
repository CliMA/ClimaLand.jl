import Plots
import JLD2
using LaTeXStrings: latexstring
using Printf: @sprintf

ENV["GKSwstype"] = "nul" # avoid displaying plots

include(joinpath(@__DIR__, "compute_convergence.jl"))

float_str(x) = @sprintf "%.4f" x
pow_str(x) = "10^{$(@sprintf "%.1f" log10(x))}"
function si_str(x)
    if isnan(x) || x in (0, Inf, -Inf)
        return string(x)
    end
    exponent = floor(Int, log10(x))
    mantissa = x / 10.0^exponent
    return "$(float_str(mantissa)) \\times 10^{$exponent}"
end

function summarize_convergence(test_name, test_results, default_dt, numerical_reference_info)
    rms_solution_str = "\\textrm{rms}\\_\\textrm{solution}"
    rms_error_str = "\\textrm{rms}\\_\\textrm{error}"
    average_rms_solution_str = "\\textrm{average}\\_$rms_solution_str"
    average_rms_error_str = "\\textrm{average}\\_$rms_error_str"
    default_dt_str = "\$dt = $(pow_str(default_dt))\$"
    reference_str =
        isnothing(numerical_reference_info) ? "Y_{analytic}(t)[\\textrm{index}]" : "Y_{ref}(t)[\\textrm{index}]"
    rms_solution_def_str = "\\textrm{RMS}_{\\textrm{index}}\\{Y(t)[\\textrm{index}]\\}"
    rms_error_def_str = "\\textrm{RMS}_{\\textrm{index}}\\{Y(t)[\\textrm{index}] - $reference_str\\}"
    average_rms_solution_def_str = "\\textrm{RMS}_t\\{$rms_solution_str(t)\\}"
    average_rms_error_def_str = "\\textrm{RMS}_t\\{$rms_error_str(t)\\}"
    footnote = """
        Terminology:
            \$$rms_solution_str(t) = $rms_solution_def_str\$
            \$$rms_error_str(t) = $rms_error_def_str\$
            \$$average_rms_solution_str = $average_rms_solution_def_str\$
            \$$average_rms_error_str = $average_rms_error_def_str\$"""
    if !isnothing(numerical_reference_info)
        (ref_alg_str, ref_dt, ref_average_rms_error) = numerical_reference_info
        ref_error_def_str = "Y_{ref}(t)[\\textrm{index}] - Y_{analytic}(t)[\\textrm{index}]"
        ref_average_rms_error_def_str = "\\textrm{RMS}_{t,\\,\\textrm{index}}\\{$ref_error_def_str\\}"
        footnote = "$footnote\n\n\nThe numerical reference solution \$Y_{ref}\$ is \
                    computed using $ref_alg_str with\n\$dt = $(pow_str(ref_dt)),\\ \
                    \\textrm{and}\\ $ref_average_rms_error_def_str = \
                    $(si_str(ref_average_rms_error))\$"
    end

    linestyles = (:solid, :dash, :dot, :dashdot, :dashdotdot)
    marker_kwargs = (; markershape = :circle, markeralpha = 0.5, markerstrokewidth = 0)
    plot_kwargs = (;
        legendposition = :outerright,
        legendtitlefontpointsize = 8,
        palette = :glasbey_bw_minc_20_maxl_70_n256,
        size = (1000, 2000), # size in px
        leftmargin = 60Plots.px,
        rightmargin = 0Plots.px,
        topmargin = 0Plots.px,
        bottommargin = 30Plots.px,
    )

    plot1_min = Inf
    plot1_max = -Inf
    plot1 = Plots.plot(;
        title = "Convergence Orders",
        xaxis = (latexstring("dt"), :log10),
        yaxis = (latexstring(average_rms_error_str), :log10),
        legendtitle = "Convergence Order (99% CI)",
        plot_kwargs...,
    )

    plot2a_min = Inf
    plot2a_max = -Inf
    plot2a = Plots.plot(;
        title = latexstring("Solutions with $default_dt_str"),
        xaxis = (latexstring("t"),),
        yaxis = (latexstring(rms_solution_str),),
        legendtitle = latexstring(average_rms_solution_str),
        plot_kwargs...,
    )

    plot2b_min = Inf
    plot2b_max = -Inf
    plot2b = Plots.plot(;
        title = latexstring("Errors with $default_dt_str"),
        xaxis = (latexstring("t"),),
        yaxis = (latexstring(rms_error_str), :log10),
        legendtitle = latexstring(average_rms_error_str),
        plot_kwargs...,
    )

    alg_strs = collect(keys(test_results))
    for alg_str in all_subtypes(ClimaTimeSteppers.AbstractAlgorithmName)
        if !(string(alg_str) in alg_strs)
            @warn "Convergence of $alg_str has not been tested with $test_name"
        end
    end
    sort_by(alg_str) = (
        test_results[alg_str]["predicted_order"],
        -round(test_results[alg_str]["average_rms_errors"][end]; sigdigits = 2),
        alg_str,
    ) # sort by predicted order, then by rounded error at large dt, then by name
    sort!(alg_strs; by = sort_by)
    for alg_str in alg_strs
        predicted_order = test_results[alg_str]["predicted_order"]
        predicted_super_convergence = test_results[alg_str]["predicted_super_convergence"]
        sampled_dts = test_results[alg_str]["sampled_dts"]
        average_rms_errors = test_results[alg_str]["average_rms_errors"]
        default_dt_times = test_results[alg_str]["default_dt_times"]
        default_dt_solutions = test_results[alg_str]["default_dt_solutions"]
        default_dt_errors = test_results[alg_str]["default_dt_errors"]

        linestyle = linestyles[(predicted_order - 1) % length(linestyles) + 1]
        order, order_uncertainty = convergence_order(sampled_dts, average_rms_errors, 0.99)
        order_str = "$(float_str(order)) \\pm $(float_str(order_uncertainty))"

        # Ignore large values from algorithms without reasonable order estimates.
        ignore_max = isnan(order_uncertainty) || order_uncertainty > 10

        plot1_min = min(plot1_min, minimum(filter(!isnan, average_rms_errors); init = Inf))
        ignore_max || (plot1_max = max(plot1_max, maximum(average_rms_errors)))
        plot1_label =
            predicted_super_convergence ? "$alg_str: \$$order_str\\ \\ \\ \\textbf{\\textit{SC}}\$" :
            "$alg_str: \$$order_str\$"
        Plots.plot!(
            plot1,
            sampled_dts,
            average_rms_errors;
            label = latexstring(plot1_label),
            linestyle,
            marker_kwargs...,
        )

        plot2a_min = min(plot2a_min, minimum(filter(!isnan, default_dt_solutions)))
        ignore_max || (plot2a_max = max(plot2a_max, maximum(default_dt_solutions)))
        Plots.plot!(
            plot2a,
            default_dt_times,
            default_dt_solutions;
            label = latexstring("$alg_str: \$$(si_str(rms(default_dt_solutions)))\$"),
            linestyle,
        )

        # Remove all 0s from default_dt_errors because they cannot be plotted on
        # a logarithmic scale.
        plot2b_min = min(plot2b_min, minimum(filter(!iszero, filter(!isnan, default_dt_errors))))
        ignore_max || (plot2b_max = max(plot2b_max, maximum(default_dt_errors)))
        default_dt_errors .= max.(default_dt_errors, eps(0.0))
        Plots.plot!(
            plot2b,
            default_dt_times,
            default_dt_errors;
            label = latexstring("$alg_str: \$$(si_str(rms(default_dt_errors)))\$"),
            linestyle,
        )
    end

    # Add a factor of 2 for padding in the log plots.
    Plots.plot!(plot1; ylim = (plot1_min / 2, plot1_max * 2))
    Plots.plot!(plot2b; ylim = (plot2b_min / 2, plot2b_max * 2))

    # Add 3% of the range for padding in the linear plot.
    plot2a_padding = (plot2a_max - plot2a_min) * 0.03
    Plots.plot!(plot2a; ylim = (plot2a_min - plot2a_padding, plot2a_max + plot2a_padding))

    footnote_plot = Plots.plot(;
        title = latexstring(footnote),
        titlelocation = :left,
        titlefontsize = 12,
        axis = nothing,
        framestyle = :none,
        margin = 0Plots.px,
        bottommargin = 10Plots.px,
    )

    all_plots = (plot1, plot2a, plot2b, footnote_plot)
    n_plots = length(all_plots)
    heights₀ = repeat([1 / (n_plots - 1)], n_plots - 1)
    footnote_height = eps()
    heights = [heights₀ .- footnote_height / n_plots..., footnote_height]
    layout = Plots.grid(n_plots, 1; heights)
    full_plot = Plots.plot(
        all_plots...;
        plot_title = "Convergence analysis for \"$test_name\"",
        layout = layout,
        fontfamily = "Computer Modern",
    )

    mkpath("output")
    file_suffix = lowercase(replace(test_name, " " => "_"))
    Plots.savefig(full_plot, joinpath("output", "convergence_$file_suffix.png"))
end

convergence_results_filenames = filter(readdir("output"; join = true)) do filename
    startswith(basename(filename), "convergence_") && endswith(filename, ".jld2")
end
all_convergence_results = map(JLD2.load_object, convergence_results_filenames)
representative_test_data = merge(all_convergence_results...)
test_names = collect(keys(representative_test_data))
for test_name in test_names, convergence_results in all_convergence_results
    test_name in keys(convergence_results) || continue
    @assert convergence_results[test_name]["default_dt"] == representative_test_data[test_name]["default_dt"]
    @assert convergence_results[test_name]["numerical_reference_info"] ==
            representative_test_data[test_name]["numerical_reference_info"]
end # Check that each test has a consistent default dt and reference solution
for test_name in test_names
    all_alg_results = map(all_convergence_results) do convergence_results
        test_name in keys(convergence_results) ? convergence_results[test_name]["all_alg_results"] : Dict()
    end
    summarize_convergence(
        test_name,
        merge(all_alg_results...),
        representative_test_data[test_name]["default_dt"],
        representative_test_data[test_name]["numerical_reference_info"],
    )
end

#= TO TEST CONVERGENCE OF ALL ALGORITHMS:
include("docs/src/dev/compute_convergence.jl")
for alg_name in all_subtypes(AbstractAlgorithmName)
    empty!(ARGS)
    push!(ARGS, "--alg", string(alg_name))
    @info "Testing convergence of $alg_name"
    include("docs/src/dev/report_gen_alg.jl")
end
include("docs/src/dev/summarize_convergence.jl")
=#

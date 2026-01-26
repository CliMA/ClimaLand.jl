import Plots
import JLD2
using ClimaCorePlots
using LinearAlgebra: norm
using PrettyTables: pretty_table, ft_printf

title_str(name) = titlecase(replace(string(name), '_' => ' '))

# Generates the Table 1 from
# "Optimization-based limiters for the spectral element method" by Guba et al.,
# and also plots the values used to generate the table.
function limiter_summary(sol_dicts, alg_strs)
    table_rows = []
    mkpath("output")
    for alg_str in alg_strs
        plots = []
        plot_kwargs = (;
            clims = (0, 1),
            color = :diverging_rainbow_bgymr_45_85_c67_n256,
            colorbar = false,
            guide = "",
            margin = 10Plots.px,
        )
        for use_limiter in (false, true), use_hyperdiffusion in (false, true)
            solution = sol_dicts[alg_str, use_hyperdiffusion, use_limiter]
            initial_q = solution[1].ρq ./ solution[1].ρ
            final_q = solution[end].ρq ./ solution[end].ρ
            names = propertynames(initial_q)

            if isempty(plots)
                for name in names
                    push!(plots, Plots.plot(initial_q.:($name); plot_kwargs..., title = title_str(name)))
                end
            end
            for name in names
                push!(plots, Plots.plot(final_q.:($name); plot_kwargs..., title = ""))
            end

            for name in names
                ϕ₀ = initial_q.:($name)
                ϕ = final_q.:($name)
                Δϕ₀ = maximum(ϕ₀) - minimum(ϕ₀)
                ϕ_error = ϕ .- ϕ₀
                table_row = [
                    alg_str;;
                    string(use_limiter);;
                    string(use_hyperdiffusion);;
                    title_str(name);;
                    max(0, -(minimum(ϕ) - minimum(ϕ₀)) / Δϕ₀);;
                    max(0, (maximum(ϕ) - maximum(ϕ₀)) / Δϕ₀);;
                    map(p -> norm(ϕ_error, p) / norm(ϕ₀, p), (1, 2, Inf))...
                ]
                push!(table_rows, table_row)
            end
        end
        colorbar_plot = Plots.scatter(
            [0];
            plot_kwargs...,
            colorbar = true,
            framestyle = :none,
            legend_position = :none,
            margin = 0Plots.px,
            markeralpha = 0,
            zcolor = [0],
        )
        plot = Plots.plot(
            plots...,
            colorbar_plot;
            layout = (Plots.@layout [Plots.grid(5, 3) a{0.1w}]),
            plot_title = "Tracer specific humidity for $alg_str (Initial, \
                          Final, Final w/ Hyperdiffusion, Final w/ Limiter, \
                          Final w/ Hyperdiffusion & Limiter)",
            size = (1600, 2000),
        )
        Plots.savefig(plot, joinpath("output", "limiter_summary_$alg_str.png"))
    end
    table = pretty_table(
        vcat(table_rows...);
        header = [
            "Algorithm",
            "Limiter",
            "Hyperdiffusion",
            "Tracer Name",
            "Max Undershoot",
            "Max Overshoot",
            "1-Norm Error",
            "2-Norm Error",
            "∞-Norm Error",
        ],
        crop = :none,
        body_hlines = collect(3:3:(length(table_rows) - 1)),
        formatters = ft_printf("%.4e"),
    )
    println(table)
end

limiter_results_filenames = filter(readdir("output"; join = true)) do filename
    startswith(basename(filename), "limiter_") && endswith(filename, ".jld2")
end
limiter_results = merge(map(JLD2.load_object, limiter_results_filenames)...)
alg_strs = unique(map(key -> key[1], collect(keys(limiter_results))))
limiter_summary(limiter_results, alg_strs)

#= TO ANALYZE LIMITERS WITH ARS343 AND SSP333:
for alg_str in ("ARS343", "SSP333"), use_limiter in (false, true), use_hyperdiffusion in (false, true)
    empty!(ARGS)
    push!(ARGS, "--alg", alg_str)
    push!(ARGS, "--use_limiter", string(use_limiter))
    push!(ARGS, "--use_hyperdiffusion", string(use_hyperdiffusion))
    @info "Running deformational flow with $alg_str, limiter = $use_limiter, hyperdif = $use_hyperdiffusion"
    include("docs/src/dev/limiter_analysis.jl")
end
include("docs/src/dev/limiter_summary.jl")
=#

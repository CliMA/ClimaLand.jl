"""
Run plot_eki_diagnostics.jl for all output folders matching a given
site_id and settingsdesc under the Neon_calibration output directory.

For each matched folder, finds the last iteration to set EKI_PATH,
then runs the diagnostic plot script.

Usage: julia --project=.buildkite experiments/calibrate_neon/plot_all_eki_diagnostics.jl
"""

using Dates

site_id      = "NEON-ster"
settingsdesc = "NeonExpTo005_10layers_newSWCInit_main20260504_dt180s"
base_dir     = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_calibration"
plot_script  = "/home/evametz/Github/ClimaLand/ClimaLand.jl/experiments/calibrate_neon/plot_eki_diagnostics.jl"
ENV["CALL_DEPTH"] = "0.02" #in X.XXm

N_ITERATIONS  = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
SPINUP_DAYS   = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
CALL_DEPTH    = get(ENV, "CALL_DEPTH", "0.02")  # in X.XXm
Caldepthnumstr = replace(string(CALL_DEPTH), "." => "_")

# Set to e.g. "output_2" to fix a specific output folder, or `nothing` to use the last one.
const output_n_override = nothing

function find_last_output(settings_dir)
    out_dirs = filter(readdir(settings_dir)) do d
        isdir(joinpath(settings_dir, d)) && (d == "output" || occursin(r"^output_\d+$", d))
    end
    isempty(out_dirs) && return nothing
    # numeric ordering: bare "output" → 0, "output_N" → N
    output_idx(d) = d == "output" ? 0 : parse(Int, d[8:end])
    return joinpath(settings_dir, out_dirs[argmax(output_idx.(out_dirs))])
end

function find_last_iteration(output_1_dir)
    iter_dirs = filter(readdir(output_1_dir)) do d
        isdir(joinpath(output_1_dir, d)) && startswith(d, "iteration_")
    end
    isempty(iter_dirs) && return nothing
    println("  Found iteration dirs: ", iter_dirs)
    return joinpath(output_1_dir, maximum(sort(iter_dirs)))
end


function collect_output_paths(base, sid, sdesc)
    site_dir = joinpath(base, sid)
    isdir(site_dir) || error("Site directory not found: $site_dir")

    paths = Tuple{String, String, String, String, String}[]  # (output_path, spinup, caldepth, n_it, date_range)

    for date_dir in readdir(site_dir; join=true)
        isdir(date_dir) || continue
        date_range = basename(date_dir)

        for spinup_dir in readdir(date_dir; join=true)
            isdir(spinup_dir) || continue
            spinup = basename(spinup_dir)
            spinup == "SpinUP-$(SPINUP_DAYS)d" || continue

            for depth_dir in readdir(spinup_dir; join=true)
                isdir(depth_dir) || continue
                caldepth = basename(depth_dir)
                caldepth == "CalDepth-$(Caldepthnumstr)M" || continue

                for nit_dir in readdir(depth_dir; join=true)
                    isdir(nit_dir) || continue
                    n_it = basename(nit_dir)
                    n_it == "$(N_ITERATIONS)-It" || continue

                    settings_dir = joinpath(nit_dir, sdesc)
                    isdir(settings_dir) || continue

                    has_output = any(readdir(settings_dir)) do d
                        isdir(joinpath(settings_dir, d)) &&
                            (d == "output" || occursin(r"^output_\d+$", d))
                    end
                    has_output || continue

                    push!(paths, (settings_dir, spinup, caldepth, n_it, date_range))
                end
            end
        end
    end
    return paths
end


function parse_final_parameters(path)
    params = Pair{String,String}[]
    for line in eachline(path)
        m = match(r"^\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+?)\s*$", line)
        m === nothing && continue
        push!(params, m.captures[1] => m.captures[2])
    end
    return params
end

function write_param_csv(rows, param_names, csv_path)
    open(csv_path, "w") do io
        println(io, join(["site_id", "time", param_names..., "parameterpath"], ","))
        for r in rows
            vals = [get(r.params, n, "") for n in param_names]
            println(io, join([r.site_id, r.date_range, vals..., r.parameterpath], ","))
        end
    end
end

function run_plots()
    entries = collect_output_paths(base_dir, site_id, settingsdesc)
    isempty(entries) && error("No output folders found for $site_id / $settingsdesc")

    println("Found $(length(entries)) output folder(s) for $site_id / $settingsdesc\n")

    param_rows = NamedTuple{(:site_id,:date_range,:params,:parameterpath),
                            Tuple{String,String,Dict{String,String},String}}[]

    for (output_path, spinup, caldepth, n_it, date_range) in entries
        println("=" ^ 70)
        println("Processing: $date_range | $spinup | $caldepth | $n_it")

        output_n = if output_n_override !== nothing
            joinpath(output_path, output_n_override)
        else
            find_last_output(output_path)
        end
        if output_n === nothing || !isdir(output_n)
            println("  WARNING: no output_N dir found, skipping.")
            continue
        end
        println("  Output dir: $(basename(output_n))")

        last_iter_dir = find_last_iteration(output_n)
        if last_iter_dir === nothing
            println("  WARNING: no iteration dirs found, skipping.")
            continue
        end

        eki_file = joinpath(last_iter_dir, "eki_file.jld2")
        if !isfile(eki_file)
            println("  WARNING: eki_file.jld2 not found at $eki_file, skipping.")
            continue
        end

        println("  Last iteration: $(basename(last_iter_dir))")
        println("  EKI file: $eki_file")

        param_file = joinpath(output_n, "final_parameter_means.txt")
        if isfile(param_file)
            pairs = parse_final_parameters(param_file)
            push!(param_rows, (site_id = site_id, date_range = date_range,
                               params = Dict(pairs), parameterpath = param_file))
        else
            println("  NOTE: $param_file not found, skipping param row.")
        end

        ENV["NEON_SITE_ID"]     = site_id
        ENV["CALL_OUTPUT_PATH"] = output_path
        ENV["CALL_OUTPUT_1"]    = output_n
        ENV["CALL_DEPTH"]       = CALL_DEPTH
        ENV["CALL_EKI_PATH"]    = eki_file

        t0 = time()
        try
            include(plot_script)
            println("  Done in $(round(time()-t0; digits=1)) s")
        catch err
            println("  FAILED after $(round(time()-t0; digits=1)) s")
            showerror(stdout, err, catch_backtrace())
            println()
        end
    end

    if !isempty(param_rows)
        param_names = String[]
        seen = Set{String}()
        for r in param_rows, n in keys(r.params)
            n in seen || (push!(param_names, n); push!(seen, n))
        end
        comparison_dir = joinpath(base_dir, site_id, "Parameter_comparison")
        mkpath(comparison_dir)
        stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        csv_path = joinpath(comparison_dir, "$(settingsdesc)_$(stamp).csv")
        write_param_csv(param_rows, param_names, csv_path)
        println("\nWrote parameter CSV → $csv_path")
    else
        println("\nNo parameter files collected, no CSV written.")
    end

    println("\nAll done.")
end

run_plots()

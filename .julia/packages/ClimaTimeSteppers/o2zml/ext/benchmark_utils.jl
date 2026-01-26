#####
##### BenchmarkTools's trial utils
#####

get_summary(trial, trial_step = nothing) = (;
    # Using some BenchmarkTools internals :/
    mem = BenchmarkTools.prettymemory(trial.memory),
    mem_val = trial.memory,
    nalloc = trial.allocs,
    t_min = BenchmarkTools.prettytime(minimum(trial.times)),
    t_max = BenchmarkTools.prettytime(maximum(trial.times)),
    t_mean = BenchmarkTools.prettytime(StatsBase.mean(trial.times)),
    t_mean_val = StatsBase.mean(trial.times),
    t_med = BenchmarkTools.prettytime(StatsBase.median(trial.times)),
    n_samples = length(trial),
    percentage = isnothing(trial_step) ? -1 : minimum(trial.times) / minimum(trial_step.times) * 100,
)

function tabulate_summary(summary; n_calls_per_step)
    summary_keys = collect(keys(summary))
    mem = map(k -> summary[k].mem, summary_keys)
    nalloc = map(k -> summary[k].nalloc, summary_keys)
    t_mean = map(k -> summary[k].t_mean, summary_keys)
    t_min = map(k -> summary[k].t_min, summary_keys)
    t_max = map(k -> summary[k].t_max, summary_keys)
    t_med = map(k -> summary[k].t_med, summary_keys)
    n_samples = map(k -> summary[k].n_samples, summary_keys)
    percentage = map(k -> summary[k].percentage, summary_keys)

    func_names = if isnothing(n_calls_per_step)
        map(k -> string(k), collect(keys(summary)))
    else
        @info "(#)x entries have been multiplied by corresponding factors in order to compute percentages"
        map(k -> string(k, " ($(n_calls_per_step[k])x)"), collect(keys(summary)))
    end
    table_data = hcat(func_names, mem, nalloc, t_min, t_max, t_mean, t_med, n_samples, percentage)

    header = (
        ["Function", "Memory", "allocs", "Time", "Time", "Time", "Time", "N-samples", "step! percentage"],
        [" ", "estimate", "estimate", "min", "max", "mean", "median", "", ""],
    )

    PrettyTables.pretty_table(
        table_data;
        header,
        crop = :none,
        alignment = vcat(:l, repeat([:r], length(header[1]) - 1)),
    )
end

get_trial(f::Nothing, args, name; device, with_cu_prof = :bprofile, trace = false, crop = false, hcrop = nothing) =
    nothing
function get_trial(f, args, name; device, with_cu_prof = :bprofile, trace = false, crop = false, hcrop = nothing)
    f(args...) # compile first
    b = if device isa ClimaComms.CUDADevice
        BenchmarkTools.@benchmarkable CUDA.@sync $f($(args)...)
    else
        BenchmarkTools.@benchmarkable $f($(args)...)
    end
    sample_limit = 10
    println("--------------- Benchmarking/profiling $name...")
    trial = BenchmarkTools.run(b, samples = sample_limit)
    if device isa ClimaComms.CUDADevice
        p = if with_cu_prof == :bprofile
            CUDA.@bprofile trace = trace f(args...)
        else
            CUDA.@profile trace = trace f(args...)
        end
        if crop
            println(p) # crops by default
        else
            # use "COLUMNS" to set how many horizontal characters to crop:
            # See https://github.com/ronisbr/PrettyTables.jl/issues/11#issuecomment-2145550354
            envs = isnothing(hcrop) ? () : ("COLUMNS" => hcrop,)
            withenv(envs...) do
                io = IOContext(stdout, :crop => :horizontal, :limit => true, :displaysize => displaysize())
                show(io, p)
            end
            println()
        end
    end
    println()
    return trial
end

get_W(i::CTS.DistributedODEIntegrator) = hasproperty(i.cache, :W) ? i.cache.W : i.cache.newtons_method_cache.j
get_W(i) = i.cache.W
f_args(i, f::CTS.ForwardEulerODEFunction) = (copy(i.u), i.u, i.p, i.t, i.dt)
f_args(i, f) = (similar(i.u), i.u, i.p, i.t)

r_args(i, f::CTS.ForwardEulerODEFunction) = (copy(i.u), copy(i.u), i.u, i.p, i.t, i.dt)
r_args(i, f) = (similar(i.u), similar(i.u), i.u, i.p, i.t)

implicit_args(i::CTS.DistributedODEIntegrator) = f_args(i, i.sol.prob.f.T_imp!)
implicit_args(i) = f_args(i, i.f.f1)
remaining_args(i::CTS.DistributedODEIntegrator) = r_args(i, i.sol.prob.f.T_exp_T_lim!)
remaining_args(i) = r_args(i, i.f.f2)
wfact_fun(i) = implicit_fun(i).Wfact
implicit_fun(i::CTS.DistributedODEIntegrator) = i.sol.prob.f.T_imp!
implicit_fun(i) = i.sol.prob.f.f1
remaining_fun(i::CTS.DistributedODEIntegrator) = i.sol.prob.f.T_exp_T_lim!
remaining_fun(i) = i.sol.prob.f.f2

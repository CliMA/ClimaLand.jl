# To run
#=
using HashArrayMappedTries, PkgBenchmark
result = benchmarkpkg(HashArrayMappedTries)
export_markdown("perf.md", result)
=#
using BenchmarkTools
using HashArrayMappedTries

function create_dict(::Type{Dict}, n) where Dict
    dict=Dict{Int, Int}()
    for i in 1:n
        dict[i] = i
    end
    dict
end

function HashArrayMappedTries.insert(dict::Base.Dict{K, V}, key::K, v::V) where {K,V}
    dict = copy(dict)
    dict[key] = v
    return dict
end

function HashArrayMappedTries.delete(dict::Base.Dict{K}, key::K) where K
    dict = copy(dict)
    delete!(dict, key)
    return dict
end

function create_persistent_dict(::Type{Dict}, n) where Dict
    dict = Dict{Int, Int}()
    for i in 1:n
        dict = insert(dict, i, i)
    end
    dict
end

# PkgBenchmark ignores evals=1 so this leads to invalid results. 
# https://github.com/JuliaCI/BenchmarkTools.jl/issues/328

function create_benchmark(::Type{Dict}) where Dict
    group = BenchmarkGroup()
    group["creation, size=0"] = @benchmarkable create_dict($Dict, 0)
    group["creation (Persistent), size=0"] = @benchmarkable create_persistent_dict($Dict, 0)
    # group["setindex!, size=0"] = @benchmarkable dict[1] = 1 setup=(dict=create_dict($Dict, 0)) evals=1
    group["insert, size=0"] = @benchmarkable insert(dict, 1, 1) setup=(dict=create_persistent_dict($Dict, 0))
    for i in 0:14
        N = 2^i
        group["creation, size=$N"]              = @benchmarkable create_dict($Dict, $N)
        group["creation (Persistent), size=$N"] = @benchmarkable create_persistent_dict($Dict, $N)
        # group["setindex!, size=$N"]             = @benchmarkable dict[$N+1] = $N+1        setup=(dict=create_dict($Dict, $N)) evals=1
        group["getindex, size=$N"]              = @benchmarkable dict[$N]                 setup=(dict=create_dict($Dict, $N))
        # group["delete!, size=$N"]               = @benchmarkable delete!(dict, $N)        setup=(dict=create_dict($Dict, $N)) evals=1

        # Persistent
        group["insert, size=$N"]                = @benchmarkable insert(dict, $N+1, $N+1) setup=(dict=create_persistent_dict($Dict, $N))
        group["delete, size=$N"]                = @benchmarkable delete(dict, $N)         setup=(dict=create_persistent_dict($Dict, $N))
    end
    return group
end

const SUITE = BenchmarkGroup()
SUITE["Base.Dict"] = create_benchmark(Dict)
SUITE["HAMT"] = create_benchmark(HAMT)

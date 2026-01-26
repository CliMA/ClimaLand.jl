using Dates
using BenchmarkTools
using Statistics
using UUIDs
using Pkg

function compute(n)
    t0 = DateTime(1900, 1, 1) .+ Dates.Second.(0:(n - 1))
    t1 = DateTime(2000, 1, 1) .+ Dates.Second.(0:(n - 1))

    diff = t1 - reverse(t0)

    return mean(Dates.value.(Dates.Millisecond.(diff)) ./ 1000)
end

println("julia: ", VERSION)

pkg_name = "Dates"
m = Pkg.Operations.Context().env.manifest
println("Dates: ", m[findfirst(v -> v.name == pkg_name, m)].version)

n = 1_000_000
#n = 100_000
println("mean_total_seconds: ", compute(n))

bm = run(@benchmarkable compute(n) samples = 100)

println("min time: ", minimum(bm.times / 1.0e9))

open("julia-Dates.txt", "w") do f
    for t in bm.times
        println(f, t / 1.0e9)
    end
end

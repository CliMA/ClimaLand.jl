using SymbolicIndexingInterface
using SymbolicIndexingInterface: parameter_timeseries
using Test

struct MyDiffEqArray
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end
SymbolicIndexingInterface.current_time(mda::MyDiffEqArray) = mda.t
SymbolicIndexingInterface.state_values(mda::MyDiffEqArray) = mda.u
SymbolicIndexingInterface.is_timeseries(::Type{MyDiffEqArray}) = Timeseries()

ps = ones(3)
@test_throws ArgumentError ParameterTimeseriesCollection((ones(3), 2ones(3)), ps)

a_timeseries = MyDiffEqArray(collect(0:0.1:0.9), [[2.5i, sin(0.2i)] for i in 1:10])
b_timeseries = MyDiffEqArray(collect(0:0.25:0.9), [[3.5i, log(1.3i)] for i in 1:4])
c_timeseries = MyDiffEqArray(collect(0:0.17:0.90), [[4.3i] for i in 1:5])
collection = (a_timeseries, b_timeseries, c_timeseries)
ptc = ParameterTimeseriesCollection(collection, ps)

@test collect(eachindex(ptc)) == [1, 2, 3]
@test [x for x in ptc] == [a_timeseries, b_timeseries, c_timeseries]
@test length(ptc) == 3
@test parent(ptc) === collection
@test parameter_values(ptc) === ps

for i in 1:3
    @test ptc[i] === collection[i]
    @test parameter_timeseries(ptc, i) == collection[i].t
    for j in eachindex(collection[i].u[1])
        pti = ParameterTimeseriesIndex(i, j)
        @test ptc[pti] == getindex.(collection[i].u, j)
        for k in eachindex(collection[i].u)
            rhs = collection[i].u[k][j]
            @test ptc[pti, CartesianIndex(k)] == rhs
            @test ptc[pti, k] == rhs
            @test ptc[i, k] == collection[i].u[k]
            @test ptc[i, k, j] == rhs
            @test parameter_values(ptc, pti, k) == rhs
        end
        allidxs = eachindex(collection[i].u)
        for subidx in [:, rand(allidxs, 3), rand(Bool, length(allidxs))]
            rhs = getindex.(collection[i].u[subidx], j)
            @test ptc[pti, subidx] == rhs
            @test ptc[i, subidx, j] == rhs
            @test parameter_values(ptc, pti, subidx) == rhs
        end
    end
end

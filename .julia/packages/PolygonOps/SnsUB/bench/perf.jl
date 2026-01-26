using BenchmarkTools
using PolygonOps
using StaticArrays
using Test

const suite = BenchmarkGroup()
suite["inpolygon"] = BenchmarkGroup()

circle(t) = SVector(cos(t), sin(t))

circle_poly = map(circle, 0:(2pi/1000000):2pi)

#@show length(circle_poly)

# using Plots
# display(plot(map(x -> x[1], circle_poly),map(x -> x[2], circle_poly)))

# end == begin
push!(circle_poly, circle_poly[1])
for algo in (HaoSun(),HormannAgathos())
    @test inpolygon(SVector(0.1,0),circle_poly,algo) == 1
    @test inpolygon(SVector(0.,2.),circle_poly,algo) == 0
    @test inpolygon(SVector(2.,0.),circle_poly,algo) == 0
    @test inpolygon(SVector(1.,0.),circle_poly,algo) == -1
    @test inpolygon(SVector(0.,-2.),circle_poly,algo) == 0
end

suite["inpolygon"]["circle, HaoSun"] = @benchmarkable inpolygon(SVector(0,0),circle_poly, HaoSun())
suite["inpolygon"]["circle, HormanAgathos"] = @benchmarkable inpolygon(SVector(0,0),circle_poly, HormannAgathos())


results = run(suite)

for trial in results
    ctx = IOContext(stdout, :verbose => true, :compact => false)
    println(ctx)
    println(ctx, trial.first)
    println(ctx, trial.second)
end

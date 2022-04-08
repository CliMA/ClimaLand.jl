using ClimaCore
using Test
using StaticArrays
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
import ClimaLSM:
    prognostic_vars, prognostic_types, auxiliary_vars, auxiliary_types, name
using ClimaLSM.Domains: HybridBox, Column, Point
using ClimaLSM.Domains: coordinates

struct Model{FT, D} <: AbstractModel{FT}
    domain::D
end
ClimaLSM.name(m::Model) = :foo
ClimaLSM.prognostic_vars(m::Model) = (:a, :b)
ClimaLSM.auxiliary_vars(m::Model) = (:d, :e)

ClimaLSM.prognostic_types(m::Model{FT}) where {FT} = (FT, SVector{2, FT})

ClimaLSM.auxiliary_types(m::Model{FT}) where {FT} = (FT, SVector{2, FT})




FT = Float64
zmin = FT(1.0)
zmax = FT(2.0)
zlim = (zmin, zmax)
nelements = 5



@testset "Column Domain Vars" begin
    column = Column(; zlim = zlim, nelements = nelements)
    m = Model{FT, typeof(column)}(column)
    Y, p, coords = initialize(m)
    @test parent(Y.foo.a) == zeros(FT, 5, 1)
    @test parent(Y.foo.b) == zeros(FT, 5, 2)
    @test parent(p.foo.d) == zeros(FT, 5, 1)
    @test parent(p.foo.e) == zeros(FT, 5, 2)
end

@testset "Point Domain Vars" begin
    point = Point(; z_sfc = zlim[1])
    m = Model{FT, typeof(point)}(point)
    Y, p, coords = initialize(m)
    @test Y.foo.a == [0.0]
    @test Y.foo.b == [zero(SVector{2, FT})]
    @test p.foo.d == [0.0]
    @test p.foo.e == [zero(SVector{2, FT})]
end

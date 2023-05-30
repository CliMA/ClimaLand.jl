using ClimaCore
using Test
using StaticArrays
using ClimaLSM
import ClimaLSM:
    AbstractModel,
    AbstractImExModel,
    AbstractExpModel,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    auxiliary_types,
    name
using ClimaLSM.Domains: HybridBox, Column, Point
using ClimaLSM.Domains: coordinates
include("../src/Vegetation/component_models.jl")

@testset "Default model" begin
    struct DefaultModel{FT} <: AbstractModel{FT} end
    dm = DefaultModel{Float32}()
    @test ClimaLSM.prognostic_vars(dm) == ()
    @test ClimaLSM.prognostic_types(dm) == ()
    @test ClimaLSM.auxiliary_vars(dm) == ()
    @test ClimaLSM.auxiliary_types(dm) == ()

    x = [0, 1, 2, 3]
    dm_exp_tendency! = make_exp_tendency(dm)
    @test dm_exp_tendency!(x[1], x[2], x[3], x[4]) == nothing
    @test x == [0, 1, 2, 3]

    dm_compute_exp_tendency! = make_compute_exp_tendency(dm)
    @test dm_compute_exp_tendency!(x[1], x[2], x[3], x[4]) == nothing
    @test x == [0, 1, 2, 3]

    dm_imp_tendency! = make_imp_tendency(dm)
    @test dm_imp_tendency!(x[1], x[2], x[3], x[4]) == nothing
    @test x == [0, 1, 2, 3]

    dm_compute_imp_tendency! = make_compute_imp_tendency(dm)
    @test dm_compute_imp_tendency!(x[1], x[2], x[3], x[4]) == nothing
    @test x == [0, 1, 2, 3]

    dm_update_aux! = make_update_aux(dm)
    @test dm_update_aux!(x[1], x[2], x[3]) == nothing
    @test x == [0, 1, 2, 3]
end

@testset "Default ImEx model" begin
    struct DefaultImExModel{FT} <: AbstractImExModel{FT} end
    dm_imex = DefaultImExModel{Float32}()

    x = [0, 1, 2, 3]
    dm_imp_tendency! = make_imp_tendency(dm_imex)
    @test dm_imp_tendency!(x[1], x[2], x[3], x[4]) == nothing
    @test x == [0, 1, 2, 3]

    dm_compute_imp_tendency! = make_compute_imp_tendency(dm_imex)
    @test dm_compute_imp_tendency!(x[1], x[2], x[3], x[4]) == nothing
    @test x == [0, 1, 2, 3]

end

@testset "Default canopy component" begin
    struct DefaultCanopyComponent{FT} <: AbstractCanopyComponent{FT} end
    dcc = DefaultCanopyComponent{Float32}()
    @test ClimaLSM.prognostic_vars(dcc) == ()
    @test ClimaLSM.prognostic_types(dcc) == ()
    @test ClimaLSM.auxiliary_vars(dcc) == ()
    @test ClimaLSM.auxiliary_types(dcc) == ()

    x = [0, 1, 2, 3]
    dcc_compute_exp_tendency! = make_compute_exp_tendency(dcc, nothing)
    @test dcc_compute_exp_tendency!(x[1], x[2], x[3], x[4]) == nothing
    @test x == [0, 1, 2, 3]

    @test_throws MethodError make_compute_imp_tendency(dcc)
end

struct Model{FT, D} <: AbstractModel{FT}
    domain::D
end
ClimaLSM.name(m::Model) = :foo
ClimaLSM.prognostic_vars(m::Model) = (:a, :b)
ClimaLSM.auxiliary_vars(m::Model) = (:d, :e)

ClimaLSM.prognostic_types(m::Model{FT}) where {FT} = (FT, SVector{2, FT})

ClimaLSM.auxiliary_types(m::Model{FT}) where {FT} = (FT, SVector{2, FT})

@testset "Column Domain Vars" begin
    FT = Float64
    zmin = FT(1.0)
    zmax = FT(2.0)
    zlim = (zmin, zmax)
    nelements = 5

    column = Column(; zlim = zlim, nelements = nelements)
    m = Model{FT, typeof(column)}(column)
    Y, p, coords = initialize(m)
    @test parent(Y.foo.a) == zeros(FT, 5, 1)
    @test parent(Y.foo.b) == zeros(FT, 5, 2)
    @test parent(p.foo.d) == zeros(FT, 5, 1)
    @test parent(p.foo.e) == zeros(FT, 5, 2)
end

@testset "Point Domain Vars" begin
    FT = Float64
    zmin = FT(1.0)
    zmax = FT(2.0)
    zlim = (zmin, zmax)
    nelements = 5

    point = Point(; z_sfc = zlim[1])
    m = Model{FT, typeof(point)}(point)
    Y, p, coords = initialize(m)
    @test Y.foo.a[] == 0.0
    @test Y.foo.b[] == zero(SVector{2, FT})
    @test p.foo.d[] == 0.0
    @test p.foo.e[] == zero(SVector{2, FT})
end

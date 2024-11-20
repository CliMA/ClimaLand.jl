import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Test
using StaticArrays
using ClimaLand
import ClimaLand:
    AbstractModel,
    AbstractImExModel,
    AbstractExpModel,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    auxiliary_types,
    prognostic_domain_names,
    auxiliary_domain_names,
    name
using ClimaLand.Domains: HybridBox, Column, Point
using ClimaLand.Domains: coordinates
using ClimaLand.Canopy: AbstractCanopyComponent

FT = Float32
@testset "Default model, FT = $FT" begin
    struct DefaultModel{FT} <: AbstractModel{FT} end
    ClimaLand.name(::DefaultModel) = :default
    dm = DefaultModel{FT}()
    @test ClimaLand.prognostic_vars(dm) == ()
    @test ClimaLand.prognostic_types(dm) == ()
    @test ClimaLand.prognostic_domain_names(dm) == ()
    @test ClimaLand.auxiliary_vars(dm) == ()
    @test ClimaLand.auxiliary_types(dm) == ()
    @test ClimaLand.auxiliary_domain_names(dm) == ()

    tendency_args = ((default = [1],), 2, 3, 4)
    dm_exp_tendency! = make_exp_tendency(dm)
    @test dm_exp_tendency!(tendency_args...) == [0]
    @test tendency_args == ((default = [0],), 2, 3, 4)

    tendency_args = ((default = [1],), 2, 3, 4)
    dm_compute_exp_tendency! = make_compute_exp_tendency(dm)
    @test dm_compute_exp_tendency!(tendency_args...) == [0]
    @test tendency_args == ((default = [0],), 2, 3, 4)

    tendency_args = ((default = [1],), 2, 3, 4)
    dm_imp_tendency! = make_imp_tendency(dm)
    @test dm_imp_tendency!(tendency_args...) == [0]
    @test tendency_args == ((default = [0],), 2, 3, 4)

    tendency_args = ((default = [1],), 2, 3, 4)
    dm_compute_imp_tendency! = make_compute_imp_tendency(dm)
    @test dm_compute_imp_tendency!(tendency_args...) == [0]
    @test tendency_args == ((default = [0],), 2, 3, 4)

    x = [1, 2, 3]
    dm_update_aux! = make_update_aux(dm)
    @test isnothing(dm_update_aux!(x...))
    @test x == [1, 2, 3]

    @test ClimaLand.get_drivers(dm) == ()
    @test ClimaLand.add_drivers_to_cache((;), dm, nothing) == (;)
end

@testset "Default ImEx model, FT = $FT" begin
    struct DefaultImExModel{FT} <: AbstractImExModel{FT} end
    dm_imex = DefaultImExModel{FT}()
    ClimaLand.name(::DefaultImExModel) = :default_imex

    tendency_args = ((default_imex = [1],), 2, 3, 4)
    dm_imp_tendency! = make_imp_tendency(dm_imex)
    @test dm_imp_tendency!(tendency_args...) == [0]
    @test tendency_args == ((default_imex = [0],), 2, 3, 4)

    tendency_args = ((default_imex = [1],), 2, 3, 4)
    dm_compute_imp_tendency! = make_compute_imp_tendency(dm_imex)
    @test dm_compute_imp_tendency!(tendency_args...) == [0]
    @test tendency_args == ((default_imex = [0],), 2, 3, 4)
end

@testset "Default canopy component" begin
    FT = Float64
    struct DefaultCanopy{FT} <: AbstractImExModel{FT} end
    ClimaLand.name(::DefaultCanopy) = :canopy
    dc = DefaultCanopy{FT}()

    struct DefaultCanopyComponent{FT} <: AbstractCanopyComponent{FT} end
    ClimaLand.name(::DefaultCanopyComponent) = :component
    dcc = DefaultCanopyComponent{FT}()

    @test ClimaLand.prognostic_vars(dcc) == ()
    @test ClimaLand.prognostic_types(dcc) == ()
    @test ClimaLand.prognostic_domain_names(dcc) == ()

    @test ClimaLand.auxiliary_vars(dcc) == ()
    @test ClimaLand.auxiliary_types(dcc) == ()
    @test ClimaLand.auxiliary_domain_names(dcc) == ()

    x = [(canopy = (component = [1],),), 1, 2, 3]
    dcc_compute_exp_tendency! = make_compute_exp_tendency(dcc, dc)
    @test dcc_compute_exp_tendency!(x[1], x[2], x[3], x[4]) == nothing

    @test_throws MethodError make_compute_imp_tendency(dcc)
end

@testset "Model with no prognostic vars, FT = $FT" begin
    struct Model{FT, D} <: AbstractModel{FT}
        domain::D
    end

    zmin = FT(1.0)
    zmax = FT(2.0)
    zlim = (zmin, zmax)
    nelements = 5

    column = Column(; zlim = zlim, nelements = nelements)
    m = Model{FT, typeof(column)}(column)

    # `initialize` should fail for a model with no prognostic variables
    @test_throws AssertionError initialize(m)
end

@testset "Variables, FT = $FT" begin
    struct Foo{FT, D} <: AbstractModel{FT}
        domain::D
    end
    ClimaLand.name(m::Foo) = :foo
    ClimaLand.prognostic_vars(m::Foo) = (:a, :b)
    ClimaLand.auxiliary_vars(m::Foo) = (:d, :e)

    ClimaLand.prognostic_types(m::Foo{FT}) where {FT} = (FT, SVector{2, FT})
    ClimaLand.auxiliary_types(m::Foo{FT}) where {FT} = (FT, SVector{2, FT})

    ClimaLand.prognostic_domain_names(m::Foo) = (:surface, :surface)
    ClimaLand.auxiliary_domain_names(m::Foo) = (:surface, :subsurface)

    zmin = FT(1.0)
    zmax = FT(2.0)
    zlim = (zmin, zmax)
    nelements = 5

    column = Column(; zlim = zlim, nelements = nelements)
    m = Foo{FT, typeof(column)}(column)
    Y, p, coords = initialize(m)
    @test Array(parent(Y.foo.a)) == zeros(FT, 1)
    @test Array(parent(Y.foo.b)) == zeros(FT, 2)
    @test Array(parent(p.foo.d)) == zeros(FT, 1)
    @test Array(parent(p.foo.e)) == zeros(FT, 5, 2)
    @testset "Boundary variables default, FT = $FT" begin
        struct BC <: ClimaLand.AbstractBC end
        bc = BC()

        # Test that bc variables added if invoked
        @test boundary_vars(bc, ClimaLand.TopBoundary()) ==
              (:top_bc, :top_bc_wvec)
        @test boundary_vars(bc, ClimaLand.BottomBoundary()) ==
              (:bottom_bc, :bottom_bc_wvec)
        @test boundary_var_domain_names(bc, ClimaLand.TopBoundary()) ==
              (:surface, :surface)
        @test boundary_var_domain_names(bc, ClimaLand.BottomBoundary()) ==
              (:surface, :surface)
        @test boundary_var_types(m, bc, ClimaLand.TopBoundary()) ==
              (FT, ClimaCore.Geometry.WVector{FT})
        @test boundary_var_types(m, bc, ClimaLand.BottomBoundary()) ==
              (FT, ClimaCore.Geometry.WVector{FT})
    end
end

using Test
import ClimaLand:
    name,
    prognostic_types,
    auxiliary_types,
    prognostic_vars,
    auxiliary_vars,
    auxiliary_domain_names,
    prognostic_domain_names,
    prognostic_domain_names,
    lsm_aux_vars,
    lsm_aux_types,
    lsm_aux_domain_names,
    make_compute_exp_tendency,
    make_update_boundary_fluxes,
    make_update_aux,
    make_set_initial_cache,
    land_components
using ClimaLand

for FT in (Float32, Float64)
    @testset "LSM prognostic tests, FT = $FT" begin
        struct DummyModel1{FT} <: ClimaLand.AbstractModel{FT}
            domain::Any
        end
        ClimaLand.name(::DummyModel1) = :m1
        ClimaLand.prognostic_vars(::DummyModel1) = (:a,)
        ClimaLand.prognostic_types(::DummyModel1{FT}) where {FT} = (FT,)
        ClimaLand.prognostic_domain_names(::DummyModel1) = (:surface,)
        ClimaLand.auxiliary_vars(::DummyModel1) = (:b,)
        ClimaLand.auxiliary_types(::DummyModel1{FT}) where {FT} = (FT,)
        ClimaLand.auxiliary_domain_names(::DummyModel1) = (:surface,)
        struct DummyModel2{FT} <: ClimaLand.AbstractModel{FT}
            domain::Any
        end

        ClimaLand.name(::DummyModel2) = :m2
        ClimaLand.prognostic_vars(::DummyModel2) = (:c, :d)
        ClimaLand.prognostic_types(::DummyModel2{FT}) where {FT} = (FT, FT)
        ClimaLand.prognostic_domain_names(::DummyModel2) =
            (:surface, :subsurface)

        struct DummyModel{FT} <: ClimaLand.AbstractLandModel{FT}
            m1::Any
            m2::Any
        end
        ClimaLand.land_components(::DummyModel) = (:m1, :m2)
        ClimaLand.lsm_aux_vars(::DummyModel) = (:i1,)
        ClimaLand.lsm_aux_types(::DummyModel{FT}) where {FT} = (FT,)
        ClimaLand.lsm_aux_domain_names(::DummyModel) = (:surface,)



        function ClimaLand.make_update_aux(::DummyModel1{FT}) where {FT}
            function update_aux!(p, Y, t)
                p.m1.b .= FT(10.0)
            end
            return update_aux!
        end

        function ClimaLand.make_update_boundary_fluxes(
            ::DummyModel{FT},
        ) where {FT}
            function update_aux!(p, Y, t)
                p.i1 .= FT(5.0)
            end
            return update_aux!
        end

        function ClimaLand.make_compute_exp_tendency(
            ::DummyModel1{FT},
        ) where {FT}
            function compute_exp_tendency!(dY, Y, p, t)
                dY.m1.a .= p.m1.b
            end
            return compute_exp_tendency!
        end

        function ClimaLand.make_compute_exp_tendency(
            ::DummyModel2{FT},
        ) where {FT}
            function compute_exp_tendency!(dY, Y, p, t)
                dY.m2.c .= p.i1
                dY.m2.d .= FT(-1.0)
            end
            return compute_exp_tendency!
        end
        d = ClimaLand.Domains.Column(; zlim = (FT(-2), FT(0)), nelements = 5)
        m1 = DummyModel1{FT}(d)
        m2 = DummyModel2{FT}(d)
        m = DummyModel{FT}(m1, m2)
        Y, p, cds = ClimaLand.initialize(m)
        @test propertynames(p) == (:i1, :m1, :m2)
        exp_tendency! = make_exp_tendency(m)
        dY = similar(Y)
        exp_tendency!(dY, Y, p, FT(0))
        @test all(parent(p.i1) .== FT(5))
        @test all(parent(p.m1.b) .== FT(10))
        @test all(parent(dY.m1.a) .== FT(10))
        @test all(parent(dY.m2.c) .== FT(5))
        @test all(parent(dY.m2.d) .== FT(-1))

    end

    @testset "LSM aux tests, FT = $FT" begin
        struct DummyModel3{FT} <: ClimaLand.AbstractModel{FT}
            domain::Any
        end
        ClimaLand.name(::DummyModel3) = :m1
        ClimaLand.auxiliary_vars(::DummyModel3) = (:a, :b)
        ClimaLand.auxiliary_types(::DummyModel3{FT}) where {FT} = (FT, FT)
        ClimaLand.auxiliary_domain_names(::DummyModel3) = (:surface, :surface)
        ClimaLand.prognostic_vars(::DummyModel3) = (:c,)
        ClimaLand.prognostic_types(::DummyModel3{FT}) where {FT} = (FT,)
        ClimaLand.prognostic_domain_names(::DummyModel3) = (:surface,)

        struct DummyModel4{FT} <: ClimaLand.AbstractModel{FT}
            domain::Any
        end

        ClimaLand.name(::DummyModel4) = :m2
        ClimaLand.auxiliary_vars(::DummyModel4) = (:d, :e)
        ClimaLand.auxiliary_types(::DummyModel4{FT}) where {FT} = (FT, FT)
        ClimaLand.auxiliary_domain_names(::DummyModel4) = (:surface, :surface)
        ClimaLand.prognostic_vars(::DummyModel4) = (:f,)
        ClimaLand.prognostic_types(::DummyModel4{FT}) where {FT} = (FT,)
        ClimaLand.prognostic_domain_names(::DummyModel4) = (:surface,)

        struct DummyModelB{FT} <: ClimaLand.AbstractLandModel{FT}
            m1::Any
            m2::Any
        end
        ClimaLand.land_components(::DummyModelB) = (:m1, :m2)

        d = ClimaLand.Domains.Point(; z_sfc = FT(0.0))
        m1 = DummyModel3{FT}(d)
        m2 = DummyModel4{FT}(d)
        m = DummyModelB{FT}(m1, m2)
        Y, p, cds = ClimaLand.initialize(m)
        @test propertynames(p) == (:m1, :m2)
        function ClimaLand.make_update_aux(::DummyModel3{FT}) where {FT}
            function update_aux!(p, Y, t)
                p.m1.b .= FT(10.0)
            end
            return update_aux!
        end

        function ClimaLand.make_update_aux(::DummyModel4{FT}) where {FT}
            function update_aux!(p, Y, t)
                p.m2.d .= FT(10.0)
                p.m2.e .= FT(10.0)
            end
            return update_aux!
        end

        function ClimaLand.make_set_initial_cache(m::DummyModel3{FT}) where {FT}
            update_cache! = ClimaLand.make_update_cache(m)
            function set_initial_cache!(p, Y, t)
                p.m1.a .= FT(2.0)
                update_cache!(p, Y, t)
            end
            return set_initial_cache!
        end

        # The scenario here is that model 1 has a single prescribed but constant
        # variable (a), with another that could get updated each step (b).
        # DummyModel4 has only variables that get updated each step.
        # Test that the land model function properly calls the individual
        # model's functions
        set_initial_cache! = ClimaLand.make_set_initial_cache(m)
        set_initial_cache!(p, Y, FT(0.0))
        @test all(parent(p.m1.a) .== FT(2))
        @test all(parent(p.m1.b) .== FT(10))
        @test all(parent(p.m2.d) .== FT(10))
        @test all(parent(p.m2.e) .== FT(10))
    end
end

using Test
import ClimaLSM:
    name,
    domain_name,
    prognostic_types,
    auxiliary_types,
    prognostic_vars,
    auxiliary_vars,
    make_update_aux,
    make_set_initial_aux_state,
    land_components
using ClimaLSM


@testset "LSM aux tests" begin

    struct DummyModel1{FT} <: ClimaLSM.AbstractModel{FT}
        domain::Any
    end
    ClimaLSM.name(::DummyModel1) = :m1
    ClimaLSM.domain_name(::DummyModel1) = :surface
    ClimaLSM.auxiliary_vars(::DummyModel1) = (:a, :b)
    ClimaLSM.auxiliary_types(::DummyModel1{FT}) where {FT} = (FT, FT)

    struct DummyModel2{FT} <: ClimaLSM.AbstractModel{FT}
        domain::Any
    end

    ClimaLSM.name(::DummyModel2) = :m2
    ClimaLSM.domain_name(::DummyModel2) = :surface
    ClimaLSM.auxiliary_vars(::DummyModel2) = (:c, :d)
    ClimaLSM.auxiliary_types(::DummyModel2{FT}) where {FT} = (FT, FT)

    struct DummyModel{FT} <: ClimaLSM.AbstractLandModel{FT}
        m1::Any
        m2::Any
    end
    ClimaLSM.land_components(::DummyModel) = (:m1, :m2)


    FT = Float32
    d = ClimaLSM.Domains.Point(; z_sfc = FT(0.0))
    m1 = DummyModel1{FT}(d)
    m2 = DummyModel2{FT}(d)
    m = DummyModel{FT}(m1, m2)
    Y, p, cds = ClimaLSM.initialize(m)
    @test propertynames(p) == (:m1, :m2)
    function ClimaLSM.make_update_aux(::DummyModel1{FT}) where {FT}
        function update_aux!(p, Y, t)
            p.m1.b .= FT(10.0)
        end
        return update_aux!
    end

    function ClimaLSM.make_update_aux(::DummyModel2{FT}) where {FT}
        function update_aux!(p, Y, t)
            p.m2.c .= FT(10.0)
            p.m2.d .= FT(10.0)
        end
        return update_aux!
    end

    function ClimaLSM.make_set_initial_aux_state(m::DummyModel1{FT}) where {FT}
        update_aux! = ClimaLSM.make_update_aux(m)
        function set_initial_aux_state!(p, Y, t)
            p.m1.a .= FT(2.0)
            update_aux!(p, Y, t)
        end
        return set_initial_aux_state!
    end

    # The scenario here is that model 1 has a single prescribed but constant
    # variable (a), with another that could get updated each step (b).
    # DummyModel2 has only variables that get updated each step.
    # Test that the land model function properly calls the individual
    # model's functions
    set_initial_aux_state! = ClimaLSM.make_set_initial_aux_state(m)
    set_initial_aux_state!(p, Y, FT(0.0))
    @test all(parent(p.m1.a) .== FT(2))
    @test all(parent(p.m1.b) .== FT(10))
    @test all(parent(p.m2.c) .== FT(10))
    @test all(parent(p.m2.d) .== FT(10))

end

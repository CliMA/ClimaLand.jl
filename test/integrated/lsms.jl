using Test
import ClimaLSM:
    name,
    prognostic_types,
    auxiliary_types,
    prognostic_vars,
    auxiliary_vars,
    auxiliary_domain_names,
    prognostic_domain_names,
    prognostic_domain_names,
    interaction_vars,
    interaction_types,
    interaction_domain_names,
    make_compute_exp_tendency,
    make_interactions_update_aux,
    make_update_aux,
    make_set_initial_aux_state,
    land_components
using ClimaLSM

@testset "LSM prognostic tests" begin

    struct DummyModel1{FT} <: ClimaLSM.AbstractModel{FT}
        domain::Any
    end
    ClimaLSM.name(::DummyModel1) = :m1
    ClimaLSM.prognostic_vars(::DummyModel1) = (:a,)
    ClimaLSM.prognostic_types(::DummyModel1{FT}) where {FT} = (FT,)
    ClimaLSM.prognostic_domain_names(::DummyModel1) = (:surface,)
    ClimaLSM.auxiliary_vars(::DummyModel1) = (:b,)
    ClimaLSM.auxiliary_types(::DummyModel1{FT}) where {FT} = (FT,)
    ClimaLSM.auxiliary_domain_names(::DummyModel1) = (:surface,)
    struct DummyModel2{FT} <: ClimaLSM.AbstractModel{FT}
        domain::Any
    end

    ClimaLSM.name(::DummyModel2) = :m2
    ClimaLSM.prognostic_vars(::DummyModel2) = (:c, :d)
    ClimaLSM.prognostic_types(::DummyModel2{FT}) where {FT} = (FT, FT)
    ClimaLSM.prognostic_domain_names(::DummyModel2) = (:surface, :subsurface)

    struct DummyModel{FT} <: ClimaLSM.AbstractLandModel{FT}
        m1::Any
        m2::Any
    end
    ClimaLSM.land_components(::DummyModel) = (:m1, :m2)
    ClimaLSM.interaction_vars(::DummyModel) = (:i1,)
    ClimaLSM.interaction_types(::DummyModel{FT}) where {FT} = (FT,)
    ClimaLSM.interaction_domain_names(::DummyModel) = (:surface,)



    function ClimaLSM.make_update_aux(::DummyModel1{FT}) where {FT}
        function update_aux!(p, Y, t)
            p.m1.b .= FT(10.0)
        end
        return update_aux!
    end

    function ClimaLSM.make_interactions_update_aux(::DummyModel{FT}) where {FT}
        function update_aux!(p, Y, t)
            p.i1 .= FT(5.0)
        end
        return update_aux!
    end

    function ClimaLSM.make_compute_exp_tendency(::DummyModel1{FT}) where {FT}
        function compute_exp_tendency!(dY, Y, p, t)
            dY.m1.a .= p.m1.b
        end
        return compute_exp_tendency!
    end

    function ClimaLSM.make_compute_exp_tendency(::DummyModel2{FT}) where {FT}
        function compute_exp_tendency!(dY, Y, p, t)
            dY.m2.c .= p.i1
            dY.m2.d .= FT(-1.0)
        end
        return compute_exp_tendency!
    end
    FT = Float32
    d = ClimaLSM.Domains.Column(; zlim = (FT(-2), FT(0)), nelements = 5)
    m1 = DummyModel1{FT}(d)
    m2 = DummyModel2{FT}(d)
    m = DummyModel{FT}(m1, m2)
    Y, p, cds = ClimaLSM.initialize(m)
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


@testset "LSM aux tests" begin

    struct DummyModel3{FT} <: ClimaLSM.AbstractModel{FT}
        domain::Any
    end
    ClimaLSM.name(::DummyModel3) = :m1
    ClimaLSM.auxiliary_vars(::DummyModel3) = (:a, :b)
    ClimaLSM.auxiliary_types(::DummyModel3{FT}) where {FT} = (FT, FT)
    ClimaLSM.auxiliary_domain_names(::DummyModel3) = (:surface, :surface)

    struct DummyModel4{FT} <: ClimaLSM.AbstractModel{FT}
        domain::Any
    end

    ClimaLSM.name(::DummyModel4) = :m2
    ClimaLSM.auxiliary_vars(::DummyModel4) = (:c, :d)
    ClimaLSM.auxiliary_types(::DummyModel4{FT}) where {FT} = (FT, FT)
    ClimaLSM.auxiliary_domain_names(::DummyModel4) = (:surface, :surface)

    struct DummyModelB{FT} <: ClimaLSM.AbstractLandModel{FT}
        m1::Any
        m2::Any
    end
    ClimaLSM.land_components(::DummyModelB) = (:m1, :m2)


    FT = Float32
    d = ClimaLSM.Domains.Point(; z_sfc = FT(0.0))
    m1 = DummyModel3{FT}(d)
    m2 = DummyModel4{FT}(d)
    m = DummyModelB{FT}(m1, m2)
    Y, p, cds = ClimaLSM.initialize(m)
    @test propertynames(p) == (:m1, :m2)
    function ClimaLSM.make_update_aux(::DummyModel3{FT}) where {FT}
        function update_aux!(p, Y, t)
            p.m1.b .= FT(10.0)
        end
        return update_aux!
    end

    function ClimaLSM.make_update_aux(::DummyModel4{FT}) where {FT}
        function update_aux!(p, Y, t)
            p.m2.c .= FT(10.0)
            p.m2.d .= FT(10.0)
        end
        return update_aux!
    end

    function ClimaLSM.make_set_initial_aux_state(m::DummyModel3{FT}) where {FT}
        update_aux! = ClimaLSM.make_update_aux(m)
        function set_initial_aux_state!(p, Y, t)
            p.m1.a .= FT(2.0)
            update_aux!(p, Y, t)
        end
        return set_initial_aux_state!
    end

    # The scenario here is that model 1 has a single prescribed but constant
    # variable (a), with another that could get updated each step (b).
    # DummyModel4 has only variables that get updated each step.
    # Test that the land model function properly calls the individual
    # model's functions
    set_initial_aux_state! = ClimaLSM.make_set_initial_aux_state(m)
    set_initial_aux_state!(p, Y, FT(0.0))
    @test all(parent(p.m1.a) .== FT(2))
    @test all(parent(p.m1.b) .== FT(10))
    @test all(parent(p.m2.c) .== FT(10))
    @test all(parent(p.m2.d) .== FT(10))

end

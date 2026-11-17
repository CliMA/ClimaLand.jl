using Test
import ClimaComms
ClimaComms.@import_required_backends
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
    land_components,
    get_model_callbacks
using ClimaLand
import ClimaCore
import ClimaLand.Parameters as LP


@testset "LSM prognostic tests" begin
    struct DummyModel1{FT} <: ClimaLand.AbstractModel{FT}
        domain::Any
    end
    struct DummyModel2{FT} <: ClimaLand.AbstractModel{FT}
        domain::Any
    end
    struct DummyModel{FT} <: ClimaLand.AbstractLandModel{FT}
        m1::Any
        m2::Any
    end
    for FT in (Float32, Float64)

        ClimaLand.name(::DummyModel1{FT}) = :m1
        ClimaLand.prognostic_vars(::DummyModel1{FT}) = (:a,)
        ClimaLand.prognostic_types(::DummyModel1{FT}) = (FT,)
        ClimaLand.prognostic_domain_names(::DummyModel1{FT}) = (:surface,)
        ClimaLand.auxiliary_vars(::DummyModel1{FT}) = (:b,)
        ClimaLand.auxiliary_types(::DummyModel1{FT}) = (FT,)
        ClimaLand.auxiliary_domain_names(::DummyModel1{FT}) = (:surface,)
        ClimaLand.get_model_callbacks(::DummyModel1) = (:foo,)# this should be a callback, but we are just testing that is is added correctly
        ClimaLand.name(::DummyModel2{FT}) = :m2
        ClimaLand.prognostic_vars(::DummyModel2{FT}) = (:c, :d)
        ClimaLand.prognostic_types(::DummyModel2{FT}) = (FT, FT)
        ClimaLand.prognostic_domain_names(::DummyModel2{FT}) =
            (:surface, :subsurface)

        ClimaLand.land_components(::DummyModel{FT}) = (:m1, :m2)
        ClimaLand.lsm_aux_vars(::DummyModel{FT}) = (:i1,)
        ClimaLand.lsm_aux_types(::DummyModel{FT}) = (FT,)
        ClimaLand.lsm_aux_domain_names(::DummyModel{FT}) = (:surface,)

        function ClimaLand.make_update_aux(::DummyModel1{FT})
            function update_aux!(p, Y, t)
                p.m1.b .= FT(10.0)
            end
            return update_aux!
        end

        function ClimaLand.make_update_boundary_fluxes(::DummyModel{FT})
            function update_aux!(p, Y, t)
                p.i1 .= FT(5.0)
            end
            return update_aux!
        end

        function ClimaLand.make_compute_exp_tendency(::DummyModel1{FT})
            function compute_exp_tendency!(dY, Y, p, t)
                dY.m1.a .= p.m1.b
            end
            return compute_exp_tendency!
        end

        function ClimaLand.make_compute_exp_tendency(::DummyModel2{FT})
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
        @test get_model_callbacks(m) == (:foo,)
        @test get_model_callbacks(m2) == ()
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
end

@testset "LSM aux tests" begin
    struct DummyModel3{FT} <: ClimaLand.AbstractModel{FT}
        domain::Any
    end
    struct DummyModel4{FT} <: ClimaLand.AbstractModel{FT}
        domain::Any
    end
    struct DummyModelB{FT} <: ClimaLand.AbstractLandModel{FT}
        m1::Any
        m2::Any
    end
    for FT in (Float32, Float64)
        ClimaLand.name(::DummyModel3{FT}) = :m1
        ClimaLand.auxiliary_vars(::DummyModel3{FT}) = (:a, :b)
        ClimaLand.auxiliary_types(::DummyModel3{FT}) = (FT, FT)
        ClimaLand.auxiliary_domain_names(::DummyModel3{FT}) =
            (:surface, :surface)
        ClimaLand.prognostic_vars(::DummyModel3{FT}) = (:c,)
        ClimaLand.prognostic_types(::DummyModel3{FT}) = (FT,)
        ClimaLand.prognostic_domain_names(::DummyModel3{FT}) = (:surface,)
        ClimaLand.name(::DummyModel4{FT}) = :m2
        ClimaLand.auxiliary_vars(::DummyModel4{FT}) = (:d, :e)
        ClimaLand.auxiliary_types(::DummyModel4{FT}) = (FT, FT)
        ClimaLand.auxiliary_domain_names(::DummyModel4{FT}) =
            (:surface, :surface)
        ClimaLand.prognostic_vars(::DummyModel4{FT}) = (:f,)
        ClimaLand.prognostic_types(::DummyModel4{FT}) = (FT,)
        ClimaLand.prognostic_domain_names(::DummyModel4{FT}) = (:surface,)
        ClimaLand.land_components(::DummyModelB{FT}) = (:m1, :m2)

        d = ClimaLand.Domains.Point(; z_sfc = FT(0.0))
        m1 = DummyModel3{FT}(d)
        m2 = DummyModel4{FT}(d)
        m = DummyModelB{FT}(m1, m2)
        Y, p, cds = ClimaLand.initialize(m)
        @test propertynames(p) == (:m1, :m2)

        function ClimaLand.make_update_aux(::DummyModel3{FT})
            function update_aux!(p, Y, t)
                p.m1.a .= FT(2.0)
                p.m1.b .= FT(10.0)
            end
            return update_aux!
        end

        function ClimaLand.make_update_aux(::DummyModel4{FT})
            function update_aux!(p, Y, t)
                p.m2.d .= FT(10.0)
                p.m2.e .= FT(10.0)
            end
            return update_aux!
        end

        set_initial_cache! = ClimaLand.make_set_initial_cache(m)
        set_initial_cache!(p, Y, FT(0.0))
        @test all(parent(p.m1.a) .== FT(2))
        @test all(parent(p.m1.b) .== FT(10))
        @test all(parent(p.m2.d) .== FT(10))
        @test all(parent(p.m2.e) .== FT(10))
    end
end

@testset "Litter under vegetation" begin
    for FT in (Float32, Float64)
        toml_dict = LP.create_toml_dict(FT)
        domain = ClimaLand.Domains.Plane(;
            xlim = FT.((0, 3)),
            ylim = FT.((0, 1)),
            nelements = (3, 1),
            periodic = (true, true),
        )
        surface_space = domain.space.surface
        x = ClimaCore.Fields.coordinate_field(surface_space).x
        # bare (x < 1), grassland (1 < x < 2) and forest (x > 2) points
        max_lai = @. ifelse(x < 1, FT(0), FT(4))
        tree_share = @. ifelse(x < 2, FT(0), FT(1))
        litter = ClimaLand.vegetation_litter(
            surface_space,
            toml_dict;
            tree_share,
            max_lai,
        )
        @test litter isa ClimaLand.Soil.SoilLitter{FT}
        @test litter.κ == toml_dict["litter_thermal_conductivity"]
        expected = @. ifelse(
            x < 1,
            FT(0),
            ifelse(
                x < 2,
                toml_dict["litter_thickness_grass"],
                toml_dict["litter_thickness_tree"],
            ),
        )
        @test parent(litter.thickness) == parent(expected)
    end
end

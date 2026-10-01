using Test
using Dates
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.Domains: Point, Plane
import ClimaLand.Parameters as LP

"""
    set_state!(Y, p, canopy, t0; kwargs...)

Sets a non-trivial canopy state, with carbon pools and climate means overridden by
`kwargs`, and updates the cache.
"""
function set_state!(Y, p, canopy, t0; kwargs...)
    FT = eltype(Y.canopy.energy.T)
    Y.canopy.hydraulics.ϑ_l .= canopy.hydraulics.parameters.ν / 2
    Y.canopy.energy.T .= FT(290.5)
    if canopy.biomass isa Canopy.PrognosticCarbonModel
        state = (;
            C_sugar = 0.3,
            C_leaf = 0.2,
            C_stem = 5,
            C_root = 1,
            T_annual = 290,
            P_annual = 1,
            kwargs...,
        )
        for name in keys(state)
            getproperty(Y.canopy.biomass, name) .= FT(state[name])
        end
    end
    set_initial_cache! = make_set_initial_cache(canopy)
    set_initial_cache!(p, Y, t0)
    # Acclimated P-model capacities, so that GPP and Rd are not zero
    ClimaLand.Simulations.set_canopy_component_initial_conditions!(
        Y,
        p,
        canopy.photosynthesis,
        canopy,
    )
    set_initial_cache!(p, Y, t0)
    return nothing
end

point(x) = first(Array(parent(x)))

for FT in (Float32, Float64)
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    longlat = (FT(-118.14), FT(34.15))
    pt = Point(; z_sfc = FT(0), longlat)
    plane = Plane(;
        xlim = (FT(0), FT(10)),
        ylim = (FT(0), FT(10)),
        nelements = (2, 2),
        longlat,
    )
    t0 = 0.0
    start_date = DateTime(2005)

    radiation = PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput((t) -> eltype(t)(300)),
        TimeVaryingInput((t) -> eltype(t)(300)),
        start_date;
        cosθs = (t, s) -> default_cos_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            latitude = FT(40),
            longitude = FT(-120),
        ),
        toml_dict,
    )
    precip = TimeVaryingInput((t) -> eltype(t)(0))
    atmos = ClimaLand.PrescribedAtmosphere(
        precip,
        precip,
        TimeVaryingInput((t) -> eltype(t)(290)),
        TimeVaryingInput((t) -> eltype(t)(2)),
        TimeVaryingInput((t) -> eltype(t)(0.011)),
        TimeVaryingInput((t) -> eltype(t)(101325)),
        start_date,
        FT(3),
        toml_dict,
    )
    forcing = (; atmos, radiation, ground = PrescribedGroundConditions{FT}())
    LAI = TimeVaryingInput((t) -> FT(3))
    lai_model = Canopy.PrescribedBiomassModel{FT}(;
        LAI,
        SAI = FT(1),
        RAI = FT(1),
        rooting_depth = FT(0.5),
        height = FT(2),
    )
    biomass = Canopy.PrognosticCarbonModel{FT}(lai_model, toml_dict)
    (; M_C, a) = biomass.parameters

    canopy = Canopy.CanopyModel{FT}(pt, forcing, LAI, toml_dict; biomass)
    Y, p, _ = initialize(canopy)
    dY = similar(Y)
    exp_tendency! = make_compute_exp_tendency(canopy)

    @testset "Carbon balance of the pools, FT = $FT" begin
        @test canopy.autotrophic_respiration isa
              Canopy.PoolBasedAutotrophicRespirationModel
        canopy_plane =
            Canopy.CanopyModel{FT}(plane, forcing, LAI, toml_dict; biomass)
        Y_plane, p_plane, _ = initialize(canopy_plane)
        for (model, Y_i, p_i) in
            ((canopy, Y, p), (canopy_plane, Y_plane, p_plane))
            set_state!(Y_i, p_i, model, t0)
            dY_i = similar(Y_i)
            make_compute_exp_tendency(model)(dY_i, Y_i, p_i, t0)
            (; Rm, Rg, Ra, S, L_leaf, L_stem, L_root) =
                p_i.canopy.biomass.carbon
            GPP = M_C .* Canopy.get_GPP(p_i, model.photosynthesis)
            dC = @. dY_i.canopy.biomass.C_sugar +
               dY_i.canopy.biomass.C_leaf +
               dY_i.canopy.biomass.C_stem +
               dY_i.canopy.biomass.C_root
            residual = @. dC - (GPP - Ra - (L_leaf + L_stem + L_root))
            scale = maximum(abs, parent(GPP)) + maximum(abs, parent(dC))
            @test maximum(abs, parent(residual)) <= 10 * eps(FT) * scale
            @test all(parent(GPP) .> 0)
            @test all(parent(S) .> 0)
            @test parent(Ra) ≈ parent(Rm) .+ parent(Rg)
            @test parent(Rg) ≈ (1 - a) .* parent(S)
            @test parent(p_i.canopy.autotrophic_respiration.Ra) ≈
                  parent(Ra) ./ M_C
            @test all(parent(p_i.canopy.biomass.cVeg) .≈ FT(6.5))
            @test ClimaLand.initialize_jacobian(Y_i) isa Any
        end
    end

    @testset "GPP and LAI do not depend on the pools, FT = $FT" begin
        canopy_without = Canopy.CanopyModel{FT}(
            pt,
            forcing,
            LAI,
            toml_dict;
            biomass = lai_model,
        )
        Y_without, p_without, _ = initialize(canopy_without)
        set_state!(Y_without, p_without, canopy_without, t0)
        set_state!(Y, p, canopy, t0)
        for name in (:leaf, :stem, :root)
            @test getproperty(p_without.canopy.biomass.area_index, name) ==
                  getproperty(p.canopy.biomass.area_index, name)
        end
        @test Canopy.get_GPP(p_without, canopy.photosynthesis) ==
              Canopy.get_GPP(p, canopy.photosynthesis)
    end

    @testset "Maintenance respiration, FT = $FT" begin
        Rm() = point(p.canopy.biomass.carbon.Rm)
        # Only the sapwood and roots have a Q10: the leaf term is Rd itself.
        set_state!(Y, p, canopy, t0; C_sugar = 1e3, C_stem = 0, C_root = 0)
        Rd = point(Canopy.get_Rd_canopy(p, canopy.photosynthesis))
        @test Rd > 0
        @test Rm() ≈ M_C * Rd

        # A healthy sugar pool barely limits respiration; an empty one stops it
        set_state!(Y, p, canopy, t0; C_sugar = 1e3)
        R_ample = Rm()
        set_state!(Y, p, canopy, t0; C_sugar = 0.1)
        @test Rm() ≈ R_ample rtol = 1e-3
        set_state!(Y, p, canopy, t0; C_sugar = 0)
        @test Rm() == 0
        @test all(parent(p.canopy.autotrophic_respiration.Ra) .== 0)
        exp_tendency!(dY, Y, p, t0)
        @test all(parent(dY.canopy.biomass.C_sugar) .>= 0)
    end

    @testset "Climate dependence of the stem pool, FT = $FT" begin
        (; map_half_woody, n_map_woody, q_τ_stem, T_ref_τ_stem) =
            biomass.parameters
        w(MAP) = Canopy.woody_fraction(MAP, map_half_woody, n_map_woody)
        @test w(map_half_woody) ≈ FT(0.5)
        @test w(FT(0)) == 0
        @test w(FT(0.2)) < w(FT(1)) < w(FT(3)) < 1
        @test Canopy.woody_fraction(FT(0.2), FT(0), n_map_woody) == 1

        s(MAT) = Canopy.tau_stem_scale(MAT, T_ref_τ_stem, q_τ_stem)
        @test s(T_ref_τ_stem) == 1
        @test s(T_ref_τ_stem + 10) == 1
        @test s(T_ref_τ_stem - 20) ≈ q_τ_stem^2
        @test s(FT(0)) == Canopy.MAX_TAU_STEM_SCALE
        @test Canopy.tau_stem_scale(FT(250), T_ref_τ_stem, FT(1)) == 1

        # Without precipitation, allocation withheld from the stem goes to roots
        set_state!(Y, p, canopy, t0; P_annual = 0)
        exp_tendency!(dY, Y, p, t0)
        (; S, L_stem, L_root) = p.canopy.biomass.carbon
        (; f_leaf_c3, f_leaf_c4) = biomass.parameters
        fc3 = Canopy.get_fractional_c3(p, canopy)
        f_leaf = @. Canopy.blend(f_leaf_c3, f_leaf_c4, fc3)
        @test dY.canopy.biomass.C_stem == @. -L_stem
        @test parent(dY.canopy.biomass.C_root) ≈
              parent(@. a * (1 - f_leaf) * S - L_root)
    end

    @testset "Equilibrium pools, FT = $FT" begin
        # Drivers at a point, which do not depend on the pools
        set_state!(Y, p, canopy, t0)
        GPP = point(Canopy.get_GPP(p, canopy.photosynthesis))
        Rd = point(Canopy.get_Rd_canopy(p, canopy.photosynthesis))
        fc3 = point(Canopy.get_fractional_c3(p, canopy))
        (; Q10, T_ref) = biomass.parameters
        f_T = Q10^((point(Y.canopy.energy.T) - T_ref) / 10)
        MAT = point(Y.canopy.biomass.T_annual)
        MAP = point(Y.canopy.biomass.P_annual)
        eq = Canopy.equilibrium_carbon_pools(
            biomass.parameters,
            GPP,
            Rd,
            f_T,
            MAT,
            MAP,
            fc3,
        )
        @test eq.C_stem > eq.C_root > eq.C_leaf > 0

        # With the pools at equilibrium, and ample sugar, the model's own fluxes
        # balance the sugar pool and the stem pool
        set_state!(Y, p, canopy, t0; C_sugar = 1e3, eq...)
        (; τ_leaf, f_leaf_c3, f_leaf_c4, f_stem_c3, f_stem_c4) =
            biomass.parameters
        (; map_half_woody, n_map_woody) = biomass.parameters
        S = eq.C_leaf / (a * Canopy.blend(f_leaf_c3, f_leaf_c4, fc3) * τ_leaf)
        Rm = point(p.canopy.biomass.carbon.Rm)
        @test M_C * GPP - Rm ≈ S rtol = sqrt(eps(FT))
        f_stem =
            Canopy.blend(f_stem_c3, f_stem_c4, fc3) *
            Canopy.woody_fraction(MAP, map_half_woody, n_map_woody)
        @test a * f_stem * S ≈ point(p.canopy.biomass.carbon.L_stem) rtol =
            sqrt(eps(FT))

        # No pools where leaf respiration exceeds GPP
        @test Canopy.equilibrium_carbon_pools(
            biomass.parameters,
            Rd,
            GPP + Rd,
            f_T,
            MAT,
            MAP,
            fc3,
        ).C_stem == 0
    end

    @testset "cveg diagnostic reads the pools, FT = $FT" begin
        set_state!(Y, p, canopy, t0)
        out = ClimaCore.Fields.zeros(canopy.domain.space.surface)
        ClimaLand.Diagnostics.compute_vegetation_carbon!(out, Y, p, t0, canopy)
        @test out == p.canopy.biomass.cVeg
    end

    @testset "Initial conditions, FT = $FT" begin
        set_state!(Y, p, canopy, t0)
        ClimaLand.Simulations.set_canopy_component_initial_conditions!(
            Y,
            p,
            canopy.biomass,
            canopy,
        )
        for pool in (:C_sugar, :C_leaf, :C_stem, :C_root)
            @test all(iszero, parent(getproperty(Y.canopy.biomass, pool)))
        end
        @test Y.canopy.biomass.T_annual == p.drivers.T
        # The climatological precipitation of the site, in m yr^-1
        @test all(FT(0.1) .< parent(Y.canopy.biomass.P_annual) .< FT(2))
    end

    @testset "Component compatibility, FT = $FT" begin
        @test_throws AssertionError Canopy.CanopyModel{FT}(
            pt,
            forcing,
            LAI,
            toml_dict;
            biomass,
            autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT}(
                toml_dict,
            ),
        )
        @test_throws AssertionError Canopy.CanopyModel{FT}(
            pt,
            forcing,
            LAI,
            toml_dict;
            autotrophic_respiration = Canopy.PoolBasedAutotrophicRespirationModel{
                FT,
            }(),
        )
    end

    @testset "Wrapping ZhouOptimalLAIModel, FT = $FT" begin
        zhou = Canopy.ZhouOptimalLAIModel{FT}(pt, toml_dict)
        wrapped = Canopy.PrognosticCarbonModel{FT}(zhou, toml_dict)
        @test issubset(
            ClimaLand.prognostic_vars(zhou),
            ClimaLand.prognostic_vars(wrapped),
        )
        canopy_z =
            Canopy.CanopyModel{FT}(pt, forcing, toml_dict; biomass = zhou)
        canopy_w =
            Canopy.CanopyModel{FT}(pt, forcing, toml_dict; biomass = wrapped)
        Y_z, p_z, _ = initialize(canopy_z)
        Y_w, p_w, _ = initialize(canopy_w)
        for (canopy, Y, p) in ((canopy_z, Y_z, p_z), (canopy_w, Y_w, p_w))
            set_state!(Y, p, canopy, t0)
            ClimaLand.Simulations.set_canopy_component_initial_conditions!(
                Y,
                p,
                canopy.biomass,
                canopy,
            )
            make_set_initial_cache(canopy)(p, Y, t0)
        end
        @test p_z.canopy.biomass.area_index.leaf ==
              p_w.canopy.biomass.area_index.leaf
        @test Canopy.get_GPP(p_z, canopy_z.photosynthesis) ==
              Canopy.get_GPP(p_w, canopy_w.photosynthesis)
        fc3 = Canopy.get_fractional_c3(p_w, canopy_w)
        @test all(parent(fc3 .== Canopy.get_fractional_c3(p_z, canopy_z)))
        @test all(parent(fc3 .!= 1))
        ρ_m_liq = LP.ρ_m_liq(earth_param_set)
        @test parent(Y_w.canopy.biomass.P_annual) ≈
              parent(Y_z.canopy.biomass.precip_annual) ./ ρ_m_liq

        diags_z, diags_w = String[], String[]
        ClimaLand.Diagnostics.add_diagnostics!(diags_z, canopy_z, zhou)
        ClimaLand.Diagnostics.add_diagnostics!(diags_w, canopy_w, wrapped)
        @test !isempty(diags_z)
        @test diags_w == diags_z
    end
end

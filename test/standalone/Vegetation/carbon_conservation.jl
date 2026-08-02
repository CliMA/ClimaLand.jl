using Test
using Dates
using StaticArrays
using Insolation
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.Domains: Point, Plane

import ClimaLand
import ClimaLand.Parameters as LP

# Carbon conservation for the prognostic pools. Summing the four pool tendencies
# must reproduce d(ΣC)/dt = GPP - Ra - litter exactly, where Ra = Rm + Rg is the
# carbon model's own autotrophic respiration and litter is the sum of the three
# turnover fluxes. If allocation, growth respiration or turnover is
# double-counted or dropped, this residual is what catches it.
for FT in (Float32, Float64)
    cmax = FT(10)
    cmin = FT(0)
    nelems = 10
    longlat = (FT(-118.14), FT(34.15))
    pt = Point(; z_sfc = FT(0), longlat)
    plane = Plane(;
        xlim = (cmin, cmax),
        ylim = (cmin, cmax),
        nelements = (nelems, nelems),
        longlat,
    )
    domains = [pt, plane]

    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)

    t0 = 0.0
    start_date = DateTime(2005)

    SW_d = TimeVaryingInput((t) -> eltype(t)(300.0))
    LW_d = TimeVaryingInput((t) -> eltype(t)(300.0))
    cos_zenith_angle =
        (t, s) -> default_cos_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            latitude = FT(40.0),
            longitude = FT(-120.0),
        )
    radiation = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        cosθs = cos_zenith_angle,
        toml_dict = toml_dict,
    )

    precip = TimeVaryingInput((t) -> eltype(t)(0))
    T_atmos = TimeVaryingInput((t) -> eltype(t)(290.0))
    u_atmos = TimeVaryingInput((t) -> eltype(t)(2.0))
    q_atmos = TimeVaryingInput((t) -> eltype(t)(0.011))
    h_atmos = FT(3)
    P_atmos = TimeVaryingInput((t) -> eltype(t)(101325))
    atmos = ClimaLand.PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        toml_dict,
    )
    ground = PrescribedGroundConditions{FT}()

    LAI_value = FT(3)
    LAI = TimeVaryingInput(t -> LAI_value)
    lai_model = Canopy.PrescribedBiomassModel{FT}(;
        LAI,
        SAI = FT(1),
        RAI = FT(1),
        rooting_depth = FT(0.5),
        height = FT(2),
    )
    biomass = Canopy.PrognosticCarbonModel{FT}(lai_model, toml_dict)

    @testset "Carbon pool conservation, FT = $FT" begin
        for domain in domains
            hydraulics = Canopy.PlantHydraulicsModel{FT}(domain, toml_dict;)
            canopy = ClimaLand.Canopy.CanopyModel{FT}(
                domain,
                (; radiation, atmos, ground),
                LAI,
                toml_dict;
                hydraulics,
                biomass,
            )

            Y, p, cds = initialize(canopy)
            Y.canopy.hydraulics.ϑ_l .= canopy.hydraulics.parameters.ν / 2
            Y.canopy.energy.T .= FT(290.5)
            # Non-trivial pools, so allocation, respiration and turnover are all
            # active. A zero state would satisfy the balance trivially.
            Y.canopy.biomass.C_sugar .= FT(0.3)
            Y.canopy.biomass.C_leaf .= FT(0.2)
            Y.canopy.biomass.C_stem .= FT(5.0)
            Y.canopy.biomass.C_root .= FT(1.0)

            set_initial_cache! = make_set_initial_cache(canopy)
            set_initial_cache!(p, Y, t0)

            dY = similar(Y)
            exp_tendency! = make_compute_exp_tendency(canopy)
            exp_tendency!(dY, Y, p, t0)

            M_C = FT(0.012011)
            GPP_mol = Canopy.get_GPP(p, canopy.photosynthesis)
            GPP = @. M_C * GPP_mol
            dsum = @. dY.canopy.biomass.C_sugar +
               dY.canopy.biomass.C_leaf +
               dY.canopy.biomass.C_stem +
               dY.canopy.biomass.C_root
            litter = @. p.canopy.biomass.carbon.L_leaf +
               p.canopy.biomass.carbon.L_stem +
               p.canopy.biomass.carbon.L_root
            residual = @. dsum - (GPP - p.canopy.biomass.carbon.Ra - litter)

            # The pool fluxes are O(1e-8) kg C m^-2 s^-1, so an absolute
            # tolerance scaled to the largest term is the meaningful check.
            scale = maximum(abs.(parent(GPP))) + maximum(abs.(parent(dsum)))
            @test maximum(abs.(parent(residual))) <=
                  10 * eps(FT) * max(scale, eps(FT))

            # Ra must be exactly the sum of its two parts.
            @test all(
                parent(
                    @. abs(
                        p.canopy.biomass.carbon.Ra - (
                            p.canopy.biomass.carbon.Rm +
                            p.canopy.biomass.carbon.Rg
                        ),
                    )
                ) .<= eps(FT),
            )

            # Growth respiration is the (1-a) share of what allocation draws.
            a = canopy.biomass.parameters.a
            @test all(
                parent(
                    @. abs(
                        p.canopy.biomass.carbon.Rg -
                        (1 - a) * p.canopy.biomass.carbon.S,
                    )
                ) .<= 10 * eps(FT),
            )

            # cVeg is the sum of the four pools.
            @test all(
                parent(
                    @. abs(
                        p.canopy.biomass.cVeg - (
                            Y.canopy.biomass.C_sugar +
                            Y.canopy.biomass.C_leaf +
                            Y.canopy.biomass.C_stem +
                            Y.canopy.biomass.C_root
                        ),
                    )
                ) .<= 10 * eps(FT),
            )
        end
    end

    # Rule 1: adding the carbon model must not change LAI or GPP. The pools are
    # a pure follower in phase 1, so a run with them must reproduce the area
    # indices and photosynthesis of a run without them.
    @testset "Carbon model does not change LAI or GPP, FT = $FT" begin
        for domain in domains
            hydraulics = Canopy.PlantHydraulicsModel{FT}(domain, toml_dict;)
            canopy_without = ClimaLand.Canopy.CanopyModel{FT}(
                domain,
                (; radiation, atmos, ground),
                LAI,
                toml_dict;
                hydraulics,
                biomass = lai_model,
            )
            canopy_with = ClimaLand.Canopy.CanopyModel{FT}(
                domain,
                (; radiation, atmos, ground),
                LAI,
                toml_dict;
                hydraulics,
                biomass,
            )

            states = map((canopy_without, canopy_with)) do canopy
                Y, p, _ = initialize(canopy)
                Y.canopy.hydraulics.ϑ_l .=
                    canopy.hydraulics.parameters.ν / 2
                Y.canopy.energy.T .= FT(290.5)
                if canopy.biomass isa Canopy.PrognosticCarbonModel
                    Y.canopy.biomass.C_sugar .= FT(0.3)
                    Y.canopy.biomass.C_leaf .= FT(0.2)
                    Y.canopy.biomass.C_stem .= FT(5.0)
                    Y.canopy.biomass.C_root .= FT(1.0)
                end
                make_set_initial_cache(canopy)(p, Y, t0)
                (Y, p, canopy)
            end
            (_, p_without, canopy_a) = states[1]
            (_, p_with, canopy_b) = states[2]

            @test parent(p_without.canopy.biomass.area_index.leaf) ==
                  parent(p_with.canopy.biomass.area_index.leaf)
            @test parent(p_without.canopy.biomass.area_index.stem) ==
                  parent(p_with.canopy.biomass.area_index.stem)
            @test parent(p_without.canopy.biomass.area_index.root) ==
                  parent(p_with.canopy.biomass.area_index.root)
            @test parent(Canopy.get_GPP(p_without, canopy_a.photosynthesis)) ==
                  parent(Canopy.get_GPP(p_with, canopy_b.photosynthesis))
            @test parent(p_without.canopy.autotrophic_respiration.Ra) ==
                  parent(p_with.canopy.autotrophic_respiration.Ra)
        end
    end

    # An exhausted sugar pool must not be driven negative by maintenance
    # respiration. The 4-site battery hit exactly this at the Sahara point,
    # where GPP is zero all year: Rm kept drawing on an empty pool and C_sugar
    # reached -0.075 kg C m^-2.
    @testset "Maintenance respiration stops at an empty sugar pool, FT = $FT" begin
        hydraulics = Canopy.PlantHydraulicsModel{FT}(pt, toml_dict;)
        canopy = ClimaLand.Canopy.CanopyModel{FT}(
            pt,
            (; radiation, atmos, ground),
            LAI,
            toml_dict;
            hydraulics,
            biomass,
        )
        Y, p, _ = initialize(canopy)
        Y.canopy.hydraulics.ϑ_l .= canopy.hydraulics.parameters.ν / 2
        Y.canopy.energy.T .= FT(290.5)
        Y.canopy.biomass.C_leaf .= FT(0.2)
        Y.canopy.biomass.C_stem .= FT(5.0)
        Y.canopy.biomass.C_root .= FT(1.0)
        make_set_initial_cache(canopy)(p, Y, t0)
        dY = similar(Y)
        exp_tendency! = make_compute_exp_tendency(canopy)

        # Ample sugar: respiration runs essentially unthrottled.
        Y.canopy.biomass.C_sugar .= FT(0.3)
        make_set_initial_cache(canopy)(p, Y, t0)
        exp_tendency!(dY, Y, p, t0)
        Rm_ample = first(Array(parent(p.canopy.biomass.carbon.Rm)))
        @test Rm_ample > 0

        # Empty sugar: respiration must vanish, so the pool cannot go negative.
        Y.canopy.biomass.C_sugar .= FT(0)
        make_set_initial_cache(canopy)(p, Y, t0)
        exp_tendency!(dY, Y, p, t0)
        @test first(Array(parent(p.canopy.biomass.carbon.Rm))) == 0
        @test all(parent(dY.canopy.biomass.C_sugar) .>= 0)

        # The throttle must not meaningfully touch a healthy pool: Rm at a
        # normal sugar level must match the unthrottled limit to within 0.1%.
        Y.canopy.biomass.C_sugar .= FT(1e6)
        make_set_initial_cache(canopy)(p, Y, t0)
        exp_tendency!(dY, Y, p, t0)
        Rm_unthrottled = first(Array(parent(p.canopy.biomass.carbon.Rm)))
        @test abs(Rm_ample - Rm_unthrottled) <= 1e-3 * Rm_unthrottled
    end

    # Stage 2: Ra taken from the pools instead of from prescribed area indices.
    @testset "PoolBasedAutotrophicRespirationModel, FT = $FT" begin
        hydraulics = Canopy.PlantHydraulicsModel{FT}(pt, toml_dict;)
        canopy = ClimaLand.Canopy.CanopyModel{FT}(
            pt,
            (; radiation, atmos, ground),
            LAI,
            toml_dict;
            hydraulics,
            biomass,
            autotrophic_respiration = Canopy.PoolBasedAutotrophicRespirationModel{
                FT,
            }(),
        )
        Y, p, _ = initialize(canopy)
        Y.canopy.hydraulics.ϑ_l .= canopy.hydraulics.parameters.ν / 2
        Y.canopy.energy.T .= FT(290.5)
        Y.canopy.biomass.C_sugar .= FT(0.3)
        Y.canopy.biomass.C_leaf .= FT(0.2)
        Y.canopy.biomass.C_stem .= FT(5.0)
        Y.canopy.biomass.C_root .= FT(1.0)
        make_set_initial_cache(canopy)(p, Y, t0)

        M_C = FT(0.012011)
        # Ra reported to the canopy must be exactly the carbon model's own Ra,
        # in mol CO2 m^-2 s^-1. Anything else means the pools pay a different
        # respiration than the canopy reports.
        @test all(
            parent(
                @. abs(
                    p.canopy.autotrophic_respiration.Ra -
                    p.canopy.biomass.carbon.Ra / M_C,
                )
            ) .<= 10 * eps(FT),
        )
        @test all(parent(p.canopy.autotrophic_respiration.Ra) .> 0)

        # With every pool empty and no leaf area to respire, Ra must vanish.
        # The JULES model cannot do this: it respires prescribed SAI and RAI,
        # which do not go to zero over bare ground.
        Y.canopy.biomass.C_sugar .= FT(0)
        Y.canopy.biomass.C_leaf .= FT(0)
        Y.canopy.biomass.C_stem .= FT(0)
        Y.canopy.biomass.C_root .= FT(0)
        make_set_initial_cache(canopy)(p, Y, t0)
        @test all(parent(p.canopy.autotrophic_respiration.Ra) .== 0)
    end

    # Stage 3: litter must reach the soil carbon pool without being diluted or
    # amplified by the vertical discretisation. Each shape is normalised by its
    # own column integral, so the integral of the input returns the surface
    # litter flux exactly whatever the layer thickness or rooting depth.
    @testset "Litter input conserves the surface flux, FT = $FT" begin
        col = ClimaLand.Domains.Column(;
            zlim = FT.((-15, 0)),
            longlat = (FT(-92), FT(39)),
            nelements = 15,
            dz_tuple = FT.((3, 0.05)),
        )
        z = col.fields.z
        litter_depth = biomass.parameters.soil_litter_depth
        L_surface, L_root = FT(3e-8), FT(2e-8)
        for rooting_depth in (FT(0.5), FT(2.0))
            root_shape = @. Canopy.root_distribution(z, rooting_depth)
            surf_shape = @. Canopy.root_distribution(z, litter_depth)
            rn = ClimaCore.Fields.zeros(col.space.surface)
            sn = ClimaCore.Fields.zeros(col.space.surface)
            ClimaCore.Operators.column_integral_definite!(rn, root_shape)
            ClimaCore.Operators.column_integral_definite!(sn, surf_shape)
            input = @. L_surface * surf_shape / sn + L_root * root_shape / rn
            total = ClimaCore.Fields.zeros(col.space.surface)
            ClimaCore.Operators.column_integral_definite!(total, input)
            got = first(Array(parent(total)))
            @test abs(got - (L_surface + L_root)) <=
                  100 * eps(FT) * (L_surface + L_root)
            # The un-normalised shape does NOT integrate to 1 on this grid: a
            # 5 cm e-folding depth is under-resolved by a 5 cm top layer. That
            # is exactly why the normalisation is needed.
            @test first(Array(parent(sn))) < FT(0.95)
        end
    end

    # Stem turnover lengthens toward the cold. The mean annual temperature this
    # reads is seeded from air temperature at initialisation; left at zero it
    # would be 0 K and the scaling q^((T_ref-MAT)/10) would start at 2^28.
    @testset "Woody allocation scales with mean annual precipitation, FT = $FT" begin
        p_c = biomass.parameters
        h, n = p_c.map_half_woody, p_c.n_map_woody
        # Half its maximum exactly at map_half, by construction.
        @test Canopy.woody_fraction(h, h, n) ≈ FT(0.5) atol = FT(1e-6)
        # Monotone in precipitation, and bounded by (0, 1).
        @test Canopy.woody_fraction(FT(0), h, n) == FT(0)
        @test Canopy.woody_fraction(FT(0.1), h, n) <
              Canopy.woody_fraction(FT(1.0), h, n) <
              Canopy.woody_fraction(FT(5.0), h, n) <
              FT(1)
        # A wet column keeps essentially all of its woody allocation, so the
        # ramp cannot quietly suppress forests.
        @test Canopy.woody_fraction(FT(3.0), h, n) > FT(0.9)
        # half <= 0 disables the mechanism exactly, which is what lets the
        # constant-allocation control be run through the same code path.
        @test Canopy.woody_fraction(FT(0.2), FT(0), n) == FT(1)
        @test Canopy.woody_fraction(FT(0), FT(0), n) == FT(1)
        # Negative precipitation cannot arise physically, but a running sum can
        # dip below zero transiently; it must not produce a negative fraction.
        @test Canopy.woody_fraction(FT(-1), h, n) == FT(0)
    end

    @testset "Stem turnover scales with mean annual temperature, FT = $FT" begin
        p_c = biomass.parameters
        @test Canopy.tau_stem_scale(FT(283), p_c.T_ref_τ_stem, p_c.q_τ_stem) ==
              FT(1)
        @test Canopy.tau_stem_scale(FT(298), p_c.T_ref_τ_stem, p_c.q_τ_stem) ==
              FT(1)          # flat above T_ref: warm sites untouched
        @test Canopy.tau_stem_scale(FT(263), p_c.T_ref_τ_stem, p_c.q_τ_stem) ≈
              FT(4) atol = FT(1e-5)   # 20 K colder, q=2 => 2^2
        # q = 1 must recover a constant turnover exactly.
        @test Canopy.tau_stem_scale(FT(250), p_c.T_ref_τ_stem, FT(1)) == FT(1)

        hydraulics = Canopy.PlantHydraulicsModel{FT}(pt, toml_dict;)
        canopy = ClimaLand.Canopy.CanopyModel{FT}(
            pt,
            (; radiation, atmos, ground),
            LAI,
            toml_dict;
            hydraulics,
            biomass,
        )
        # The scale is capped, so even a completely unseeded state - which is
        # what a missed initial condition looks like - gives a long but finite
        # turnover rather than q^28.
        @test Canopy.tau_stem_scale(FT(0), p_c.T_ref_τ_stem, p_c.q_τ_stem) ==
              FT(Canopy.MAX_TAU_STEM_SCALE)

        # The initial condition seeds MAT from air temperature. Call it directly:
        # it runs from the simulation's initial-state pass, not from the cache.
        Y, p, _ = initialize(canopy)
        Y.canopy.hydraulics.ϑ_l .= canopy.hydraulics.parameters.ν / 2
        Y.canopy.energy.T .= FT(290.5)
        make_set_initial_cache(canopy)(p, Y, t0)
        ClimaLand.Simulations.set_canopy_component_initial_conditions!(
            Y,
            p,
            canopy.biomass,
            canopy,
        )
        @test all(parent(Y.canopy.biomass.T_annual) .> FT(200))
        @test all(isfinite.(parent(Y.canopy.biomass.T_annual)))
    end

    # Every prognostic variable needs a Jacobian block, or the implicit solver
    # refuses to build with "A does not have any entries at the following keys".
    # The tendency tests above never reach this: they call the tendency directly
    # rather than constructing a simulation, so the pools have to be checked
    # against the solver explicitly.
    @testset "Carbon pools have Jacobian entries, FT = $FT" begin
        hydraulics = Canopy.PlantHydraulicsModel{FT}(pt, toml_dict;)
        canopy = ClimaLand.Canopy.CanopyModel{FT}(
            pt,
            (; radiation, atmos, ground),
            LAI,
            toml_dict;
            hydraulics,
            biomass,
        )
        Y, _, _ = initialize(canopy)
        for pool in (:C_sugar, :C_leaf, :C_stem, :C_root, :T_annual)
            @test hasproperty(Y.canopy.biomass, pool)
        end
        @test ClimaLand.initialize_jacobian(Y) isa Any
    end

    # The carbon model must also compose with the prognostic LAI model, which
    # carries nine time-integrated prognostic variables of its own. Wrapping it
    # must preserve all of them and keep forwarding the C3/C4 competition, or
    # the Zhou LAI - and with it GPP - would silently change.
    @testset "Composition with ZhouOptimalLAIModel, FT = $FT" begin
        domain = pt
        zhou = Canopy.ZhouOptimalLAIModel{FT}(
            domain,
            toml_dict;
            SAI = FT(1.0),
            RAI = FT(1.0),
            rooting_depth = FT(0.5),
            height = FT(2.0),
        )
        wrapped = Canopy.PrognosticCarbonModel{FT}(zhou, toml_dict)

        zhou_vars = ClimaLand.prognostic_vars(zhou)
        wrapped_vars = ClimaLand.prognostic_vars(wrapped)
        @test all(v -> v in wrapped_vars, zhou_vars)
        # The wrapper adds exactly the four pools plus its own climate means;
        # naming them beats a magic count, which silently needs editing every
        # time a state variable is added.
        added = (:C_sugar, :C_leaf, :C_stem, :C_root, :T_annual, :P_annual)
        @test all(v -> v in wrapped_vars, added)
        @test length(wrapped_vars) == length(zhou_vars) + length(added)
        @test isempty(setdiff(wrapped_vars, (zhou_vars..., added...)))
        @test length(ClimaLand.prognostic_types(wrapped)) ==
              length(wrapped_vars)
        @test length(ClimaLand.prognostic_domain_names(wrapped)) ==
              length(wrapped_vars)
        # LAI itself stays prognostic under the wrapper.
        @test :LAI in wrapped_vars
        @test all(
            v -> v in ClimaLand.auxiliary_vars(wrapped),
            ClimaLand.auxiliary_vars(zhou),
        )
        @test length(ClimaLand.auxiliary_types(wrapped)) ==
              length(ClimaLand.auxiliary_vars(wrapped))
        @test wrapped.height == zhou.height
        @test wrapped.rooting_depth == zhou.rooting_depth

        # The wrapper must expose the wrapped model's diagnostics too. Anything
        # dispatching on the biomass model silently selects a default when a
        # wrapper is introduced, and here that surfaces as the output_vars
        # assertion in default_diagnostics rejecting a legitimate variable.
        zhou_diags, wrapped_diags = String[], String[]
        canopy_for_diag = ClimaLand.Canopy.CanopyModel{FT}(
            pt,
            (; radiation, atmos, ground),
            LAI,
            toml_dict;
            hydraulics = Canopy.PlantHydraulicsModel{FT}(pt, toml_dict;),
            biomass,
        )
        ClimaLand.Diagnostics.add_diagnostics!(
            zhou_diags,
            canopy_for_diag,
            zhou,
        )
        ClimaLand.Diagnostics.add_diagnostics!(
            wrapped_diags,
            canopy_for_diag,
            wrapped,
        )
        @test !isempty(zhou_diags)
        @test all(d -> d in wrapped_diags, zhou_diags)
        @test "fc3" in wrapped_diags
    end
end

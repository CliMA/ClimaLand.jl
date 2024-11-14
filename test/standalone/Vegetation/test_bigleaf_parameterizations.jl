using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaLand.Canopy

import ClimaLand
import ClimaLand.Parameters as LP


for FT in (Float32, Float64)
    @testset "Optimality Photosynthesis model Parameterizations, FT = $FT" begin
        earth_param_set = LP.LandParameters(FT)
        P = FT(101325) #Pa
        ci = FT(56.337) / P # convert to mol/mol
        oi = FT(21225.1557) / P# convert to mol/mol
        APAR = FT(1.8967 * 1e-6) # convert from μmol to mol
        θj = FT(0.85) #unitless
        ϕ = FT(0.514) # unitless
        Γstar = FT(4.332) / P # convert from Pa to mol/mol
        Kc = FT(41.03) / P# convert from Pa to mol/mol
        Ko = FT(28210) / P# convert from Pa to mol/mol
        c = FT(0.05336251)
        c_light = FT(LP.light_speed(earth_param_set))
        planck_h = FT(LP.planck_constant(earth_param_set))
        N_a = FT(LP.avogadro_constant(earth_param_set))
        Jmax, Vcmax = ClimaLand.Canopy.optimality_max_photosynthetic_rates(
            APAR,
            θj,
            ϕ,
            oi,
            ci,
            Γstar,
            Kc,
            Ko,
            c,
        )
        # We are comparing to output from another code that didnt use exactly
        # the same inputs, so we cant expect machine precision agreement
        @test abs(Jmax - FT(5.683554489145192e-7)) < 1e-10
        @test abs(Vcmax - FT(1.8572441998848603e-7)) < 1e-10
        @test typeof(Jmax) == FT
        @test typeof(Vcmax) == FT
        params = OptimalityFarquharParameters(FT)
        @test params.is_c3 ≈ 1.0
        model = OptimalityFarquharModel(params)
        @test ClimaLand.auxiliary_vars(model) == (:An, :GPP, :Rd, :Vcmax25)
        @test ClimaLand.auxiliary_types(model) == (FT, FT, FT, FT)
        @test ClimaLand.auxiliary_domain_names(model) ==
              (:surface, :surface, :surface, :surface)
        Rd = zeros(FT, 1)
        An = similar(Rd)
        Vcmax25 = similar(An)
        T = FT(280)
        β = FT(1)
        medlyn_factor = FT(10.0)
        c_co2 = FT(4.8e-4)
        R = FT(LP.gas_constant(earth_param_set))
        f_abs = FT(0.5)
        par_d = FT(1)
        λ_γ_PAR = FT(5e-7)
        energy_per_mole_photon_par = planck_h * c_light / λ_γ_PAR * N_a
        ClimaLand.Canopy.update_photosynthesis!(
            Rd,
            An,
            Vcmax25,
            model,
            T,
            f_abs,
            β,
            medlyn_factor,
            c_co2,
            R,
            energy_per_mole_photon_par,
            par_d,
        )
        @test Rd[1] != 0.0
        @test Vcmax25[1] != 0.0
        @test An[1] != 0.0
    end

    @testset "Big Leaf Parameterizations, FT = $FT" begin
        earth_param_set = LP.LandParameters(FT)
        # Test with defaults
        ARparams = AutotrophicRespirationParameters(FT)
        RTparams = BeerLambertParameters(FT)
        RT = BeerLambertModel{FT}(RTparams)
        is_c3 = FT(1) # set the photosynthesis mechanism to C3
        photosynthesisparams = FarquharParameters(FT, is_c3)
        stomatal_g_params = MedlynConductanceParameters(FT)

        LAI = FT(5.0) # m2 (leaf) m-2 (ground)
        RAI = FT(1.0)
        SAI = FT(1.0)
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        c = FT(LP.light_speed(earth_param_set))
        h = FT(LP.planck_constant(earth_param_set))
        N_a = FT(LP.avogadro_constant(earth_param_set))
        λ = FT(5e-7) # m (500 nm)
        energy_per_photon = h * c / λ

        # Drivers
        T = FT(290) # K
        P = FT(101250) #Pa
        q = FT(0.02)
        VPD = ClimaLand.vapor_pressure_deficit(T, P, q, thermo_params)#Pa
        p_l = FT(-2e6) # Pa
        ca = FT(4.11e-4) # mol/mol
        R = FT(LP.gas_constant(earth_param_set))
        θs = FT.(Array(0:0.1:(π / 2)))
        SW(θs) = cos.(θs) * FT.(500) # W/m^2
        G = compute_G(RTparams.G_Function, θs)
        K_c = extinction_coeff.(G, θs)
        α_soil_PAR = FT(0.2)
        output =
            canopy_sw_rt_beer_lambert.(
                RTparams.Ω,
                RTparams.α_PAR_leaf,
                LAI,
                K_c,
                α_soil_PAR,
            )
        FAPAR = getproperty.(output, :abs)
        FTPAR = getproperty.(output, :trans)
        FRPAR = getproperty.(output, :refl)
        @test all(@. FTPAR ≈ exp(-K_c * LAI * RTparams.Ω))
        @test all(@. FRPAR ≈ FT(1) - FAPAR - FTPAR * (1 - α_soil_PAR))

        @test all(
            @. FAPAR ≈
               (1 - RTparams.α_PAR_leaf) .* (1 - exp(-K_c * LAI * RTparams.Ω)) *
               (1 - α_soil_PAR)
        )
        To = photosynthesisparams.To
        Vcmax =
            photosynthesisparams.Vcmax25 *
            arrhenius_function(T, To, R, photosynthesisparams.ΔHVcmax)
        Kc = MM_Kc(
            photosynthesisparams.Kc25,
            photosynthesisparams.ΔHkc,
            T,
            To,
            R,
        )
        Ko = MM_Ko(
            photosynthesisparams.Ko25,
            photosynthesisparams.ΔHko,
            T,
            To,
            R,
        )
        Γstar = co2_compensation(
            photosynthesisparams.Γstar25,
            photosynthesisparams.ΔHΓstar,
            T,
            To,
            R,
        )

        @test photosynthesisparams.Kc25 *
              arrhenius_function(T, To, R, photosynthesisparams.ΔHkc) ≈ Kc
        @test photosynthesisparams.Ko25 *
              arrhenius_function(T, To, R, photosynthesisparams.ΔHko) ≈ Ko
        @test photosynthesisparams.Γstar25 *
              arrhenius_function(T, To, R, photosynthesisparams.ΔHΓstar) ≈ Γstar
        @test photosynthesisparams.Vcmax25 *
              arrhenius_function(T, To, R, photosynthesisparams.ΔHVcmax) ≈ Vcmax

        m_t = medlyn_term(stomatal_g_params.g1, T, P, q, thermo_params)

        @test m_t == 1 + stomatal_g_params.g1 / sqrt(VPD)
        ci = intercellular_co2(ca, Γstar, m_t)
        @test ci == ca * (1 - 1 / m_t)
        @test intercellular_co2(ca, FT(1), m_t) == FT(1)

        Ac = rubisco_assimilation(
            photosynthesisparams.is_c3,
            Vcmax,
            ci,
            Γstar,
            Kc,
            Ko,
            photosynthesisparams.oi,
        )
        @test Ac ==
              Vcmax * (ci - Γstar) /
              (ci + Kc * (1 + photosynthesisparams.oi / Ko))
        Jmax = max_electron_transport(
            photosynthesisparams.Vcmax25,
            photosynthesisparams.ΔHJmax,
            T,
            To,
            R,
        )
        @test Jmax ==
              photosynthesisparams.Vcmax25 *
              FT(exp(1)) *
              arrhenius_function(T, To, R, photosynthesisparams.ΔHJmax)
        APAR = FT(1)
        J =
            electron_transport.(
                APAR,
                Jmax,
                photosynthesisparams.θj,
                photosynthesisparams.ϕ,
            ) # mol m-2 s-1
        @test all(
            @.(
                photosynthesisparams.θj * J^2 -
                (photosynthesisparams.ϕ * APAR / 2 + Jmax) * J +
                photosynthesisparams.ϕ * APAR / 2 * Jmax < eps(FT)
            )
        )

        Aj = light_assimilation.(Ref(photosynthesisparams.is_c3), J, ci, Γstar)
        @test all(@.(Aj == J * (ci - Γstar) / (4 * (ci + 2 * Γstar))))
        β = moisture_stress(
            p_l,
            photosynthesisparams.sc,
            photosynthesisparams.pc,
        )
        @test β ==
              (1 + exp(photosynthesisparams.sc * photosynthesisparams.pc)) / (
            1 + exp(photosynthesisparams.sc * (p_l - photosynthesisparams.pc))
        )
        #    C4 tests
        is_c3 = 0.0
        @test rubisco_assimilation(
            is_c3,
            Vcmax,
            ci,
            Γstar,
            Kc,
            Ko,
            photosynthesisparams.oi,
        ) == Vcmax
        @test light_assimilation(is_c3, J, ci, Γstar) == J

        Rd = dark_respiration(
            photosynthesisparams.Vcmax25,
            β,
            photosynthesisparams.f,
            photosynthesisparams.ΔHRd,
            T,
            To,
            R,
        )
        @test Rd ≈
              photosynthesisparams.Vcmax25 *
              β *
              photosynthesisparams.f *
              arrhenius_function(T, To, R, photosynthesisparams.ΔHRd)
        An = net_photosynthesis.(Ac, Aj, Rd, β)
        stomatal_conductance =
            medlyn_conductance.(
                stomatal_g_params.g0,
                stomatal_g_params.Drel,
                m_t,
                An,
                ca,
            )
        @test all(
            @.(
                stomatal_conductance ≈ (
                    stomatal_g_params.g0 +
                    stomatal_g_params.Drel * m_t * (An / ca)
                )
            )
        )
        GPP = compute_GPP.(An, K_c, LAI, RTparams.Ω) # mol m-2 s-1
        @test all(
            @.(
                GPP ≈
                An * (1 - exp(-K_c * LAI * RTparams.Ω)) / (K_c * RTparams.Ω)
            )
        )

        @test all(
            @.(
                upscale_leaf_conductance(stomatal_conductance, LAI, T, R, P) ≈
                stomatal_conductance * LAI * R * T / P
            )
        )

        # Tests for Autotrophic Respiration parameterisation
        h_canopy = FT(1.0) # h planck defined above
        Nl, Nr, Ns = nitrogen_content(
            ARparams.ne,
            photosynthesisparams.Vcmax25,
            LAI,
            SAI,
            RAI,
            ARparams.ηsl,
            h_canopy,
            ARparams.σl,
            ARparams.μr,
            ARparams.μs,
        )
        Rpm = plant_respiration_maintenance(Rd, β, Nl, Nr, Ns)
        Rg = plant_respiration_growth.(ARparams.Rel, An, Rpm)

        @test Nl ==
              photosynthesisparams.Vcmax25 / ARparams.ne * ARparams.σl * LAI
        @test Nr ==
              ARparams.μr * photosynthesisparams.Vcmax25 / ARparams.ne *
              ARparams.σl *
              RAI
        @test Ns ≈
              ARparams.μs * photosynthesisparams.Vcmax25 / ARparams.ne *
              ARparams.ηsl *
              h_canopy *
              LAI *
              ClimaLand.heaviside(SAI)# == gives a very small error
        @test Rpm == Rd * (β + (Nr + Ns) / Nl)
        @test all(@.(Rg ≈ ARparams.Rel * (An - Rpm)))

    end
end

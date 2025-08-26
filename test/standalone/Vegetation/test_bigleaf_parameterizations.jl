using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaLand.Canopy

import ClimaLand
import ClimaLand.Parameters as LP


for FT in (Float32, Float64)

    @testset "Big Leaf Parameterizations, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        earth_param_set = LP.LandParameters(toml_dict)
        # Test with defaults
        ARparams = AutotrophicRespirationParameters(toml_dict)
        RTparams = BeerLambertParameters(FT)
        RT = BeerLambertModel{FT}(RTparams)
        is_c3 = FT(1) # set the photosynthesis mechanism to C3
        photosynthesisparams = FarquharParameters(FT, is_c3)
        stomatal_g_params = MedlynConductanceParameters(toml_dict)

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
        VPD =
            max(ClimaLand.vapor_pressure_deficit(T, P, q, thermo_params), FT(0))#Pa
        p_l = FT(-2e6) # Pa
        ca = FT(4.11e-4) # mol/mol
        R = FT(LP.gas_constant(earth_param_set))
        θs = FT.(Array(0:0.1:(π / 2)))
        SW(θs) = cos.(θs) * FT.(500) # W/m^2
        K = extinction_coeff.(RTparams.G_Function, cos.(θs))
        α_soil_PAR = FT(0.2)
        output =
            canopy_sw_rt_beer_lambert.(
                RTparams.G_Function,
                cos.(θs),
                RTparams.Ω,
                RTparams.α_PAR_leaf,
                LAI,
                α_soil_PAR,
            )
        FAPAR = getproperty.(output, :abs)
        FTPAR = getproperty.(output, :trans)
        FRPAR = getproperty.(output, :refl)
        @test all(@. FTPAR ≈ exp(-K * LAI * RTparams.Ω))
        @test all(@. FRPAR ≈ FT(1) - FAPAR - FTPAR * (1 - α_soil_PAR))

        @test all(
            @. FAPAR ≈
               (1 - RTparams.α_PAR_leaf) .* (1 - exp(-K * LAI * RTparams.Ω)) *
               (1 - α_soil_PAR)
        )
        To = photosynthesisparams.To
        Vcmax =
            photosynthesisparams.Vcmax25 *
            arrhenius_function(T, To, R, photosynthesisparams.ΔHVcmax)
        Kc = ClimaLand.Canopy.MM_Kc(
            photosynthesisparams.Kc25,
            photosynthesisparams.ΔHkc,
            T,
            To,
            R,
        )
        Ko = ClimaLand.Canopy.MM_Ko(
            photosynthesisparams.Ko25,
            photosynthesisparams.ΔHko,
            T,
            To,
            R,
        )
        Γstar = ClimaLand.Canopy.co2_compensation_farquhar(
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

        @test m_t == 1 + stomatal_g_params.g1 / sqrt(sqrt(eps(FT)) + VPD)
        ci = ClimaLand.Canopy.intercellular_co2_farquhar(ca, Γstar, m_t)
        @test ci == ca * (1 - 1 / m_t)
        @test ClimaLand.Canopy.intercellular_co2_farquhar(ca, FT(1), m_t) ==
              FT(1)
        Kmm = @. Kc * (1 + photosynthesisparams.oi / Ko)
        Ac = ClimaLand.Canopy.rubisco_assimilation_farquhar(
            photosynthesisparams.is_c3,
            Vcmax,
            ci,
            Γstar,
            Kmm,
        )
        @test Ac ==
              Vcmax * (ci - Γstar) /
              (ci + Kc * (1 + photosynthesisparams.oi / Ko))
        Jmax = ClimaLand.Canopy.max_electron_transport_farquhar(
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
            ClimaLand.Canopy.electron_transport_farquhar.(
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

        Aj =
            ClimaLand.Canopy.light_assimilation_farquhar.(
                Ref(photosynthesisparams.is_c3),
                J,
                ci,
                Γstar,
            )
        @test all(@.(Aj == J * (ci - Γstar) / (4 * (ci + 2 * Γstar))))
        β = compute_tuzet_moisture_stress(
            p_l,
            photosynthesisparams.pc,
            photosynthesisparams.sc,
        )
        @test β ==
              (1 + exp(photosynthesisparams.sc * photosynthesisparams.pc)) / (
            1 + exp(photosynthesisparams.sc * (p_l - photosynthesisparams.pc))
        )
        #    C4 tests

        is_c3 = 0.0
        (;
            Q10,
            s1,
            s2,
            s3,
            s4,
            s5,
            s6,
            Vcmax25,
            ΔHVcmax,
            E,
            fC4,
            fC3,
            To,
            ΔHRd,
        ) = photosynthesisparams
        @test ClimaLand.Canopy.compute_Vcmax_farquhar(
            is_c3,
            Vcmax25,
            T,
            R,
            To,
            ΔHVcmax,
            Q10,
            s1,
            s2,
            s3,
            s4,
        ) ==
              Vcmax25 * Q10^((T - To) / 10) / (1 + exp(s1 * (T - s2))) /
              (1 + exp(s3 * (s4 - T)))
        Kmm = @. Kc * (1 + photosynthesisparams.oi / Ko)
        @test ClimaLand.Canopy.rubisco_assimilation_farquhar(
            is_c3,
            Vcmax,
            ci,
            Γstar,
            Kmm,
        ) == Vcmax
        @test ClimaLand.Canopy.light_assimilation_farquhar(
            is_c3,
            J,
            ci,
            Γstar,
            APAR,
            E,
        ) == APAR * E

        Rd = ClimaLand.Canopy.dark_respiration_farquhar(
            is_c3,
            Vcmax25,
            β,
            T,
            R,
            To,
            fC3,
            ΔHRd,
            Q10,
            s5,
            s6,
            fC4,
        )
        @test Rd ≈
              photosynthesisparams.Vcmax25 * β * fC4 * Q10^((T - To) / 10) /
              (1 + exp(s5 * (T - s6)))
        A = ClimaLand.Canopy.gross_photosynthesis.(Ac, Aj)
        An = ClimaLand.Canopy.net_photosynthesis.(A .* β, Rd)
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

        @test all(
            @.(
                ClimaLand.Canopy.conductance_molar_flux_to_m_per_s(
                    stomatal_conductance,
                    T,
                    R,
                    P,
                ) ≈ stomatal_conductance * R * T / P
            )
        )

        # Tests for Autotrophic Respiration parameterisation
        h_canopy = FT(1.0) # h planck defined above
        Nl, Nr, Ns = ClimaLand.Canopy.nitrogen_content(
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
        Rpm = ClimaLand.Canopy.plant_respiration_maintenance(Rd, β, Nl, Nr, Ns)
        Rg = ClimaLand.Canopy.plant_respiration_growth.(ARparams.Rel, An, Rpm)

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

using Test
using ClimaLand


for FT in (Float32, Float64)
    @testset "Roundoff bug, FT = $FT" begin
        Ω = FT(0.69)
        χl = FT(0.1)
        G_Function = ClimaLand.Canopy.CLMGFunction(χl)
        α_PAR_leaf = FT(0.1)
        λ_γ_PAR = FT(5e-7)
        λ_γ_NIR = FT(1.65e-6)
        τ_PAR_leaf = FT(0.05)
        α_NIR_leaf = FT(0.45)
        τ_NIR_leaf = FT(0.25)
        ϵ_canopy = FT(0.97)
        ts_params = ClimaLand.Canopy.TwoStreamParameters(
            α_PAR_leaf,
            τ_PAR_leaf,
            α_NIR_leaf,
            τ_NIR_leaf,
            ϵ_canopy,
            Ω,
            λ_γ_PAR,
            λ_γ_NIR,
            UInt64(20),
            G_Function,
            FT(0.05),
        )

        BL_params = ClimaLand.Canopy.BeerLambertParameters(
            α_PAR_leaf,
            α_NIR_leaf,
            ϵ_canopy,
            Ω,
            λ_γ_PAR,
            λ_γ_NIR,
            G_Function,
            FT(0.05),
        )

        RTmodels = [
            ClimaLand.Canopy.TwoStreamModel{FT, typeof(ts_params)}(ts_params),
            ClimaLand.Canopy.BeerLambertModel{FT, typeof(BL_params)}(BL_params),
        ]

        PAR = FT(1.0)

        LAI = eps(FT)
        θs = FT(π / 2)
        K = ClimaLand.Canopy.extinction_coeff(G_Function, θs)
        α_soil_PAR = FT(0.2)
        frac_diff = FT(0.0)

        RT = RTmodels[1]
        canopy_par_ts = @. ClimaLand.Canopy.plant_absorbed_pfd(
            RT,
            PAR,
            RT.parameters.α_PAR_leaf,
            RT.parameters.τ_PAR_leaf,
            LAI,
            K,
            θs,
            α_soil_PAR,
            frac_diff,
        )[1]
        # INC - OUT = CANOPY_ABS + (1-α_soil)*CANOPY_TRANS
        # OUT = REFL
        @test canopy_par_ts.abs < eps(FT)
        @test canopy_par_ts.refl ≈ FT(α_soil_PAR)
        @test canopy_par_ts.trans ≈ FT(1)

        RT = RTmodels[2]
        canopy_par_bl = @. ClimaLand.Canopy.plant_absorbed_pfd(
            RT,
            PAR,
            RT.parameters.α_PAR_leaf,
            LAI,
            K,
            α_soil_PAR,
        )[1]
        # INC - OUT = CANOPY_ABS + (1-α_soil)*CANOPY_TRANS
        # OUT = REFL
        @test canopy_par_bl.abs < eps(FT)
        @test canopy_par_bl.refl ≈ FT(α_soil_PAR)
        @test canopy_par_bl.trans ≈ FT(1)

    end
end

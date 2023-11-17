using Test
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: HybridBox, Column, obtain_surface_domain
using ClimaLSM.Soil
using ClimaLSM.Pond

for FT in (Float32, Float64)
    @testset "Pond soil Multi Column LSM integration test, FT = $FT" begin
        function precipitation(t)
            if t < 20
                precip = -1e-8
            else
                precip = t < 100 ? -5e-5 : 0
            end
            return precip
        end
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten(; α = vg_α, n = vg_n)
        θ_r = FT(0)
        zmax = FT(0)
        zmin = FT(-1)
        nelems = 20
        lsm_domain = HybridBox(;
            xlim = (FT(0), FT(1)),
            ylim = (FT(0), FT(1)),
            zlim = (zmin, zmax),
            nelements = (1, 1, nelems),
            npolynomial = 1,
            periodic = (true, true),
        )
        soil_ps =
            Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)
        soil_args = (; domain = lsm_domain, parameters = soil_ps)
        surface_water_args = (; domain = obtain_surface_domain(lsm_domain))

        land_args = (; precip = precipitation)

        land = LandHydrology{FT}(;
            land_args = land_args,
            soil_model_type = Soil.RichardsModel{FT},
            soil_args = soil_args,
            surface_water_model_type = Pond.PondModel{FT},
            surface_water_args = surface_water_args,
        )
        Y, p, coords = initialize(land)
        # test that the dss buffers are correctly added
        @test propertynames(p) == (
            :soil_infiltration,
            :soil,
            :surface_water,
            :dss_buffer_3d,
            :dss_buffer_2d,
        )
        @test typeof(p.dss_buffer_3d) == typeof(
            ClimaCore.Spaces.create_dss_buffer(
                ClimaCore.Fields.zeros(land.soil.domain.space.subsurface),
            ),
        )
        @test typeof(p.dss_buffer_2d) == typeof(
            ClimaCore.Spaces.create_dss_buffer(
                ClimaCore.Fields.zeros(land.surface_water.domain.space.surface),
            ),
        )
        function init_soil!(Ysoil, coords, params)
            function hydrostatic_profile(
                z::FT,
                params::RichardsParameters{FT},
            ) where {FT}
                (; ν, hydrology_cm, θ_r) = params
                (; α, n) = hydrology_cm
                m = 1 - 1 / n
                #unsaturated zone only, assumes water table starts at z_∇
                z_∇ = FT(-1)# matches zmin
                S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
                ϑ_l = S * (ν - θ_r) + θ_r
                return FT(ϑ_l)
            end
            Ysoil.soil.ϑ_l .= hydrostatic_profile.(coords.z, Ref(params))
        end
        init_soil!(Y, coords.subsurface, land.soil.parameters)
        # initialize the pond height to zero
        Y.surface_water.η .= FT(0)
        dY = similar(Y)

        exp_tendency! = make_exp_tendency(land)
        t = Float64(0)
        dt = Float64(1)
        # set aux state values to those corresponding with Y(t=0)
        set_initial_aux_state! = make_set_initial_aux_state(land)
        set_initial_aux_state!(p, Y, t)
        # integrate
        for step in 1:10
            exp_tendency!(dY, Y, p, t)
            t = t + dt
            @. Y.surface_water.η = dY.surface_water.η * dt + Y.surface_water.η
        end
        # it running is the test

        # Infiltration at point test

        η = FT.([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
        i_c = FT.([2.0, 2.0, -2.0, -2.0, 2.0, 2.0])
        P = FT.([3.0, -1.0, 3.0, -3.0, 1.0, 3.0])
        iap = FT.([2.0, -1.0, -2.0, -3.0, 2.0, 2.0])
        @test infiltration_at_point.(η, i_c, P) ≈ iap

        land_prognostic_vars = prognostic_vars(land)
        land_prognostic_types = prognostic_types(land)
        land_auxiliary_vars = auxiliary_vars(land)
        land_auxiliary_types = auxiliary_types(land)
        # prognostic vars
        @test land_prognostic_vars.soil == prognostic_vars(land.soil)
        @test land_prognostic_vars.surface_water ==
              prognostic_vars(land.surface_water)
        # auxiliary vars
        @test land_auxiliary_vars.soil == auxiliary_vars(land.soil)
        @test land_auxiliary_vars.surface_water ==
              auxiliary_vars(land.surface_water)
        @test land_auxiliary_vars.interactions ==
              ClimaLSM.interaction_vars(land)
        # prognostic types
        @test land_prognostic_types.soil == prognostic_types(land.soil)
        @test land_prognostic_types.surface_water ==
              prognostic_types(land.surface_water)
        # auxiliary types
        @test land_auxiliary_types.soil == auxiliary_types(land.soil)
        @test land_auxiliary_types.surface_water ==
              auxiliary_types(land.surface_water)
        @test land_auxiliary_types.interactions ==
              ClimaLSM.interaction_types(land)
    end


    @testset "Pond soil Single Column LSM integration test, FT = $FT" begin
        function precipitation(t)
            if t < 20
                precip = -1e-8
            else
                precip = t < 100 ? -5e-5 : 0.0
            end
            return precip
        end
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten(; α = vg_α, n = vg_n)
        θ_r = FT(0)
        zmax = FT(0)
        zmin = FT(-1)
        nelems = 20
        lsm_domain = Column(; zlim = (zmin, zmax), nelements = nelems)

        soil_ps =
            Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)
        soil_args = (; domain = lsm_domain, parameters = soil_ps)
        surface_water_args = (; domain = obtain_surface_domain(lsm_domain))

        land_args = (; precip = precipitation)

        land = LandHydrology{FT}(;
            land_args = land_args,
            soil_model_type = Soil.RichardsModel{FT},
            soil_args = soil_args,
            surface_water_model_type = Pond.PondModel{FT},
            surface_water_args = surface_water_args,
        )
        Y, p, coords = initialize(land)
        @test propertynames(p) == (:soil_infiltration, :soil, :surface_water)
        function init_soil!(Ysoil, coords, params)
            function hydrostatic_profile(
                z::FT,
                params::RichardsParameters{FT},
            ) where {FT}
                (; ν, hydrology_cm, θ_r) = params
                (; α, n) = hydrology_cm
                m = 1 - 1 / n
                #unsaturated zone only, assumes water table starts at z_∇
                z_∇ = FT(-1)# matches zmin
                S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
                ϑ_l = S * (ν - θ_r) + θ_r
                return FT(ϑ_l)
            end
            Ysoil.soil.ϑ_l .= hydrostatic_profile.(coords.z, Ref(params))
        end
        init_soil!(Y, coords.subsurface, land.soil.parameters)
        # initialize the pond height to zero
        Y.surface_water.η .= FT(0)
        dY = similar(Y)
        exp_tendency! = make_exp_tendency(land)
        t = Float64(0)
        dt = Float64(1)
        for step in 1:10
            exp_tendency!(dY, Y, p, t)
            t = t + dt
            @. Y.surface_water.η = dY.surface_water.η * dt + Y.surface_water.η
        end

        # it running is the test
    end
end

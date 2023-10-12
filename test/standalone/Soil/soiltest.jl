using Test
using Statistics
using ClimaCore
import CLIMAParameters as CP
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil

import ClimaLSM
import ClimaLSM.Parameters as LSMP

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))


for FT in (Float32, Float64)
    @testset "Richards equation with flux BCs, FT = $FT" begin
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten(; α = vg_α, n = vg_n)
        θ_r = FT(0)
        zmax = FT(0)
        zmin = FT(-10)
        nelems = 50

        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        top_flux_bc = FluxBC((p, t) -> 0.0)
        bot_flux_bc = FluxBC((p, t) -> 0.0)
        sources = ()
        boundary_fluxes =
            (; top = (water = top_flux_bc,), bottom = (water = bot_flux_bc,))
        params =
            Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)

        soil = Soil.RichardsModel{FT}(;
            parameters = params,
            domain = soil_domain,
            boundary_conditions = boundary_fluxes,
            sources = sources,
        )

        Y, p, coords = initialize(soil)
        # Here we do not add a dss buffer: check that `initialize`
        # has not added a buffer to `p` and that it only contains the
        # `soil` variables.
        @test propertynames(p) == (:soil,)
        # specify ICs
        function init_soil!(Ysoil, z, params)
            function hydrostatic_profile(
                z::FT,
                params::RichardsParameters{FT},
            ) where {FT}
                (; ν, hydrology_cm, θ_r) = params
                (; α, m, n) = hydrology_cm
                #unsaturated zone only, assumes water table starts at z_∇
                z_∇ = FT(-10)# matches zmin
                S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
                ϑ_l = S * (ν - θ_r) + θ_r
                return FT(ϑ_l)
            end
            Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
        end

        init_soil!(Y, coords.subsurface.z, soil.parameters)

        t0 = FT(0)
        set_initial_aux_state! = make_set_initial_aux_state(soil)
        set_initial_aux_state!(p, Y, t0)
        dY = similar(Y)

        imp_tendency! = make_imp_tendency(soil)
        exp_tendency! = make_exp_tendency(soil)
        imp_tendency!(dY, Y, p, t0)
        exp_tendency!(dY, Y, p, t0)
        ClimaLSM.dss!(dY, p, t0)

        @test mean(parent(dY)) < eps(FT)
        # should be hydrostatic equilibrium at every layer, at each step:
        @test mean(parent(p.soil.ψ .+ coords.subsurface.z)[:] .+ FT(10)) <
              eps(FT)
    end

    @testset "Soil Energy and Water tendency unit tests, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        vg_m = FT(1) - FT(1) / vg_n
        hcm = vanGenuchten(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)
        κ_minerals = FT(2.5)
        κ_om = FT(0.25)
        κ_quartz = FT(8.0)
        κ_air = FT(0.025)
        κ_ice = FT(2.21)
        κ_liq = FT(0.57)
        ρp = FT(2.66 / 1e3 * 1e6)
        ρc_ds = FT(2e6 * (1.0 - ν))
        κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
        κ_dry = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
        κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
        κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)
        zmax = FT(0)
        zmin = FT(-1)
        nelems = 200
        Δz = abs(zmax - zmin) / nelems
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        top_flux_bc = FluxBC((p, t) -> 0.0)
        bot_flux_bc = FluxBC((p, t) -> 0.0)
        sources = ()
        # Use the same BCs for RRE and heat
        boundary_fluxes = (;
            top = (water = top_flux_bc, heat = top_flux_bc),
            bottom = (water = bot_flux_bc, heat = bot_flux_bc),
        )


        ###
        hyd_off_en_on = Soil.EnergyHydrologyParameters{FT}(;
            κ_dry = κ_dry,
            κ_sat_frozen = κ_sat_frozen,
            κ_sat_unfrozen = κ_sat_unfrozen,
            ρc_ds = ρc_ds,
            ν = ν,
            ν_ss_om = ν_ss_om,
            ν_ss_quartz = ν_ss_quartz,
            ν_ss_gravel = ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat = FT(0),
            S_s = S_s,
            θ_r = θ_r,
            earth_param_set = earth_param_set,
        )
        soil_heat_on = Soil.EnergyHydrology{FT}(;
            parameters = hyd_off_en_on,
            domain = soil_domain,
            boundary_conditions = boundary_fluxes,
            sources = sources,
        )

        Y, p, coords = initialize(soil_heat_on)

        # specify ICs
        function init_soil_heat!(Ysoil, z, params)
            FT = eltype(Ysoil.soil.ϑ_l)
            ν = params.ν
            Ysoil.soil.ϑ_l .= ν / 2.0
            Ysoil.soil.θ_i .= 0.0
            T = FT.(280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0)
            ρc_s = @. Soil.volumetric_heat_capacity(
                Ysoil.soil.ϑ_l,
                Ysoil.soil.θ_i,
                params,
            )
            Ysoil.soil.ρe_int = @. Soil.volumetric_internal_energy.(
                Ysoil.soil.θ_i,
                ρc_s,
                T,
                params,
            )
        end

        t0 = FT(0)
        init_soil_heat!(Y, coords.subsurface.z, soil_heat_on.parameters)
        set_initial_aux_state! = make_set_initial_aux_state(soil_heat_on)
        set_initial_aux_state!(p, Y, t0)
        exp_tendency! = make_exp_tendency(soil_heat_on)
        dY = similar(Y)
        exp_tendency!(dY, Y, p, t0)
        ClimaLSM.dss!(dY, p, t0)

        F_face = FT(0)
        κ = parent(p.soil.κ)
        F_below = -0.5 * (κ[end] + κ[end - 1]) * (-Δz + 0.5)
        dY_top = -(F_face - F_below) / Δz
        F_top = -0.5 * (κ[2] + κ[1]) * ((-1.0 + Δz) + 0.5)
        dY_bot = -(F_top - F_face) / Δz
        expected = parent(p.soil.κ)
        expected[1] = dY_bot
        expected[end] = dY_top
        @test mean(abs.(expected .- parent(dY.soil.ρe_int))) /
              median(parent(Y.soil.ρe_int)) < eps(FT)
        @test maximum(abs.(parent(dY.soil.ϑ_l))) == FT(0)
        @test maximum(abs.(parent(dY.soil.θ_i))) == FT(0)

        ###
        hyd_on_en_off = Soil.EnergyHydrologyParameters{FT}(
            κ_dry = FT(0),
            κ_sat_frozen = FT(0),
            κ_sat_unfrozen = FT(0),
            ρc_ds = ρc_ds,
            ν = ν,
            ν_ss_om = ν_ss_om,
            ν_ss_quartz = ν_ss_quartz,
            ν_ss_gravel = ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat = K_sat,
            S_s = S_s,
            θ_r = θ_r,
            earth_param_set = earth_param_set,
        )
        soil_water_on = Soil.EnergyHydrology{FT}(;
            parameters = hyd_on_en_off,
            domain = soil_domain,
            boundary_conditions = boundary_fluxes,
            sources = sources,
        )

        Y, p, coords = initialize(soil_water_on)

        # specify ICs
        function init_soil_water!(Ysoil, z, params)
            FT = eltype(Ysoil.soil.ϑ_l)
            ν = params.ν
            Ysoil.soil.ϑ_l .= ν / 2.0 .+ ν / 4.0 .* (z .+ 0.5) .^ 2.0
            Ysoil.soil.θ_i .= 0.0
            ρc_s = @. Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Ysoil.soil.θ_i,
                params,
            )
            Ysoil.soil.ρe_int = @. volumetric_internal_energy.(
                Ysoil.soil.θ_i,
                ρc_s,
                FT(288),
                params,
            )
        end

        t0 = FT(0)
        init_soil_water!(Y, coords.subsurface.z, soil_water_on.parameters)
        set_initial_aux_state! = make_set_initial_aux_state(soil_water_on)
        set_initial_aux_state!(p, Y, t0)
        exp_tendency! = make_exp_tendency(soil_water_on)
        dY = similar(Y)
        exp_tendency!(dY, Y, p, t0)
        ClimaLSM.dss!(dY, p, t0)

        function dKdθ(θ::FT)::FT where {FT}
            S = (θ - θ_r) / (ν - θ_r)
            if S < 1
                f::FT = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
                f1::FT = f^2.0 / 2.0 / S^0.5
                f2::FT =
                    2 * S^(1 / vg_m - 1 / 2) * f /
                    (1 - S^(1 / vg_m))^(1.0 - vg_m)
                return (f1 + f2) * K_sat / (ν - θ_r)
            else
                return FT(0)
            end
        end

        function dψdθ(θ::FT)::FT where {FT}
            S = (θ - θ_r) / (ν - θ_r)
            if S < 1.0
                return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
                       (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
                       S^(-1 / vg_m - 1)
            else
                return 1.0 / S_s
            end
        end


        function d2ψdθ2(θ::FT)::FT where {FT}
            S = (θ - θ_r) / (ν - θ_r)
            if S < 1.0
                return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r)^2.0 * (
                    S^(-2.0 / vg_m - 2.0) *
                    (-1 / vg_m) *
                    (1 / vg_n - 1) *
                    (S^(-1 / vg_m) - 1)^(1 / vg_n - 2) +
                    (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
                    (-1 / vg_m - 1) *
                    S^(-1 / vg_m - 2)
                )
            else
                return FT(0)
            end
        end

        function K(θ::FT)::FT where {FT}
            S = (θ - θ_r) / (ν - θ_r)
            if S < 1
                f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
                return K_sat * f^2.0 * sqrt(S)
            else
                return K_sat
            end

        end

        function dθdz(z::FT)::FT where {FT}
            return ν / 2.0 * (z + 0.5)
        end


        function d2θdz2(z::FT)::FT where {FT}
            return ν / 2.0
        end

        θ = parent(Y.soil.ϑ_l)# on the center
        θ_face = 0.5 * (θ[2:end] + θ[1:(end - 1)])
        Z = parent(coords.subsurface.z)
        Z_face = 0.5 * (Z[2:end] + Z[1:(end - 1)])
        K_face = 0.5 .* (K.(θ[2:end]) .+ K.(θ[1:(end - 1)]))
        flux_interior = (@. -K_face * (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
        flux = vcat(
            [top_flux_bc.bc(p, FT(0.0))],
            flux_interior,
            [bot_flux_bc.bc(p, FT(0.0))],
        )

        expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz
        @test mean(abs.(expected .- parent(dY.soil.ϑ_l))) / ν < 10^2 * eps(FT)
        @test maximum(abs.(parent(dY.soil.θ_i))) == FT(0)

        ρe_int_l = parent(
            Soil.volumetric_internal_energy_liq.(
                p.soil.T,
                Ref(soil_water_on.parameters),
            ),
        )
        ρe_int_l_face = 0.5 * (ρe_int_l[2:end] + ρe_int_l[1:(end - 1)])
        flux_interior =
            (@. -K_face * ρe_int_l_face * (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
        flux = vcat(
            [top_flux_bc.bc(p, FT(0.0))],
            flux_interior,
            [bot_flux_bc.bc(p, FT(0.0))],
        )
        expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz

        @test mean(abs.(expected .- parent(dY.soil.ρe_int))) /
              median(parent(Y.soil.ρe_int)) < 10^2 * eps(FT)

        ###

        hyd_off_en_off = Soil.EnergyHydrologyParameters{FT}(;
            κ_dry = FT(0),
            κ_sat_frozen = FT(0),
            κ_sat_unfrozen = FT(0),
            ρc_ds = ρc_ds,
            ν = ν,
            ν_ss_om = ν_ss_om,
            ν_ss_quartz = ν_ss_quartz,
            ν_ss_gravel = ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat = FT(0),
            S_s = S_s,
            θ_r = θ_r,
            earth_param_set = earth_param_set,
        )

        soil_both_off = Soil.EnergyHydrology{FT}(;
            parameters = hyd_off_en_off,
            domain = soil_domain,
            boundary_conditions = boundary_fluxes,
            sources = sources,
        )

        Y, p, coords = initialize(soil_both_off)

        # specify ICs
        function init_soil_off!(Ysoil, z, params)
            FT = eltype(Ysoil.soil.ϑ_l)
            ν = params.ν
            Ysoil.soil.ϑ_l .= ν / 2.0
            Ysoil.soil.θ_i .= 0.0
            T = FT.(280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0)
            ρc_s = @. Soil.volumetric_heat_capacity(
                Ysoil.soil.ϑ_l,
                Ysoil.soil.θ_i,
                params,
            )
            Ysoil.soil.ρe_int = @. Soil.volumetric_internal_energy.(
                Ysoil.soil.θ_i,
                ρc_s,
                T,
                params,
            )
        end

        t0 = FT(0)
        init_soil_off!(Y, coords.subsurface.z, soil_both_off.parameters)
        set_initial_aux_state! = make_set_initial_aux_state(soil_both_off)
        set_initial_aux_state!(p, Y, t0)
        exp_tendency! = make_exp_tendency(soil_both_off)
        dY = similar(Y)
        exp_tendency!(dY, Y, p, t0)
        ClimaLSM.dss!(dY, p, t0)

        @test maximum(abs.(parent(dY.soil.ρe_int))) == FT(0)
        @test maximum(abs.(parent(dY.soil.ϑ_l))) == FT(0)
        @test maximum(abs.(parent(dY.soil.θ_i))) == FT(0)


        ### Test with both energy and hydrology on
        hyd_on_en_on = Soil.EnergyHydrologyParameters{FT}(;
            κ_dry = κ_dry,
            κ_sat_frozen = κ_sat_frozen,
            κ_sat_unfrozen = κ_sat_unfrozen,
            ρc_ds = ρc_ds,
            ν = ν,
            ν_ss_om = ν_ss_om,
            ν_ss_quartz = ν_ss_quartz,
            ν_ss_gravel = ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat = K_sat,
            S_s = S_s,
            θ_r = θ_r,
            earth_param_set = earth_param_set,
        )

        soil_both_on = Soil.EnergyHydrology{FT}(;
            parameters = hyd_on_en_on,
            domain = soil_domain,
            boundary_conditions = boundary_fluxes,
            sources = sources,
        )

        Y, p, coords = initialize(soil_both_on)

        # specify ICs
        function init_soil_on!(Ysoil, z, params)
            FT = eltype(Ysoil.soil.ϑ_l)
            ν = params.ν
            Ysoil.soil.ϑ_l .= ν / 2.0 .+ ν / 4.0 .* (z .+ 0.5) .^ 2.0
            Ysoil.soil.θ_i .= 0.0
            T = FT.(280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0)
            ρc_s = @. Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Ysoil.soil.θ_i,
                params,
            )
            Ysoil.soil.ρe_int =
                @. volumetric_internal_energy.(Ysoil.soil.θ_i, ρc_s, T, params)
        end

        t0 = FT(0)
        init_soil_on!(Y, coords.subsurface.z, soil_both_on.parameters)
        set_initial_aux_state! = make_set_initial_aux_state(soil_both_on)
        set_initial_aux_state!(p, Y, t0)

        exp_tendency! = make_exp_tendency(soil_both_on)
        dY = similar(Y)
        exp_tendency!(dY, Y, p, t0)
        ClimaLSM.dss!(dY, p, t0)

        @test maximum(abs.(parent(dY.soil.θ_i))) == FT(0)

        θ = parent(Y.soil.ϑ_l)# on the center
        θ_face = 0.5 * (θ[2:end] + θ[1:(end - 1)])
        Z = parent(coords.subsurface.z)
        Z_face = 0.5 * (Z[2:end] + Z[1:(end - 1)])
        K_face = 0.5 .* (K.(θ[2:end]) .+ K.(θ[1:(end - 1)]))
        vf = parent(
            Soil.viscosity_factor.(
                p.soil.T,
                hyd_on_en_on.γ,
                hyd_on_en_on.γT_ref,
            ),
        )
        vf_face = 0.5 .* (vf[2:end] .+ vf[1:(end - 1)])

        flux_interior =
            (@. -K_face * vf_face * (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
        flux = vcat(
            [top_flux_bc.bc(p, FT(0.0))],
            flux_interior,
            [bot_flux_bc.bc(p, FT(0.0))],
        )

        expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz
        @test mean(abs.(expected .- parent(dY.soil.ϑ_l))) / ν < 10^2 * eps(FT)

        ρe_int_l = parent(
            Soil.volumetric_internal_energy_liq.(
                p.soil.T,
                Ref(soil_both_on.parameters),
            ),
        )
        ρe_int_l_face = 0.5 * (ρe_int_l[2:end] + ρe_int_l[1:(end - 1)])
        κ = parent(p.soil.κ)
        κ_face = 0.5 * (κ[2:end] + κ[1:(end - 1)])


        flux_interior = (@. -K_face *
            vf_face *
            ρe_int_l_face *
            (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
        flux_one = vcat(
            [top_flux_bc.bc(p, FT(0.0))],
            flux_interior,
            [bot_flux_bc.bc(p, FT(0.0))],
        )
        flux_interior = (@. -κ_face * (Z_face + 0.5))
        flux_two = vcat(
            [top_flux_bc.bc(p, FT(0.0))],
            flux_interior,
            [bot_flux_bc.bc(p, FT(0.0))],
        )
        flux = flux_one .+ flux_two
        expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz

        @test mean(abs.(expected .- parent(dY.soil.ρe_int))) /
              median(parent(Y.soil.ρe_int)) < 10^2 * eps(FT)
    end

    @testset "Phase change source term, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        ν = FT(0.495)
        K_sat = FT(0.0) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)
        κ_minerals = FT(2.5)
        κ_om = FT(0.25)
        κ_quartz = FT(8.0)
        κ_air = FT(0.025)
        κ_ice = FT(2.21)
        κ_liq = FT(0.57)
        ρp = FT(2.66 / 1e3 * 1e6)
        ρc_ds = FT(2e6 * (1.0 - ν))
        κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
        κ_dry = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
        κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
        κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)
        zmax = FT(0)
        zmin = FT(-1)
        nelems = 200
        Δz = FT(0.005)
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        top_flux_bc = FluxBC((p, t) -> 0.0)
        bot_flux_bc = FluxBC((p, t) -> 0.0)
        boundary_fluxes = (;
            top = (water = top_flux_bc, heat = top_flux_bc),
            bottom = (water = bot_flux_bc, heat = bot_flux_bc),
        )

        sources = (PhaseChange{FT}(Δz),)

        ###
        hyd_off_en_on = Soil.EnergyHydrologyParameters{FT}(;
            κ_dry = κ_dry,
            κ_sat_frozen = κ_sat_frozen,
            κ_sat_unfrozen = κ_sat_unfrozen,
            ρc_ds = ρc_ds,
            ν = ν,
            ν_ss_om = ν_ss_om,
            ν_ss_quartz = ν_ss_quartz,
            ν_ss_gravel = ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat = K_sat,
            S_s = S_s,
            θ_r = θ_r,
            earth_param_set = earth_param_set,
        )
        soil_heat_on = Soil.EnergyHydrology{FT}(;
            parameters = hyd_off_en_on,
            domain = soil_domain,
            boundary_conditions = boundary_fluxes,
            sources = sources,
        )

        Y, p, coords = initialize(soil_heat_on)

        # specify ICs
        function init_soil_heat!(Ysoil, z, params)
            FT = eltype(Ysoil.soil.ϑ_l)
            ν = params.ν
            Ysoil.soil.ϑ_l .= ν / 2.0
            Ysoil.soil.θ_i .= 0.0
            T = FT(270)
            ρc_s = @. Soil.volumetric_heat_capacity(
                Ysoil.soil.ϑ_l,
                Ysoil.soil.θ_i,
                params,
            )
            Ysoil.soil.ρe_int = @. Soil.volumetric_internal_energy.(
                Ysoil.soil.θ_i,
                ρc_s,
                T,
                params,
            )
        end

        init_soil_heat!(Y, coords.subsurface.z, soil_heat_on.parameters)
        set_initial_aux_state! = make_set_initial_aux_state(soil_heat_on)
        set_initial_aux_state!(p, Y, FT(0.0))
        dY = similar(Y)
        dY .= FT(0.0)
        source!(dY, sources[1], Y, p, soil_heat_on)
        _ρ_l = FT(LSMP.ρ_cloud_liq(soil_heat_on.parameters.earth_param_set))
        _ρ_i = FT(LSMP.ρ_cloud_ice(soil_heat_on.parameters.earth_param_set))
        @test parent(dY.soil.ϑ_l) ≈ -(_ρ_i / _ρ_l) .* parent(dY.soil.θ_i)
    end
end

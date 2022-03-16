using Test
using UnPack
using Statistics
using DifferentialEquations
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil

FT = Float64

@testset "Hydraulic functions" begin
    @test volumetric_liquid_fraction(0.1, 0.2) == 0.1
    @test volumetric_liquid_fraction(0.2, 0.1) == 0.1
end


@testset "Richards equation" begin

    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    ν = FT(0.495)
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    zmax = FT(0)
    zmin = FT(-10)
    nelems = 50

    soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelems)
    top_flux_bc = FT(0.0)
    bot_flux_bc = FT(0.0)
    sources = ()
    boundary_fluxes = FluxBC{FT}(top_flux_bc, bot_flux_bc)
    params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r)

    soil = Soil.RichardsModel{FT}(;
        param_set = params,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil)

    # specify ICs
    function init_soil!(Ysoil, z, params)
        function hydrostatic_profile(
            z::FT,
            params::RichardsParameters{FT},
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r = params
            #unsaturated zone only, assumes water table starts at z_∇
            z_∇ = FT(-10)# matches zmin
            S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
            ϑ_l = S * (ν - θ_r) + θ_r
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
    end

    init_soil!(Y, coords.z, soil.param_set)

    soil_ode! = make_ode_function(soil)

    t0 = FT(0)
    tf = FT(60)
    dt = FT(1)
    cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values)
    prob = ODEProblem(soil_ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(); dt = dt, callback = cb)

    @test sum(parent(sol.u[end]) .== parent(Y.soil.ϑ_l)) == nelems
    # should be hydrostatic equilibrium at every layer, at each step:
    @test mean(
        sum([
            parent(saved_values.saveval[k].soil.ψ .+ coords.z)[:] .+ 10.0 for
            k in 2:1:50
        ]),
    ) < 1e-10

end


@testset "Soil Energy and Water RHS unit tests" begin
    ν = FT(0.495)
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.1)
    rre_params_on =
        Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r)
    rre_params_off =
        Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, 0.0, S_s, θ_r)

    κ = FT(10.0)
    ρc_s = FT(3e6)
    heat_params_on = Soil.HeatParameters{FT}(κ, ρc_s)
    heat_params_off = Soil.HeatParameters{FT}(0.0, ρc_s)

    zmax = FT(0)
    zmin = FT(-1)
    nelems = 200
    Δz = 0.005
    soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelems)
    top_flux_bc = FT(0.0)
    bot_flux_bc = FT(0.0)

    sources = ()
    boundary_fluxes = Soil.FluxBC{FT}(top_flux_bc, bot_flux_bc)


    # Test with only heat on, (hydraulic K = 0)
    soil_heat_on = Soil.SoilEnergyHydrology{FT}(;
        rre_param_set = rre_params_off,
        heat_param_set = heat_params_on,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_heat_on)

    # specify ICs
    function init_soil_heat!(Ysoil, z, rre_params, heat_params)
        ν = rre_params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0
        Ysoil.soil.θ_i .= 0.0
        T = 280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0
        @. Ysoil.soil.ρe_int =
            volumetric_internal_energy(0.0, heat_params.ρc_s, T)
    end

    init_soil_heat!(
        Y,
        coords.z,
        soil_heat_on.rre_param_set,
        soil_heat_on.heat_param_set,
    )
    soil_ode! = make_ode_function(soil_heat_on)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    F_face = 0.0
    F_below = -κ * (-Δz + 0.5)
    dY_top = -(F_face - F_below) / Δz
    F_top = -κ * ((-1.0 + Δz) + 0.5)
    dY_bot = -(F_top - F_face) / Δz
    expected = zeros(nelems) .+ κ
    expected[1] = dY_bot
    expected[end] = dY_top
    @test mean(abs.(expected .- parent(dY.soil.ρe_int))) / 1e7 < 1e-13 # to put on same scale as water
    @test mean(abs.(parent(dY.soil.ϑ_l))) < 1e-14
    @test mean(abs.(parent(dY.soil.θ_i))) < 1e-14


    # Test with water on (but κ = 0 for heat)
    soil_water_on = Soil.SoilEnergyHydrology{FT}(;
        rre_param_set = rre_params_on,
        heat_param_set = heat_params_off,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_water_on)

    # specify ICs
    function init_soil_water!(Ysoil, z, rre_params, heat_params)
        ν = rre_params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0 .+ ν / 4.0 .* (z .+ 0.5) .^ 2.0
        Ysoil.soil.θ_i .= 0.0
        @. Ysoil.soil.ρe_int =
            volumetric_internal_energy(0.0, heat_params.ρc_s, 280.0)
    end

    init_soil_water!(
        Y,
        coords.z,
        soil_water_on.rre_param_set,
        soil_water_on.heat_param_set,
    )
    soil_ode! = make_ode_function(soil_water_on)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    function dKdθ(θ::FT)::FT where {FT}
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1
            f::FT = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
            f1::FT = f^2.0 / 2.0 / S^0.5
            f2::FT =
                2 * S^(1 / vg_m - 1 / 2) * f / (1 - S^(1 / vg_m))^(1.0 - vg_m)
            return (f1 + f2) * Ksat / (ν - θ_r)
        else
            return 0.0
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
            return 0.0
        end
    end

    function K(θ::FT)::FT where {FT}
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1
            f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
            return Ksat * f^2.0 * sqrt(S)
        else
            return Ksat
        end

    end

    function dθdz(z::FT)::FT where {FT}
        return ν / 2.0 * (z + 0.5)
    end


    function d2θdz2(z::FT)::FT where {FT}
        return ν / 2.0
    end

    θ = parent(Y.soil.ϑ_l)
    Z = parent(coords.z)
    flux_c = (@. -K(θ) * (1.0 + dψdθ(θ) * dθdz(Z)))
    flux_f = (flux_c[2:end] .+ flux_c[1:(end - 1)]) ./ 2.0
    flux_f = vcat([0.0], flux_f, [0.0])
    expected = -(flux_f[2:end] - flux_f[1:(end - 1)]) ./ Δz
    @test mean(abs.(expected .- parent(dY.soil.ϑ_l))) < 1e-13
    @test mean(abs.(parent(dY.soil.θ_i))) < 1e-14

    ρe_int_l = volumetric_internal_energy_liq(280.0)
    flux_c = (@. -K(θ) * ρe_int_l * (1.0 + dψdθ(θ) * dθdz(Z)))
    flux_f = (flux_c[2:end] .+ flux_c[1:(end - 1)]) ./ 2.0
    flux_f = vcat([0.0], flux_f, [0.0])
    expected = -(flux_f[2:end] - flux_f[1:(end - 1)]) ./ Δz

    @test mean(abs.(expected .- parent(dY.soil.ρe_int))) < 1e-6


    ### Test with both off
    soil_both_off = Soil.SoilEnergyHydrology{FT}(;
        rre_param_set = rre_params_off,
        heat_param_set = heat_params_off,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_both_off)

    # specify ICs
    function init_soil_off!(Ysoil, z, rre_params, heat_params)
        ν = rre_params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0
        Ysoil.soil.θ_i .= 0.0
        T = 280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0
        @. Ysoil.soil.ρe_int =
            volumetric_internal_energy(0.0, heat_params.ρc_s, T)
    end

    init_soil_off!(
        Y,
        coords.z,
        soil_both_off.rre_param_set,
        soil_both_off.heat_param_set,
    )
    soil_ode! = make_ode_function(soil_both_off)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    @test mean(abs.(parent(dY.soil.ρe_int))) < 1e-14
    @test mean(abs.(parent(dY.soil.ϑ_l))) < 1e-14
    @test mean(abs.(parent(dY.soil.θ_i))) < 1e-14


    ### Test with both on


    soil_both_on = Soil.SoilEnergyHydrology{FT}(;
        rre_param_set = rre_params_on,
        heat_param_set = heat_params_on,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_both_on)

    # specify ICs
    function init_soil_on!(Ysoil, z, rre_params, heat_params)
        ν = rre_params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0 .+ ν / 4.0 .* (z .+ 0.5) .^ 2.0
        Ysoil.soil.θ_i .= 0.0
        T = 280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0

        @. Ysoil.soil.ρe_int =
            volumetric_internal_energy(0.0, heat_params.ρc_s, T)
    end

    init_soil_on!(
        Y,
        coords.z,
        soil_both_on.rre_param_set,
        soil_both_on.heat_param_set,
    )
    soil_ode! = make_ode_function(soil_both_on)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    @test mean(abs.(parent(dY.soil.θ_i))) < 1e-14

    θ = parent(Y.soil.ϑ_l)
    Z = parent(coords.z)
    flux_c = (@. -K(θ) * (1.0 + dψdθ(θ) * dθdz(Z)))
    flux_f = (flux_c[2:end] .+ flux_c[1:(end - 1)]) ./ 2.0
    flux_f = vcat([0.0], flux_f, [0.0])
    expected = -(flux_f[2:end] - flux_f[1:(end - 1)]) ./ Δz
    @test mean(abs.(expected .- parent(dY.soil.ϑ_l))) < 1e-13

    temp = 280.0 .+ 0.5 .* (Z .+ 0.5) .^ 2.0
    ρe_int_l = volumetric_internal_energy_liq.(temp)
    flux_c = (@. -K(θ) * ρe_int_l * (1.0 + dψdθ(θ) * dθdz(Z)))
    flux_f = (flux_c[2:end] .+ flux_c[1:(end - 1)]) ./ 2.0
    flux_f = vcat([0.0], flux_f, [0.0])
    part_one = -(flux_f[2:end] - flux_f[1:(end - 1)]) ./ Δz

    F_face = 0.0
    F_below = -κ * (-Δz + 0.5)
    dY_top = -(F_face - F_below) / Δz
    F_top = -κ * ((-1.0 + Δz) + 0.5)
    dY_bot = -(F_top - F_face) / Δz
    part_two = zeros(nelems) .+ κ
    part_two[1] = dY_bot
    part_two[end] = dY_top

    expected = part_one + part_two
    @test mean(abs.(expected .- parent(dY.soil.ρe_int))) < 1e-6

end

using Test
using Statistics
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: PlantHydraulicsDomain
using ClimaLSM.PlantHydraulics
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64

@testset "Plant hydraulics model integration tests" begin
    K_sat_root = FT(1e-5) # (accelerated) see Kumar, 2008 and
    # Pierre's database for total global plant conductance (1/resistance) 
    # (https://github.com/yalingliu-cu/plant-strategies/blob/master/Product%20details.pdf)
    K_sat_stem = FT(1e-3)
    K_sat_leaf = FT(1e-3)
    vg_α = FT(0.24) # Fitted VG to Weibull curve parameters found in Venturas, 2018
    vg_n = FT(2)
    vg_m = FT(1) - FT(1) / vg_n
    ν = FT(0.495)
    S_s = FT(1e-3)
    root_depths = [FT(-1.0)] # m, rooting depth
    z_ground = FT(0.0)
    z_stem_top = FT(5.0)
    z_leaf_top = FT(10.0)
    z_stem_midpoint = FT(2.5)
    z_leaf_midpoint = FT(7.5)
    earth_param_set = create_lsm_parameters(FT)
    plant_hydraulics_domain = PlantHydraulicsDomain{FT}(
        root_depths,
        [z_ground, z_stem_top, z_leaf_top],
        [z_stem_midpoint, z_leaf_midpoint],
    )
    param_set =
        PlantHydraulics.PlantHydraulicsParameters{FT, typeof(earth_param_set)}(
            K_sat_root,
            K_sat_stem,
            K_sat_leaf,
            vg_α,
            vg_n,
            vg_m,
            ν,
            S_s,
            earth_param_set,
        )

    function leaf_transpiration(t::FT) where {FT}
        T = FT(0)
        T_0 = FT(0)
        T_f = FT(1e-5)
        if t < FT(60^2 * 5)
            T = T_0
        elseif t < FT(60^2 * 10)
            T = FT(-80 * T_f)
        else
            T = FT(T_f)
        end
        return T
    end

    p_soil0 = [FT(0.0)]
    transpiration =
        PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
    root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
    plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
        domain = plant_hydraulics_domain,
        parameters = param_set,
        root_extraction = root_extraction,
        transpiration = transpiration,
    )

    # Set system to hydrostatic equilibrium state by setting fluxes to zero, and setting LHS of both ODEs to 0
    function initial_rhs!(F, Y)
        T0 = FT(0)
        flux_in_stem = sum(
            flux.(
                root_depths,
                z_stem_midpoint,
                p_soil0,
                Y[1],
                vg_α,
                vg_n,
                vg_m,
                ν,
                S_s,
                K_sat_root,
                K_sat_stem,
            ),
        )
        flux_out_stem = flux(
            z_stem_midpoint,
            z_leaf_midpoint,
            Y[1],
            Y[2],
            vg_α,
            vg_n,
            vg_m,
            ν,
            S_s,
            K_sat_stem,
            K_sat_leaf,
        )

        F[1] = flux_in_stem - T0
        F[2] = flux_out_stem - T0
    end

    soln = nlsolve(initial_rhs!, [-0.03, -0.02])
    p_stem_0 = soln.zero[1]
    p_leaf_0 = soln.zero[2]

    S_l_stem_0 =
        inverse_water_retention_curve(vg_α, vg_n, vg_m, p_stem_0, ν, S_s)
    S_l_leaf_0 =
        inverse_water_retention_curve(vg_α, vg_n, vg_m, p_leaf_0, ν, S_s)

    ϑ_l_stem_0 = augmented_liquid_fraction(ν, S_l_stem_0)
    ϑ_l_leaf_0 = augmented_liquid_fraction(ν, S_l_leaf_0)

    ϑ_l_0 = [ϑ_l_stem_0, ϑ_l_leaf_0]
    Y, p, coords = initialize(plant_hydraulics)

    Y.vegetation.ϑ_l .= ϑ_l_0

    plant_hydraulics_ode! = make_ode_function(plant_hydraulics)

    hour = FT(60^2)
    day = FT(24 * 60^2)
    t0 = FT(0)
    tf = FT(hour * 24 * 4)
    dt = FT(1)

    prob = ODEProblem(plant_hydraulics_ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(), dt = dt)

    dY = similar(Y)
    plant_hydraulics_ode!(dY, Y, p, 0.0)
    @test sqrt(mean(dY.vegetation.ϑ_l .^ 2.0)) < 1e-8 # starts in equilibrium

    ϑ_l_stem = reduce(hcat, sol.u)[1, :]
    ϑ_l_leaf = reduce(hcat, sol.u)[2, :]

    S_l_stem = effective_saturation.(ν, ϑ_l_stem)
    S_l_leaf = effective_saturation.(ν, ϑ_l_leaf)

    p_stem = water_retention_curve.(vg_α, vg_n, vg_m, S_l_stem, ν, S_s)
    p_leaf = water_retention_curve.(vg_α, vg_n, vg_m, S_l_leaf, ν, S_s)

end

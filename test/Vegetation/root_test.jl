#ideas: plot transpiration(t); theta(t); p(t) see if looks normal


using Test
using Statistics
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
import CLIMAParameters as CP
using Plots

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: RootDomain
using ClimaLSM.Roots
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64

#@testset "Root model integration tests" begin
    a_root = FT(0.3968) # modified logistic function parameters fitted to VG, for plant hydraulics model
    b_root = FT(7.692)
    a_stem = FT(0.3968) # modified logistic function parameters fitted to VG, for plant hydraulics model
    b_stem = FT(7.692)
    a_leaf = FT(0.3968) # modified logistic function parameters fitted to VG, for plant hydraulics model
    b_leaf = FT(7.692)
    K_sat_root = FT(0.0443 / 3600 / 100) # units of [m s-1]
    K_sat_stem = FT(0.0443 / 3600 / 100) # units of [m s-1]
    K_sat_leaf = FT(0.0443 / 3600 / 100) # units of [m s-1]
    vg_α = FT(2.6) # units of [m-1]
    vg_n = FT(2.0)
    vg_m = FT(1) - FT(1) / vg_n
    ν = FT(0.495)
    S_s = FT(1e-3) # units of [m-1]
    root_depths = [FT(-1.0)] # m, rooting depth
    z_ground = FT(0.0)
    z_stem_top = FT(5.0)
    z_leaf_top = FT(10.0)
    z_stem_midpoint = FT(2.5)
    z_leaf_midpoint = FT(7.5)
    earth_param_set = create_lsm_parameters(FT)
    root_domain = RootDomain{FT}(root_depths, [z_ground, z_stem_top, z_leaf_top], [z_stem_midpoint, z_leaf_midpoint])
    param_set = Roots.RootsParameters{FT, typeof(earth_param_set)}(
        a_root,
        b_root,
        a_stem,
        b_stem,
        a_leaf,
        b_leaf,
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
        T = FT(0.0)
        #=
        T_0 = FT(1e-7) 
        if t < FT(500)
            T = T_0
        elseif t < FT(1000)
            T = FT(10 * (T_0 / 5) * (t - 500) / 500 + T_0)
        else
            T = FT(10 * (T_0 / 5) + T_0)
        end
        =#
        return T
    end

    p_soil0 = [FT(-2)]
    transpiration =
        PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
    root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
    roots = Roots.RootsModel{FT}(;
        domain = root_domain,
        parameters = param_set,
        root_extraction = root_extraction,
        transpiration = transpiration,
    )

    #= Set system to equilibrium state by setting LHS of both ODEs to 0
    function f!(F, Y)
        T0 = FT(1.0000e-7)
        flux_in_stem = sum(
            flux.(
                root_depths,
                z_stem_midpoint,
                p_soil0,
                Y[1],
                a_root, 
                a_stem, 
                b_root, 
                b_stem, 
                K_sat_root, 
                K_sat_stem),
        )
        flux_out_stem = flux(
            z_stem_midpoint,
            z_leaf_midpoint,
            Y[1],
            Y[2],
            a_stem, 
            a_leaf, 
            b_stem, 
            b_leaf, 
            K_sat_stem, 
            K_sat_leaf,
        )

        F[1] = flux_in_stem - T0
        F[2] = flux_out_stem - T0
    end

    soln = nlsolve(f!, [0.03,0.02])
    p_stem_0 = soln.zero[1]
    p_leaf_0 = soln.zero[2]

    @show(p_stem_0)
    @show(p_leaf_0)

    ϑ_l_stem_0 = absolute_pressure_to_ϑ_l(vg_α,vg_n,vg_m,p_stem_0,ν,S_s)
    ϑ_l_leaf_0 = absolute_pressure_to_ϑ_l(vg_α,vg_n,vg_m,p_leaf_0,ν,S_s)
   
    ϑ_l_0 = [ϑ_l_stem_0, ϑ_l_leaf_0]
    =#

    ϑ_l_0 = FT(0.494) .+ zeros(2)
    Y, p, coords = initialize(roots)

    Y.vegetation.ϑ_l .= ϑ_l_0

    root_ode! = make_ode_function(roots)

    t0 = FT(0)
    tf = FT(60 * 60.0 * 30)
    dt = FT(1)

    prob = ODEProblem(root_ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(), dt = dt)

    dY = similar(Y)
    root_ode!(dY, Y, p, 0.0)
    #@test sqrt(mean(dY.vegetation.ϑ_l .^ 2.0)) < 1e-8 # starts in equilibrium

    ϑ_l_stem = reduce(hcat, sol.u)[1, :]
    ϑ_l_leaf = reduce(hcat, sol.u)[2, :]

    ϑ_l_of_t = plot(0:dt:tf,ϑ_l_stem)
    plot!(0:dt:tf,ϑ_l_leaf)

    p_stem = ϑ_l_to_absolute_pressure.(vg_α,vg_n,vg_m,ϑ_l_stem,ν,S_s)
    p_leaf = ϑ_l_to_absolute_pressure.(vg_α,vg_n,vg_m,ϑ_l_leaf,ν,S_s)

    p_of_t = plot(0:dt:tf,p_stem)
    plot!(p_of_t,0:dt:tf,p_leaf)

    flux_in_stem = flux.(root_depths,z_stem_midpoint,p_soil0,p_stem,a_root,a_stem,b_root,b_stem,K_sat_root,K_sat_stem)
    flux_out_stem = flux.(z_stem_midpoint,z_leaf_midpoint,p_stem,p_leaf,a_stem,a_leaf,b_stem,b_leaf,K_sat_stem,K_sat_leaf)

    flux_of_t = plot(0:dt:tf,flux_in_stem)
    plot!(flux_of_t,0:dt:tf,flux_out_stem)

    savefig(ϑ_l_of_t,"ϑ_l_of_t.png")
    savefig(p_of_t,"p_of_t.png")
    savefig(flux_of_t,"flux_of_t.png")

    #=
    function f2!(F, Y)
        p_soilf = p_soil0
        Tf = FT(1e-7 .* 3.0)
        flux_in_stem = sum(
            flux.(
                root_depths,
                z_stem_midpoint,
                p_soilf,
                Y[1],
                a_root,
                a_stem,
                b_root,
                b_stem,
                K_sat_root,
                K_sat_stem,
            ),        
        )
        flux_out_stem = flux(
            z_stem_midpoint,
            z_leaf_midpoint,
            Y[1],
            Y[2],
            a_stem,
            a_leaf,
            b_stem,
            b_leaf,
            K_sat_stem,
            K_sat_leaf,
        )
        F[1] = flux_in_stem - Tf
        F[2] = flux_out_stem - Tf
        @show(flux_in_stem)
        @show(flux_out_stem)
    end

    # Check that the final state is in the new equilibrium
    soln = nlsolve(f2!, [0.03, 0.02]; ftol = 1e-10)
    p_stem_f = soln.zero[1]
    @show(p_stem_f)
    p_leaf_f = soln.zero[2]
    @show(p_leaf_f)
    @test abs(p_stem_f - p_stem[end]) < 1e-10
    @test abs(p_leaf_f - p_leaf[end]) < 1e-10

    @show(p_stem[end])
    @show(p_leaf[end])
    =#

#end

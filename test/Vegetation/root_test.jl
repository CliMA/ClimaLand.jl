using Test
using Statistics
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Midpoint
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: RootDomain
using ClimaLSM.Roots
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64

@testset "Root model integration tests" begin
    a_variables = (root = FT(13192), stem = FT(515.5605), leaf = FT(515.5605))
    b_variables = (root = FT(2.1079), stem = FT(0.9631), leaf = FT(0.9631))
    size_reservoir_moles = (stem = FT(11000.8837), leaf = FT(16766.2790))
    K_max_moles = (root = FT(12.9216), stem = FT(3.4415), leaf = FT(3.4415))
    z_root_depths = [FT(-1.0)] # m, rooting depth
    delta_z = FT(4.0) # height of compartments
    n_stem = Int64(6)
    n_leaf = Int64(5)
    earth_param_set = create_lsm_parameters(FT)

    root_domain = RootDomain(z_root_depths, n_stem, n_leaf, delta_z)
    param_set = Roots.RootsParameters{FT, typeof(earth_param_set)}(
        a_variables,
        b_variables,
        size_reservoir_moles,
        K_max_moles,
        earth_param_set,
    )

    function leaf_transpiration(t::ft) where {ft}
        mass_mole_water = ft(0.018)
        T = ft(0.0)
        T_0 = ft(0.01 / mass_mole_water)
        if t < ft(500)
            T = T_0
        elseif t < ft(1000)
            T = ft(10 * (T_0 / 5) * (t - 500) / 500 + T_0)
        else
            T = ft(10 * (T_0 / 5) * 500 / 500 + T_0)
        end
        return T
    end

    p_soil0 = [FT(-0.02)]                                                          #pressure at the bottom of the root
    transpiration =
        PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
    root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
    roots = Roots.RootsModel{FT}(;
        domain = root_domain,
        parameters = param_set,
        root_extraction = root_extraction,
        transpiration = transpiration,
    )


    # testing many pieces at once (intergration test)


    function initial_rhs!(F, Y)
        T0 = leaf_transpiration(0.0)
        for i in 1:(n_stem + n_leaf)
            if i == 1
                flow_in = sum(
                    flow.(
                        z_root_depths,
                        root_domain.compartment_heights[i],
                        p_soil0,
                        Y[i],
                        a_variables.root,
                        b_variables.root,
                        K_max_moles.root,
                    ),
                )
            else
                flow_in = flow(
                    root_domain.compartment_heights[i - 1],
                    root_domain.compartment_heights[i],
                    Y[i - 1],
                    Y[i],
                    a_variables[root_domain.compartment_labels[i]],
                    b_variables[root_domain.compartment_labels[i]],
                    K_max_moles[root_domain.compartment_labels[i]],
                )
            end
            F[i] = flow_in - T0
        end
    end

    soln = nlsolve(
        initial_rhs!,
        [-1.0, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45],
    )  #solves function and uses the numbers as guesses to get closer to answer (initial guess)
    theta_0 = p_to_theta.(soln.zero) #pressures where you are at equilibrium
    sizes =
        getproperty.(Ref(size_reservoir_moles), root_domain.compartment_labels)
    y0 = @. FT(theta_0 * sizes)
    Y, p, coords = initialize(roots)
    Y.vegetation.rwc .= y0

    root_ode! = make_ode_function(roots)

    dY = similar(Y)
    root_ode!(dY, Y, p, 0.0)
    @test sqrt(mean(dY.vegetation.rwc .^ 2.0)) < 1e-8 # starts in equilibrium




    t0 = FT(0)
    tf = FT(60 * 60.0 * 10)
    dt = FT(1)

    prob = ODEProblem(root_ode!, Y, (t0, tf), p)
    sol = solve(prob, Midpoint(), dt = dt)


    # extract the final state at the simulation end
    y_final = sol.u[end - 1].vegetation.rwc

    # convert back to pressure
    theta = y_final ./ sizes
    p_compartment = theta_to_p.(theta)


    function final_rhs!(F, Y)
        p_soilf = p_soil0
        Tf = 0.01 / 0.018 .* 3.0
        for i in 1:(n_stem + n_leaf)
            if i == 1
                flow_in = sum(
                    flow.(
                        z_root_depths,
                        root_domain.compartment_heights[i],
                        p_soilf,
                        Y[i],
                        a_variables.root,
                        b_variables.root,
                        K_max_moles.root,
                    ),
                )
            else
                flow_in = flow(
                    root_domain.compartment_heights[i - 1],
                    root_domain.compartment_heights[i],
                    Y[i - 1],
                    Y[i],
                    a_variables[root_domain.compartment_labels[i]],
                    b_variables[root_domain.compartment_labels[i]],
                    K_max_moles[root_domain.compartment_labels[i]],
                )
            end
            F[i] = flow_in - Tf
        end
    end


    # Check that the final state is in the new equilibrium
    soln = nlsolve(
        final_rhs!,
        [-1.0, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45];
        ftol = 3.01e-3,
    )
    p_f = soln.zero
    equilib = p_f .- p_compartment
    equilib = broadcast(abs, equilib)

    @test equilib[1] < 3.01e-3
    @test equilib[2] < 3.01e-3
    @test equilib[3] < 3.01e-3
    @test equilib[4] < 3.01e-3
    @test equilib[5] < 3.01e-3
    @test equilib[6] < 3.01e-3
    @test equilib[7] < 3.01e-3
    @test equilib[8] < 3.01e-3
    @test equilib[9] < 3.01e-3
    @test equilib[10] < 3.01e-3
    @test equilib[11] < 3.01e-3
end



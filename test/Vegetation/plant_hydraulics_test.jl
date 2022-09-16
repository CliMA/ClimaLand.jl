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
    K_sat = (root = FT(1e-5), stem = FT(1e-3), leaf = FT(1e-3))# (accelerated) see Kumar, 2008 and
    # Pierre's database for total global plant conductance (1/resistance) 
    # (https://github.com/yalingliu-cu/plant-strategies/blob/master/Product%20details.pdf)
    vg_α = FT(0.24) # Fitted VG to Weibull curve parameters found in Venturas, 2018
    vg_n = FT(2)
    vg_m = FT(1) - FT(1) / vg_n
    ν = FT(0.495)
    S_s = FT(1e-3)
    z_root_depths = [FT(-1.0)] # m, rooting depth
    Δz = FT(1.0) # height of compartments
    n_stem = Int64(6)
    n_leaf = Int64(5)
    earth_param_set = create_lsm_parameters(FT)

    plant_hydraulics_domain =
        PlantHydraulicsDomain(z_root_depths, n_stem, n_leaf, Δz)
    param_set =
        PlantHydraulics.PlantHydraulicsParameters{FT, typeof(earth_param_set)}(
            K_sat,
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
        for i in 1:(n_leaf + n_stem)
            if i == 1
                flux_in = sum(
                    flux.(
                        z_root_depths,
                        plant_hydraulics_domain.compartment_midpoints[i],
                        p_soil0,
                        Y[i],
                        vg_α,
                        vg_n,
                        vg_m,
                        ν,
                        S_s,
                        K_sat[:root],
                        K_sat[plant_hydraulics_domain.compartment_labels[i]],
                    ),
                )
            else
                flux_in = flux(
                    plant_hydraulics_domain.compartment_midpoints[i - 1],
                    plant_hydraulics_domain.compartment_midpoints[i],
                    Y[i - 1],
                    Y[i],
                    vg_α,
                    vg_n,
                    vg_m,
                    ν,
                    S_s,
                    K_sat[plant_hydraulics_domain.compartment_labels[i - 1]],
                    K_sat[plant_hydraulics_domain.compartment_labels[i]],
                )
            end
            F[i] = flux_in - T0
        end
    end

    soln = nlsolve(initial_rhs!, Vector(-0.03:0.01:0.07))

    S_l = inverse_water_retention_curve.(vg_α, vg_n, vg_m, soln.zero, ν, S_s)

    ϑ_l_0 = augmented_liquid_fraction.(ν, S_l)

    Y, p, coords = initialize(plant_hydraulics)

    Y.vegetation.ϑ_l .= ϑ_l_0

    plant_hydraulics_ode! = make_ode_function(plant_hydraulics)

    dY = similar(Y)
    plant_hydraulics_ode!(dY, Y, p, 0.0)
    @test sqrt(mean(dY.vegetation.ϑ_l .^ 2.0)) < 1e-8 # starts in equilibrium
end

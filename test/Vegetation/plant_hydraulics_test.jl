using Test
using Statistics
using DiffEqCallbacks
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, RK4
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
    # Parameters are the same as the ones used in the Ozark tutorial
    RAI = FT(1) # m2/m2
    SAI = FT(1) # m2/m2
    LAI = FT(1) # m2/m2
    area_index = (root = RAI, stem = SAI, leaf = LAI)
    K_sat_plant = 1.8e-8 # m/s. Typical conductivity range is [1e-8, 1e-5] m/s. See Kumar, 2008 and
    # Pierre Gentine's database for total global plant conductance (1/resistance) 
    # (https://github.com/yalingliu-cu/plant-strategies/blob/master/Product%20details.pdf)
    K_sat_root = FT(K_sat_plant) # m/s
    K_sat_stem = FT(K_sat_plant)
    K_sat_leaf = FT(K_sat_plant)
    K_sat = (root = K_sat_root, stem = K_sat_stem, leaf = K_sat_leaf)
    plant_vg_α = FT(0.002) # 1/m
    plant_vg_n = FT(4.2) # unitless
    plant_vg_m = FT(1) - FT(1) / plant_vg_n
    plant_ν = FT(0.7) # m3/m3
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    root_depths = -Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0 # 1st element is the deepest root depth 
    function root_distribution(z::T) where {T}
        return T(1.0 / 0.5) * exp(z / T(0.5)) # (1/m)
    end
    Δz = FT(1.0) # height of compartments
    n_stem = Int64(6) # number of stem elements
    n_leaf = Int64(5) # number of leaf elements
    earth_param_set = create_lsm_parameters(FT)

    plant_hydraulics_domain =
        PlantHydraulicsDomain(root_depths, n_stem, n_leaf, Δz)
    param_set =
        PlantHydraulics.PlantHydraulicsParameters{FT, typeof(earth_param_set)}(
            area_index,
            K_sat,
            plant_vg_α,
            plant_vg_n,
            plant_vg_m,
            plant_ν,
            plant_S_s,
            root_distribution,
            earth_param_set,
        )

    function leaf_transpiration(t::FT) where {FT}
        T = FT(0)
    end

    ψ_soil0 = [FT(0.0)]
    transpiration =
        PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
    root_extraction = PrescribedSoilPressure{FT}((t::FT) -> ψ_soil0)
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
                        root_depths,
                        plant_hydraulics_domain.compartment_midpoints[i],
                        ψ_soil0,
                        Y[i],
                        plant_vg_α,
                        plant_vg_n,
                        plant_vg_m,
                        plant_ν,
                        plant_S_s,
                        K_sat[:root],
                        K_sat[:stem],
                    ) .* root_distribution.(root_depths) .* (
                        vcat(root_depths, [0.0])[2:end] -
                        vcat(root_depths, [0.0])[1:(end - 1)]
                    ),
                )
            else
                flux_in = flux(
                    plant_hydraulics_domain.compartment_midpoints[i - 1],
                    plant_hydraulics_domain.compartment_midpoints[i],
                    Y[i - 1],
                    Y[i],
                    plant_vg_α,
                    plant_vg_n,
                    plant_vg_m,
                    plant_ν,
                    plant_S_s,
                    K_sat[plant_hydraulics_domain.compartment_labels[i - 1]],
                    K_sat[plant_hydraulics_domain.compartment_labels[i]],
                )
            end
            F[i] = flux_in - T0
        end
    end

    soln = nlsolve(initial_rhs!, Vector(-0.03:0.01:0.07))

    S_l =
        inverse_water_retention_curve.(
            plant_vg_α,
            plant_vg_n,
            plant_vg_m,
            soln.zero,
            plant_ν,
            plant_S_s,
        )

    ϑ_l_0 = augmented_liquid_fraction.(plant_ν, S_l)

    Y, p, coords = initialize(plant_hydraulics)

    Y.vegetation.ϑ_l .= ϑ_l_0

    plant_hydraulics_ode! = make_ode_function(plant_hydraulics)

    dY = similar(Y)
    plant_hydraulics_ode!(dY, Y, p, 0.0)
    @test sqrt(mean(dY.vegetation.ϑ_l .^ 2.0)) < 1e-8 # starts in equilibrium
end

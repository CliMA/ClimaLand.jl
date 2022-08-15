using Test
using DiffEqCallbacks
using UnPack
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using ClimaLSM
using ClimaLSM.Domains: Column, PlantHydraulicsDomain
using ClimaLSM.Soil
using ClimaLSM.PlantHydraulics
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64
@testset " Soil plant hydrology LSM integration test" begin
    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    earth_param_set = create_lsm_parameters(FT)
    K_sat_root = FT(1e-5) # Pierre's database (1e-6) (https://github.com/yalingliu-cu/plant-strategies/blob/master/Product%20details.pdf)
    K_sat_stem = FT(1e-3)
    K_sat_leaf = FT(1e-3)
    vg_α = FT(0.24)
    vg_n = FT(2)
    vg_m = FT(1) - FT(1) / vg_n
    ν = FT(0.495)
    S_s = FT(1e-3)
    z_root_depths = -Array(1:1:20.0) ./ 20.0 * 3.0 .+ 0.15 / 2.0
    z_ground = FT(0.0)
    z_stem_top = FT(6.0)
    z_leaf_top = FT(12) # height of leaf
    # currently hardcoded to match the soil coordinates. this has to
    # be fixed eventually.
    z_stem_midpoint = FT(3)
    z_leaf_midpoint = FT(9)
    plant_hydraulics_domain = PlantHydraulicsDomain{FT}(
        z_root_depths,
        [z_ground, z_stem_top, z_leaf_top],
        [z_stem_midpoint, z_leaf_midpoint],
    )
    plant_hydraulics_ps =
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

    zmin = FT(-3.0)
    zmax = FT(0.0)
    nelements = 20
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelements)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)
    soil_args = (domain = soil_domain, parameters = soil_ps)
    plant_hydraulics_args =
        (domain = plant_hydraulics_domain, parameters = plant_hydraulics_ps)
    land = SoilPlantHydrologyModel{FT}(;
        soil_model_type = Soil.RichardsModel{FT},
        soil_args = soil_args,
        vegetation_model_type = PlantHydraulics.PlantHydraulicsModel{FT},
        vegetation_args = plant_hydraulics_args,
    )
    Y, p, coords = initialize(land)
    # specify ICs
    function init_soil!(Ysoil, z, params)
        function hydrostatic_profile(
            z::FT,
            params::RichardsParameters{FT},
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r = params
            #unsaturated zone only, assumes water table starts at z_∇
            z_∇ = FT(-3)# matches zmin
            S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
            ϑ_l = S * (ν - θ_r) + θ_r
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
    end
    init_soil!(Y, coords.subsurface.z, land.soil.parameters)

    ## soil is at total ψ+z = -3.0 #m
    ## Want (ψ+z)_plant = (ψ+z)_soil 
    p_stem_0 = (-3.0 - z_ground)
    p_leaf_0 = (-3.0 - z_leaf_top)

    ϑ_l_stem_0 =
        inverse_water_retention_curve(vg_α, vg_n, vg_m, p_stem_0, ν, S_s)
    ϑ_l_leaf_0 =
        inverse_water_retention_curve(vg_α, vg_n, vg_m, p_leaf_0, ν, S_s)
    ϑ_l_0 = [ϑ_l_stem_0, ϑ_l_leaf_0]

    Y.vegetation.ϑ_l .= ϑ_l_0

    ode! = make_ode_function(land)
    t0 = FT(0)
    tf = FT(2)
    dt = FT(1)
    cb = SavingCallback((u, t, integrator) -> integrator.p, saved_values)
    prob = ODEProblem(ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(), dt = dt, callback = cb)
    #Currently just testing to make sure it runs, but need to have a better test suite.
end

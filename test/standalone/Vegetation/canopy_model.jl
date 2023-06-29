using Test
import CLIMAParameters as CP
using ClimaCore
using Thermodynamics
using Insolation
using Dates
using ClimaLSM.Canopy
using ClimaLSM
using ClimaLSM: PrescribedAtmosphere, PrescribedRadiativeFluxes
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
using ClimaLSM.Domains: Point

import ClimaLSM
import ClimaLSM.Parameters as LSMP

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
@testset "Canopy software pipes" begin
    FT = Float32
    domain = Point(; z_sfc = FT(0.0))

    RTparams = BeerLambertParameters{FT}()
    photosynthesis_params = FarquharParameters{FT}(C3();)
    stomatal_g_params = MedlynConductanceParameters{FT}()

    stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
    photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
    rt_model = BeerLambertModel{FT}(RTparams)

    earth_param_set = create_lsm_parameters(FT)
    LAI = FT(8.0) # m2 [leaf] m-2 [ground]
    z_0m = FT(2.0) # m, Roughness length for momentum - value from tall forest ChatGPT
    z_0b = FT(0.1) # m, Roughness length for scalars - value from tall forest ChatGPT
    h_c = FT(20.0) # m, canopy height
    h_sfc = FT(20.0) # m, canopy height
    h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
    shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
        LAI,
        h_c,
        z_0m,
        z_0b,
        earth_param_set,
    )
    lat = FT(0.0) # degree
    long = FT(-180) # degree

    function zenith_angle(
        t::FT,
        orbital_data;
        latitude = lat,
        longitude = long,
        insol_params = earth_param_set.insol_params,
    ) where {FT}
        return FT(
            instantaneous_zenith_angle(
                DateTime(t),
                orbital_data,
                longitude,
                latitude,
                insol_params,
            )[1],
        )
    end

    function shortwave_radiation(
        t::FT;
        latitude = lat,
        longitude = long,
        insol_params = earth_param_set.insol_params,
    ) where {FT}
        return FT(1000) # W/m^2
    end

    function longwave_radiation(t::FT) where {FT}
        return FT(200) # W/m^2
    end

    u_atmos = t -> eltype(t)(10) #m.s-1

    liquid_precip = (t) -> eltype(t)(0) # m
    snow_precip = (t) -> eltype(t)(0) # m
    T_atmos = t -> eltype(t)(290) # Kelvin
    q_atmos = t -> eltype(t)(0.001) # kg/kg
    P_atmos = t -> eltype(t)(1e5) # Pa
    h_atmos = h_int # m
    c_atmos = (t) -> eltype(t)(4.11e-4) # mol/mol
    atmos = PrescribedAtmosphere(
        liquid_precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        h_atmos;
        c_co2 = c_atmos,
    )
    radiation = PrescribedRadiativeFluxes(
        FT,
        shortwave_radiation,
        longwave_radiation;
        θs = zenith_angle,
        orbital_data = Insolation.OrbitalData(),
    )

    # Plant Hydraulics
    RAI = FT(1)
    SAI = FT(0)
    area_index = (root = RAI, stem = SAI, leaf = LAI)
    K_sat_plant = FT(1.8e-8) # m/s
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    plant_ν = FT(0.7) # m3/m3
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
    root_depths = FT.(-Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0) # 1st element is the deepest root depth
    function root_distribution(z::T) where {T}
        return T(1.0 / 0.5) * exp(z / T(0.5)) # (1/m)
    end
    param_set = PlantHydraulics.PlantHydraulicsParameters(;
        area_index = area_index,
        ν = plant_ν,
        S_s = plant_S_s,
        root_distribution = root_distribution,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )
    Δz = FT(1.0) # height of compartments
    n_stem = Int64(0) # number of stem elements
    n_leaf = Int64(1) # number of leaf elements
    compartment_centers =
        FT.(
            Vector(
                range(
                    start = Δz / 2,
                    step = Δz,
                    stop = Δz * (n_stem + n_leaf) - (Δz / 2),
                ),
            ),
        )
    compartment_faces =
        FT.(
            Vector(
                range(start = 0.0, step = Δz, stop = Δz * (n_stem + n_leaf)),
            )
        )

    ψ_soil0 = FT(0.0)
    root_extraction = PrescribedSoilPressure{FT}((t::FT) -> ψ_soil0)

    plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
        parameters = param_set,
        root_extraction = root_extraction,
        root_depths = root_depths,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_surfaces = compartment_faces,
        compartment_midpoints = compartment_centers,
    )
    canopy = ClimaLSM.Canopy.CanopyModel{FT}(;
        parameters = shared_params,
        domain = domain,
        radiative_transfer = rt_model,
        photosynthesis = photosynthesis_model,
        conductance = stomatal_model,
        hydraulics = plant_hydraulics,
        atmos = atmos,
        radiation = radiation,
    )
    Y, p, coords = ClimaLSM.initialize(canopy)

    for component in ClimaLSM.Canopy.canopy_components(canopy)
        @test propertynames(getproperty(p.canopy, component)) ==
              ClimaLSM.auxiliary_vars(getproperty(canopy, component))
        @test propertynames(getproperty(Y.canopy, component)) ==
              ClimaLSM.prognostic_vars(getproperty(canopy, component))

        @test getproperty(auxiliary_types(canopy), component) ==
              auxiliary_types(getproperty(canopy, component))
        @test getproperty(auxiliary_vars(canopy), component) ==
              auxiliary_vars(getproperty(canopy, component))
        @test getproperty(prognostic_types(canopy), component) ==
              prognostic_types(getproperty(canopy, component))
        @test getproperty(prognostic_types(canopy), component) ==
              prognostic_types(getproperty(canopy, component))
    end
    Y.canopy.hydraulics[1] = plant_ν
    exp_tendency! = make_exp_tendency(canopy)
    t0 = FT(0.0)
    dY = similar(Y)
    exp_tendency!(dY, Y, p, t0)
    (evapotranspiration, shf, lhf) =
        canopy_surface_fluxes(canopy.atmos, canopy, Y, p, t0)

    @test p.canopy.hydraulics.fa[1] ≈ evapotranspiration

    # Penman-monteith
    Δ = FT(100 * (0.444017302 + (290 - 273.15) * 0.0286064092))
    Rn = shortwave_radiation(t0)
    G = FT(0.0)
    thermo_params = canopy.parameters.earth_param_set.thermo_params
    ts_in = construct_atmos_ts(atmos, t0, thermo_params)
    ρa = Thermodynamics.air_density(thermo_params, ts_in)
    cp = cp_m(thermo_params, Thermodynamics.PhasePartition.(q_atmos(t0)))

    es =
        Thermodynamics.saturation_vapor_pressure.(
            Ref(thermo_params),
            T_atmos(t0),
            Ref(Thermodynamics.Liquid()),
        )
    ea =
        Thermodynamics.partial_pressure_vapor.(
            Ref(thermo_params),
            P_atmos(t0),
            Thermodynamics.PhasePartition.(q_atmos(t0)),
        )

    VPD = es .- ea

    conditions = surface_fluxes(atmos, canopy, Y, p, t0) #Per unit m^2 of leaf
    # here is where we adjust evaporation for stomatal conductance = 1/r_sfc
    r_ae = 1 / (conditions.Ch * abs(atmos.u(t0))) # s/m
    ga = 1 / r_ae
    γ = FT(66)
    R = FT(LSMP.gas_constant(earth_param_set))
    gs = parent(
        ClimaLSM.Canopy.upscale_leaf_conductance.(
            p.canopy.conductance.gs,
            LAI,
            T_atmos(t0),
            R,
            P_atmos(t0),
        ),
    )[1]
    Lv = FT(2453e6) #J/m^3

    ET = penman_monteith(
        Δ, # Rate of change of saturation specific humidity with air temperature. (Pa K−1)
        Rn, # Net irradiance (W m−2), the external source of energy flux
        G, # Ground heat flux (W m−2)
        ρa, # Dry air density (kg m−3)
        cp, # Specific heat capacity of air (J kg−1 K−1)
        VPD, # vapor pressure deficit (Pa)
        ga, # Conductivity of air, atmospheric conductance (m s−1)
        γ, # Psychrometric constant (γ ≈ 66 Pa K−1)
        gs, # Conductivity of stoma, surface or stomatal conductance (m s−1)
        Lv, # Volumetric latent heat of vaporization. Energy required per water volume vaporized. (Lv = 2453 MJ m−3)
    )

    @test abs(
        (parent(evapotranspiration)[1] - ET) / parent(evapotranspiration)[1],
    ) < 0.15

    @test ClimaLSM.surface_evaporative_scaling(canopy, Y, p) == FT(1.0)
    @test ClimaLSM.surface_height(canopy, Y, p) == h_sfc
    T_sfc = canopy.atmos.T(t0)
    @test ClimaLSM.surface_temperature(canopy, Y, p, t0) == T_sfc
    ρ_sfc = ClimaLSM.surface_air_density(canopy.atmos, canopy, Y, p, t0, T_sfc)
    @test ClimaLSM.surface_specific_humidity(canopy, Y, p, T_sfc, ρ_sfc) ==
          Thermodynamics.q_vap_saturation_generic.(
        Ref(thermo_params),
        T_sfc,
        ρ_sfc,
        Ref(Thermodynamics.Liquid()),
    )
    @test ρ_sfc == compute_ρ_sfc.(Ref(thermo_params), Ref(ts_in), T_sfc)
    @test ClimaLSM.domain(canopy) == :surface
end

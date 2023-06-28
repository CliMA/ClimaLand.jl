using Test
using Statistics
using NLsolve
using ClimaCore
import CLIMAParameters as CP
using ClimaLSM
using ClimaLSM.Domains: Point, Plane
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM
using Insolation
using Dates

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

@testset "Plant hydraulics parameterizations" begin
    FT = Float32
    ν = FT(0.5)
    S_s = FT(1e-2)
    K_sat = FT(1.8e-8)
    ψ63 = FT(-4 / 0.0098)
    Weibull_param = FT(4)
    a = FT(0.05 * 0.0098)
    conductivity_model = PlantHydraulics.Weibull{FT}(K_sat, ψ63, Weibull_param)
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
    S_l = FT.([0.7, 0.9, 1.0, 1.1])
    @test all(
        @. PlantHydraulics.inverse_water_retention_curve(
            retention_model,
            PlantHydraulics.water_retention_curve(retention_model, S_l, ν, S_s),
            ν,
            S_s,
        ) == S_l
    )
    ψ = PlantHydraulics.water_retention_curve.(retention_model, S_l, ν, S_s)
    @test all(
        @. PlantHydraulics.hydraulic_conductivity(conductivity_model, ψ) ≈
           min(K_sat * exp(-(ψ / ψ63)^Weibull_param), K_sat)
    )
end

@testset "Plant hydraulics model integration tests" begin
    FT = Float64
    domains = [
        Point(; z_sfc = FT(0.0)),
        Plane(;
            xlim = (0.0, 1.0),
            ylim = (0.0, 1.0),
            nelements = (2, 2),
            periodic = (true, true),
            npolynomial = 1,
        ),
    ]
    RTparams = BeerLambertParameters{FT}()
    photosynthesis_params = FarquharParameters{FT}(C3();)
    stomatal_g_params = MedlynConductanceParameters{FT}()

    stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
    photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
    rt_model = BeerLambertModel{FT}(RTparams)

    earth_param_set = create_lsm_parameters(FT)
    LAI = FT(1.0) # m2 [leaf] m-2 [ground]
    z_0m = FT(2.0) # m, Roughness length for momentum
    z_0b = FT(0.1) # m, Roughness length for scalars
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

    for domain in domains
        # Parameters are the same as the ones used in the Ozark tutorial
        Δz = FT(1.0) # height of compartments
        n_stem = Int64(5) # number of stem elements
        n_leaf = Int64(6) # number of leaf elements
        SAI = FT(1) # m2/m2
        RAI = FT(1) # m2/m2
        area_index = (root = RAI, stem = SAI, leaf = LAI)
        K_sat_plant = 1.8e-8 # m/s.
        ψ63 = FT(-4 / 0.0098) # / MPa to m
        Weibull_param = FT(4) # unitless
        a = FT(0.05 * 0.0098) # 1/m
        plant_ν = FT(0.7) # m3/m3
        plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
        conductivity_model =
            PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
        retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
        root_depths = -Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0 # 1st element is the deepest root depth 
        function root_distribution(z::T) where {T}
            return T(1.0 / 0.5) * exp(z / T(0.5)) # (1/m)
        end
        compartment_midpoints = Vector(
            range(
                start = Δz / 2,
                step = Δz,
                stop = Δz * (n_stem + n_leaf) - (Δz / 2),
            ),
        )

        compartment_surfaces =
            Vector(range(start = 0.0, step = Δz, stop = Δz * (n_stem + n_leaf)))

        param_set = PlantHydraulics.PlantHydraulicsParameters(;
            area_index = area_index,
            ν = plant_ν,
            S_s = plant_S_s,
            root_distribution = root_distribution,
            conductivity_model = conductivity_model,
            retention_model = retention_model,
        )

        function leaf_transpiration(t::FT) where {FT}
            T = FT(1e-8) # m/s
        end

        ψ_soil0 = FT(0.0)
        transpiration =
            PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
        root_extraction = PrescribedSoilPressure{FT}((t::FT) -> ψ_soil0)

        plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
            parameters = param_set,
            root_extraction = root_extraction,
            transpiration = transpiration,
            root_depths = root_depths,
            n_stem = n_stem,
            n_leaf = n_leaf,
            compartment_surfaces = compartment_surfaces,
            compartment_midpoints = compartment_midpoints,
        )
        model = ClimaLSM.Canopy.CanopyModel{FT}(;
            parameters = shared_params,
            domain = domain,
            radiative_transfer = rt_model,
            photosynthesis = photosynthesis_model,
            conductance = stomatal_model,
            hydraulics = plant_hydraulics,
            atmos = atmos,
            radiation = radiation,
        )
        # Set system to hydrostatic equilibrium state by setting fluxes to zero, and setting LHS of both ODEs to 0
        function initial_compute_exp_tendency!(F, Y)
            T0A = FT(1e-8) * area_index[:leaf]
            for i in 1:(n_leaf + n_stem)
                if i == 1
                    fa =
                        sum(
                            flux.(
                                root_depths,
                                plant_hydraulics.compartment_midpoints[i],
                                ψ_soil0,
                                Y[i],
                                PlantHydraulics.hydraulic_conductivity(
                                    conductivity_model,
                                    ψ_soil0,
                                ),
                                PlantHydraulics.hydraulic_conductivity(
                                    conductivity_model,
                                    Y[i],
                                ),
                            ) .* root_distribution.(root_depths) .* (
                                vcat(root_depths, [0.0])[2:end] -
                                vcat(root_depths, [0.0])[1:(end - 1)]
                            ),
                        ) * area_index[:root]
                else
                    fa =
                        flux(
                            plant_hydraulics.compartment_midpoints[i - 1],
                            plant_hydraulics.compartment_midpoints[i],
                            Y[i - 1],
                            Y[i],
                            PlantHydraulics.hydraulic_conductivity(
                                conductivity_model,
                                ψ_soil0,
                            ),
                            PlantHydraulics.hydraulic_conductivity(
                                conductivity_model,
                                Y[i],
                            ),
                        ) * (
                            area_index[plant_hydraulics.compartment_labels[1]] +
                            area_index[plant_hydraulics.compartment_labels[2]]
                        ) / 2
                end
                F[i] = fa - T0A
            end
        end

        soln = nlsolve(
            initial_compute_exp_tendency!,
            Vector(-0.03:0.01:0.07);
            ftol = 1e-11,
        )

        S_l =
            inverse_water_retention_curve.(
                retention_model,
                soln.zero,
                plant_ν,
                plant_S_s,
            )

        ϑ_l_0 = augmented_liquid_fraction.(plant_ν, S_l)

        Y, p, coords = initialize(model)
        dY = similar(Y)
        for i in 1:(n_stem + n_leaf)
            Y.canopy.hydraulics.ϑ_l[i] .= ϑ_l_0[i]
            p.canopy.hydraulics.ψ[i] .= NaN
            p.canopy.hydraulics.fa[i] .= NaN
            dY.canopy.hydraulics.ϑ_l[i] .= NaN
        end
        canopy_exp_tendency! = make_exp_tendency(model)
        canopy_exp_tendency!(dY, Y, p, 0.0)

        m = similar(dY.canopy.hydraulics.ϑ_l[1])
        m .= FT(0)
        for i in 1:(n_stem + n_leaf)
            @. m += sqrt(dY.canopy.hydraulics.ϑ_l[i]^2.0)
        end
        #@test maximum(parent(m)) < 1e-11 # starts in equilibrium


        # repeat using the plant hydraulics model directly
        # make sure it agrees with what we get when use the canopy model ODE
        Y, p, coords = initialize(model)
        standalone_dY = similar(Y)
        for i in 1:(n_stem + n_leaf)
            Y.canopy.hydraulics.ϑ_l[i] .= ϑ_l_0[i]
            p.canopy.hydraulics.ψ[i] .= NaN
            p.canopy.hydraulics.fa[i] .= NaN
            standalone_dY.canopy.hydraulics.ϑ_l[i] .= NaN
        end
        update_aux! = make_update_aux(model)
        standalone_exp_tendency! =
            make_compute_exp_tendency(model.hydraulics, nothing)
        update_aux!(p, Y, 0.0)
        standalone_exp_tendency!(standalone_dY, Y, p, 0.0)

        m = similar(dY.canopy.hydraulics.ϑ_l[1])
        m .= FT(0)
        for i in 1:(n_stem + n_leaf)
            @. m += sqrt(
                (
                    dY.canopy.hydraulics.ϑ_l[i] -
                    standalone_dY.canopy.hydraulics.ϑ_l[i]
                )^2.0,
            )
        end
        @test sum(parent(m)) < eps(FT)
    end
end

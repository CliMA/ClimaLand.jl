using Test
import ClimaComms
ClimaComms.@import_required_backends
using Statistics
using NLsolve
using ClimaCore
import ClimaParams
using ClimaLand
using ClimaLand.Domains: Point, Plane
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import Insolation
using Dates


for FT in (Float32, Float64)
    @testset "LAI assertions, FT = $FT" begin
        LAI = FT.([0, 1])
        SAI = FT.([0, 1])
        RAI = FT.([0, 1])
        n_stem = [0, 1]
        n_leaf = [1]
        for L in LAI
            for S in SAI
                for R in RAI
                    for n_s in n_stem
                        for n_l in n_leaf
                            area_index = (root = R, stem = S, leaf = L)
                            if n_l == 0
                                @test_throws AssertionError ClimaLand.Canopy.PlantHydraulics.lai_consistency_check(
                                    n_s,
                                    n_l,
                                    area_index,
                                )
                            end

                            if (L > eps(FT) || S > eps(FT)) && R < eps(FT)
                                @test_throws AssertionError ClimaLand.Canopy.PlantHydraulics.lai_consistency_check(
                                    n_s,
                                    n_l,
                                    area_index,
                                )
                            end

                            if S > FT(0) && n_s == 0
                                @test_throws AssertionError ClimaLand.Canopy.PlantHydraulics.lai_consistency_check(
                                    n_s,
                                    n_l,
                                    area_index,
                                )
                            end

                            if S < eps(FT) && n_s > 0
                                @test_throws AssertionError ClimaLand.Canopy.PlantHydraulics.lai_consistency_check(
                                    n_s,
                                    n_l,
                                    area_index,
                                )
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Plant hydraulics parameterizations, FT = $FT" begin
        ν = FT(0.5)
        S_s = FT(1e-2)
        K_sat = FT(1.8e-8)
        ψ63 = FT(-4 / 0.0098)
        Weibull_param = FT(4)
        a = FT(0.05 * 0.0098)
        conductivity_model =
            PlantHydraulics.Weibull{FT}(K_sat, ψ63, Weibull_param)
        retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
        S_l = FT.([0.7, 0.9, 1.0, 1.1])
        @test all(
            @. PlantHydraulics.inverse_water_retention_curve(
                retention_model,
                PlantHydraulics.water_retention_curve(
                    retention_model,
                    S_l,
                    ν,
                    S_s,
                ),
                ν,
                S_s,
            ) == S_l
        )
        ψ = PlantHydraulics.water_retention_curve.(retention_model, S_l, ν, S_s)
        # @show @. abs(PlantHydraulics.hydraulic_conductivity(conductivity_model, ψ) -
        # min(K_sat * exp(-(ψ / ψ63)^Weibull_param), K_sat))
        @test all(
            @. abs(
                PlantHydraulics.hydraulic_conductivity(conductivity_model, ψ) -
                min(K_sat * exp(-(ψ / ψ63)^Weibull_param), K_sat),
            ) < 2 * eps(FT)
        )
    end

    @testset "Plant hydraulics model integration tests, FT = $FT" begin
        default_params_filepath =
            joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
        toml_dict = LP.create_toml_dict(FT, default_params_filepath)
        domains = [
            Point(; z_sfc = FT(0.0)),
            Plane(;
                xlim = FT.((0, 1)),
                ylim = FT.((0, 1)),
                nelements = (2, 2),
                periodic = (true, true),
            ),
        ]

        AR_params = AutotrophicRespirationParameters(toml_dict)
        RTparams = BeerLambertParameters(FT)
        is_c3 = FT(1) # set the photosynthesis mechanism to C3
        photosynthesis_params = FarquharParameters(FT, is_c3)
        stomatal_g_params = MedlynConductanceParameters(toml_dict)

        AR_model = AutotrophicRespirationModel{FT}(AR_params)
        stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
        photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
        rt_model = BeerLambertModel{FT}(RTparams)

        earth_param_set = LP.LandParameters(FT)
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        LAI = (t) -> 1.0 # m2 [leaf] m-2 [ground]
        z_0m = FT(2.0) # m, Roughness length for momentum
        z_0b = FT(0.1) # m, Roughness length for scalars
        h_canopy = FT(20.0) # m, canopy height
        h_sfc = FT(20.0) # m, canopy height
        h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
        shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
            z_0m,
            z_0b,
            earth_param_set,
        )
        lat = FT(0.0) # degree
        long = FT(-180) # degree
        start_date = DateTime(2005)

        zenith_angle =
            (t, s) -> default_zenith_angle(
                t,
                s;
                insol_params = earth_param_set.insol_params,
                latitude = lat,
                longitude = long,
            )


        function shortwave_radiation(
            t;
            latitude = lat,
            longitude = long,
            insol_params = earth_param_set.insol_params,
        )
            return 1000 # W/m^2
        end

        function longwave_radiation(t)
            return 200 # W/m^2
        end

        u_atmos = t -> 10 #m.s-1

        liquid_precip = (t) -> 0 # m
        snow_precip = (t) -> 0 # m
        T_atmos = t -> 290 # Kelvin
        q_atmos = t -> 0.001 # kg/kg
        P_atmos = t -> 1e5 # Pa
        h_atmos = h_int # m
        c_atmos = (t) -> 4.11e-4 # mol/mol
        atmos = PrescribedAtmosphere(
            TimeVaryingInput(liquid_precip),
            TimeVaryingInput(snow_precip),
            TimeVaryingInput(T_atmos),
            TimeVaryingInput(u_atmos),
            TimeVaryingInput(q_atmos),
            TimeVaryingInput(P_atmos),
            start_date,
            h_atmos,
            earth_param_set;
            c_co2 = TimeVaryingInput(c_atmos),
        )
        radiation = PrescribedRadiativeFluxes(
            FT,
            TimeVaryingInput(shortwave_radiation),
            TimeVaryingInput(longwave_radiation),
            start_date;
            θs = zenith_angle,
            earth_param_set = earth_param_set,
        )
        Δz = FT(1.0) # height of compartments
        n_stem = Int64(5) # number of stem elements
        n_leaf = Int64(6) # number of leaf elements
        SAI = FT(1) # m2/m2
        RAI = FT(1) # m2/m2
        ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
            TimeVaryingInput(LAI),
            SAI,
            RAI,
        )
        K_sat_plant = 1.8e-8 # m/s.
        ψ63 = FT(-4 / 0.0098) # / MPa to m
        Weibull_param = FT(4) # unitless
        a = FT(0.05 * 0.0098) # 1/m
        plant_ν = FT(0.7) # m3/m3
        plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
        conductivity_model =
            PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
        retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
        compartment_midpoints = Vector{FT}(
            range(
                start = Δz / 2,
                step = Δz,
                stop = Δz * (n_stem + n_leaf) - (Δz / 2),
            ),
        )

        compartment_surfaces = Vector{FT}(
            range(start = 0.0, step = Δz, stop = Δz * (n_stem + n_leaf)),
        )


        function leaf_transpiration(t)
            T = FT(1e-8) # m/s
        end

        ψ_soil0 = FT(0.0)
        transpiration = PrescribedTranspiration{FT}(leaf_transpiration)

        soil_driver = PrescribedGroundConditions{FT}()

        autotrophic_parameters = AutotrophicRespirationParameters(toml_dict)
        autotrophic_respiration_model =
            AutotrophicRespirationModel{FT}(autotrophic_parameters)
        RD = FT(0.5)
        for domain in domains
            # check once with rooting_depth as float and once as field
            for rooting_depth in (RD, fill(RD, domain.space.surface))
                plant_hydraulics_param_set =
                    PlantHydraulics.PlantHydraulicsParameters(;
                        ai_parameterization = ai_parameterization,
                        ν = plant_ν,
                        S_s = plant_S_s,
                        rooting_depth = rooting_depth,
                        conductivity_model = conductivity_model,
                        retention_model = retention_model,
                    )
                plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
                    parameters = plant_hydraulics_param_set,
                    transpiration = transpiration,
                    n_stem = n_stem,
                    n_leaf = n_leaf,
                    compartment_surfaces = compartment_surfaces,
                    compartment_midpoints = compartment_midpoints,
                )
                model = ClimaLand.Canopy.CanopyModel{FT}(;
                    parameters = shared_params,
                    domain = domain,
                    autotrophic_respiration = autotrophic_respiration_model,
                    radiative_transfer = rt_model,
                    photosynthesis = photosynthesis_model,
                    conductance = stomatal_model,
                    hydraulics = plant_hydraulics,
                    boundary_conditions = Canopy.AtmosDrivenCanopyBC(
                        atmos,
                        radiation,
                        soil_driver,
                    ),
                )
                # Set system to hydrostatic equilibrium
                function initial_compute_exp_tendency!(F, Y)
                    AI = (; leaf = LAI(1.0), root = RAI, stem = SAI)
                    T0A = FT(1e-8) * AI[:leaf]
                    for i in 1:(n_leaf + n_stem)
                        if i == 1
                            fa =
                                water_flux.(
                                    -1 .* RD,
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
                                ) .* AI[:stem]
                        else
                            fa =
                                water_flux(
                                    plant_hydraulics.compartment_midpoints[i - 1],
                                    plant_hydraulics.compartment_midpoints[i],
                                    Y[i - 1],
                                    Y[i],
                                    PlantHydraulics.hydraulic_conductivity(
                                        conductivity_model,
                                        Y[i - 1],
                                    ),
                                    PlantHydraulics.hydraulic_conductivity(
                                        conductivity_model,
                                        Y[i],
                                    ),
                                ) * AI[plant_hydraulics.compartment_labels[i]]
                        end
                        F[i] = fa - T0A
                    end
                end
                #=======================
                Here we solve for the steady state of the hydraulics
                system. Then, the solution is used to check that evaluating the
                tendecy of the model also results in a steady state. This check is repeated using
                the plant hydraulics model directly.
                =======================#

                soln = nlsolve(
                    initial_compute_exp_tendency!,
                    Vector{FT}(-0.03:0.01:0.07);
                    ftol = eps(FT),
                    method = :newton,
                    iterations = 20,
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
                    Y.canopy.hydraulics.ϑ_l.:($i) .= ϑ_l_0[i]
                    p.canopy.hydraulics.ψ.:($i) .= NaN
                    p.canopy.hydraulics.fa.:($i) .= NaN
                    dY.canopy.hydraulics.ϑ_l.:($i) .= NaN
                end
                set_initial_cache! = make_set_initial_cache(model)
                set_initial_cache!(p, Y, 0.0)
                canopy_exp_tendency! = make_exp_tendency(model)
                canopy_exp_tendency!(dY, Y, p, 0.0)

                tolerance = sqrt(eps(FT))
                @test all(parent(dY.canopy.hydraulics.ϑ_l) .< tolerance) # starts in equilibrium


                # repeat using the plant hydraulics model directly
                # make sure it agrees with what we get when use the canopy model ODE
                Y, p, coords = initialize(model)
                standalone_dY = similar(Y)
                for i in 1:(n_stem + n_leaf)
                    Y.canopy.hydraulics.ϑ_l.:($i) .= ϑ_l_0[i]
                    p.canopy.hydraulics.ψ.:($i) .= NaN
                    p.canopy.hydraulics.fa.:($i) .= NaN
                    standalone_dY.canopy.hydraulics.ϑ_l.:($i) .= NaN
                end
                set_initial_cache!(p, Y, 0.0)
                standalone_exp_tendency! =
                    make_compute_exp_tendency(model.hydraulics, model)
                standalone_exp_tendency!(standalone_dY, Y, p, 0.0)
                @test all(
                    parent(
                        standalone_dY.canopy.hydraulics.ϑ_l .-
                        dY.canopy.hydraulics.ϑ_l,
                    ) .≈ FT(0),
                )
            end
        end
    end

    @testset "No plant, FT = $FT" begin
        domain = Point(; z_sfc = FT(0.0))

        default_params_filepath =
            joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
        toml_dict = LP.create_toml_dict(FT, default_params_filepath)
        AR_params = AutotrophicRespirationParameters(toml_dict)
        RTparams = BeerLambertParameters(FT)
        is_c3 = FT(1) # set the photosynthesis mechanism to C3
        photosynthesis_params = FarquharParameters(FT, is_c3)
        stomatal_g_params = MedlynConductanceParameters(toml_dict)

        AR_model = AutotrophicRespirationModel{FT}(AR_params)
        stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
        photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
        rt_model = BeerLambertModel{FT}(RTparams)

        earth_param_set = LP.LandParameters(FT)
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        LAI = FT(0.0) # m2 [leaf] m-2 [ground]
        z_0m = FT(2.0) # m, Roughness length for momentum
        z_0b = FT(0.1) # m, Roughness length for scalars
        h_canopy = FT(0.0) # m, canopy height
        h_sfc = FT(0.0) # m
        h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
        shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
            z_0m,
            z_0b,
            earth_param_set,
        )
        lat = FT(0.0) # degree
        long = FT(-180) # degree
        start_date = DateTime(2005)

        zenith_angle =
            (t, s) -> default_zenith_angle(
                t,
                s;
                insol_params = earth_param_set.insol_params,
                latitude = lat,
                longitude = long,
            )

        function shortwave_radiation(
            t;
            latitude = lat,
            longitude = long,
            insol_params = earth_param_set.insol_params,
        )
            return 1000 # W/m^2
        end

        function longwave_radiation(t)
            return 200 # W/m^2
        end

        u_atmos = t -> 10 #m.s-1

        liquid_precip = (t) -> 0 # m
        snow_precip = (t) -> 0 # m
        T_atmos = t -> 290 # Kelvin
        q_atmos = t -> 0.001 # kg/kg
        P_atmos = t -> 1e5 # Pa
        h_atmos = h_int # m
        c_atmos = (t) -> 4.11e-4 # mol/mol
        atmos = PrescribedAtmosphere(
            TimeVaryingInput(liquid_precip),
            TimeVaryingInput(snow_precip),
            TimeVaryingInput(T_atmos),
            TimeVaryingInput(u_atmos),
            TimeVaryingInput(q_atmos),
            TimeVaryingInput(P_atmos),
            start_date,
            h_atmos,
            earth_param_set;
            c_co2 = TimeVaryingInput(c_atmos),
        )
        radiation = PrescribedRadiativeFluxes(
            FT,
            TimeVaryingInput(shortwave_radiation),
            TimeVaryingInput(longwave_radiation),
            start_date;
            θs = zenith_angle,
            earth_param_set = earth_param_set,
        )

        n_stem = Int64(0) # number of stem elements
        n_leaf = Int64(1) # number of leaf elements
        SAI = FT(0) # m2/m2
        RAI = FT(0) # m2/m2
        lai_fun = t -> 0
        ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
            TimeVaryingInput(lai_fun),
            SAI,
            RAI,
        )
        K_sat_plant = 0 # m/s.
        ψ63 = FT(-4 / 0.0098) # / MPa to m
        Weibull_param = FT(4) # unitless
        a = FT(0.05 * 0.0098) # 1/m
        plant_ν = FT(0.1) # m3/m3
        plant_S_s = FT(1e-2)
        conductivity_model =
            PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
        retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
        compartment_midpoints = [h_canopy]
        compartment_surfaces = [FT(0.0), h_canopy]
        # set rooting_depth param to largest possible value to test no roots
        param_set = PlantHydraulics.PlantHydraulicsParameters(;
            ai_parameterization = ai_parameterization,
            ν = plant_ν,
            S_s = plant_S_s,
            rooting_depth = maxintfloat(FT),
            conductivity_model = conductivity_model,
            retention_model = retention_model,
        )

        transpiration = DiagnosticTranspiration{FT}()
        soil_driver = PrescribedGroundConditions{FT}()
        plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
            parameters = param_set,
            transpiration = transpiration,
            n_stem = n_stem,
            n_leaf = n_leaf,
            compartment_surfaces = compartment_surfaces,
            compartment_midpoints = compartment_midpoints,
        )

        autotrophic_parameters = AutotrophicRespirationParameters(toml_dict)
        autotrophic_respiration_model =
            AutotrophicRespirationModel(autotrophic_parameters)

        model = ClimaLand.Canopy.CanopyModel{FT}(;
            parameters = shared_params,
            domain = domain,
            autotrophic_respiration = autotrophic_respiration_model,
            radiative_transfer = rt_model,
            photosynthesis = photosynthesis_model,
            conductance = stomatal_model,
            hydraulics = plant_hydraulics,
            boundary_conditions = Canopy.AtmosDrivenCanopyBC(
                atmos,
                radiation,
                soil_driver,
            ),
        )

        Y, p, coords = initialize(model)
        dY = similar(Y)
        for i in 1:(n_stem + n_leaf)
            Y.canopy.hydraulics.ϑ_l.:($i) .= FT(0.1)
            p.canopy.hydraulics.ψ.:($i) .= NaN
            p.canopy.hydraulics.fa.:($i) .= NaN
            dY.canopy.hydraulics.ϑ_l.:($i) .= NaN
        end
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, FT(0.0))
        @test all(parent(p.canopy.hydraulics.fa) .≈ FT(0.0))
        @test all(parent(p.canopy.hydraulics.fa_roots) .≈ FT(0.0))
        @test all(parent(p.canopy.turbulent_fluxes.transpiration) .≈ FT(0.0))
        @test all(parent(p.canopy.radiative_transfer.par.abs) .≈ FT(0.0))
        exp_tend! = make_exp_tendency(model)
        exp_tend!(dY, Y, p, FT(0))
        @test all(parent(dY.canopy.hydraulics.ϑ_l.:1) .≈ FT(0.0))
    end
end

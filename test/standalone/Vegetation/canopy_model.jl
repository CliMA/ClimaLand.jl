using Test
import ClimaParams
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Thermodynamics
using Dates
using StaticArrays
using ClimaLand
using ClimaLand: PrescribedAtmosphere, PrescribedRadiativeFluxes
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Domains: Point
import Insolation
using ClimaCore.MatrixFields: @name
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams


@testset "Canopy software pipes" begin
    for FT in (Float32, Float64)
        domain = ClimaLand.Domains.SphericalSurface(;
            radius = FT(100.0),
            nelements = 10,
            npolynomial = 1,
        )
        # create a field with both 1.0s and 0.0s
        mechanism_field = ClimaCore.Fields.Field(FT, domain.space.surface)
        ClimaCore.Fields.set!(
            x -> x.coordinates.lat > 0 ? 0.0 : 1.0,
            mechanism_field,
        )
        # create one case where parameters are spatially varying and one where not
        g1_cases = (FT(790), fill(FT(790), domain.space.surface))
        Vcmax25_cases = (FT(9e-5), fill(FT(9e-5), domain.space.surface))
        mechanism_cases = (FT(1), mechanism_field)
        rooting_cases = (FT(0.5), fill(FT(0.5), domain.space.surface))
        # test default values as field
        α_PAR_leaf_cases = (FT(0.1), fill(FT(0.1), domain.space.surface))
        α_NIR_leaf_cases = (FT(0.4), fill(FT(0.4), domain.space.surface))
        ld_cases = (FT(0.5), fill(FT(0.5), domain.space.surface))
        zipped_params = zip(
            g1_cases,
            Vcmax25_cases,
            mechanism_cases,
            rooting_cases,
            α_PAR_leaf_cases,
            α_NIR_leaf_cases,
            ld_cases,
        )
        for (g1, Vcmax25, is_c3, rooting_depth, α_PAR_leaf, α_NIR_leaf, ld) in
            zipped_params
            AR_params = AutotrophicRespirationParameters(FT)
            G_Function = ConstantGFunction(ld)
            RTparams =
                BeerLambertParameters(FT; α_PAR_leaf, α_NIR_leaf, G_Function)
            photosynthesis_params = FarquharParameters(FT, is_c3; Vcmax25)
            stomatal_g_params = MedlynConductanceParameters(FT; g1)
            AR_model = AutotrophicRespirationModel{FT}(AR_params)
            stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
            photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
            rt_model = BeerLambertModel{FT}(RTparams)

            earth_param_set = LP.LandParameters(FT)
            LAI = FT(8.0) # m2 [leaf] m-2 [ground]
            z_0m = FT(2.0) # m, Roughness length for momentum - value from tall forest ChatGPT
            z_0b = FT(0.1) # m, Roughness length for scalars - value from tall forest ChatGPT
            h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
            shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
                z_0m,
                z_0b,
                earth_param_set,
            )
            lat = FT(0.0) # degree
            long = FT(-180) # degree

            function zenith_angle(
                t,
                start_date;
                latitude = lat,
                longitude = long,
                insol_params = earth_param_set.insol_params,
            )
                current_datetime = start_date + Dates.Second(round(t))
                d, δ, η_UTC =
                    FT.(
                        Insolation.helper_instantaneous_zenith_angle(
                            current_datetime,
                            start_date,
                            insol_params,
                        )
                    )
                return Insolation.instantaneous_zenith_angle(
                    d,
                    δ,
                    η_UTC,
                    longitude,
                    latitude,
                )[1]
            end

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
            start_date = DateTime(2005)
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
            )

            # Plant Hydraulics
            RAI = FT(1)
            SAI = FT(0)
            lai_fun = t -> LAI
            ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
                TimeVaryingInput(lai_fun),
                SAI,
                RAI,
            )
            K_sat_plant = FT(1.8e-8) # m/s
            ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
            Weibull_param = FT(4) # unitless, Holtzman's original c param value
            a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
            plant_ν = FT(0.7) # m3/m3
            plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
            conductivity_model =
                PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
            retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
            root_depths =
                SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0) # 1st element is the deepest root depth
            param_set = PlantHydraulics.PlantHydraulicsParameters(;
                ai_parameterization = ai_parameterization,
                ν = plant_ν,
                S_s = plant_S_s,
                rooting_depth = rooting_depth,
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
                        range(
                            start = 0.0,
                            step = Δz,
                            stop = Δz * (n_stem + n_leaf),
                        ),
                    )
                )

            ψ_soil0 = FT(0.0)
            soil_driver = PrescribedGroundConditions(FT)
            plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
                parameters = param_set,
                n_stem = n_stem,
                n_leaf = n_leaf,
                compartment_surfaces = compartment_faces,
                compartment_midpoints = compartment_centers,
            )
            canopy = ClimaLand.Canopy.CanopyModel{FT}(;
                parameters = shared_params,
                domain = domain,
                autotrophic_respiration = AR_model,
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
            drivers = ClimaLand.get_drivers(canopy)
            @test drivers == (atmos, radiation)
            Y, p, coords = ClimaLand.initialize(canopy)
            @test propertynames(p.drivers) == (
                :P_liq,
                :P_snow,
                :T,
                :P,
                :u,
                :q,
                :c_co2,
                :thermal_state,
                :SW_d,
                :LW_d,
                :θs,
            )
            # Check that structure of Y is value (will error if not)
            @test !isnothing(zero(Y))
            @test typeof(canopy.energy) == PrescribedCanopyTempModel{FT}
            @test propertynames(p) == (:canopy, :dss_buffer_2d, :drivers)
            for component in ClimaLand.Canopy.canopy_components(canopy)
                # Only hydraulics has a prognostic variable
                if component == :hydraulics
                    @test propertynames(getproperty(Y.canopy, component)) ==
                          ClimaLand.prognostic_vars(
                        getproperty(canopy, component),
                    )
                end
                @test propertynames(getproperty(p.canopy, component)) ==
                      ClimaLand.auxiliary_vars(getproperty(canopy, component))
                @test getproperty(auxiliary_types(canopy), component) ==
                      auxiliary_types(getproperty(canopy, component))
                @test getproperty(auxiliary_vars(canopy), component) ==
                      auxiliary_vars(getproperty(canopy, component))
                @test getproperty(prognostic_types(canopy), component) ==
                      prognostic_types(getproperty(canopy, component))
                @test getproperty(prognostic_types(canopy), component) ==
                      prognostic_types(getproperty(canopy, component))
            end
            parent(Y.canopy.hydraulics.ϑ_l) .= plant_ν
            set_initial_cache! = make_set_initial_cache(canopy)
            exp_tendency! = make_exp_tendency(canopy)
            imp_tendency! = ClimaLand.make_imp_tendency(canopy)
            jacobian! = ClimaLand.make_jacobian(canopy)
            # set up jacobian info
            jac_kwargs = (;
                jac_prototype = ImplicitEquationJacobian(Y),
                Wfact = jacobian!,
            )

            t0 = FT(0.0)
            dY = similar(Y)
            set_initial_cache!(p, Y, t0)
            # check that this is updated correctly:
            # @test p.canopy.autotrophic_respiration.Ra ==
            exp_tendency!(dY, Y, p, t0)
            turb_fluxes = ClimaLand.turbulent_fluxes(atmos, canopy, Y, p, t0)

            @test p.canopy.hydraulics.fa.:1 == turb_fluxes.transpiration
            @test p.canopy.energy.turbulent_fluxes.shf == turb_fluxes.shf
            @test p.canopy.energy.turbulent_fluxes.lhf == turb_fluxes.lhf
            @test p.canopy.energy.turbulent_fluxes.transpiration ==
                  turb_fluxes.transpiration
            _σ = FT(LP.Stefan(earth_param_set))
            f_abs_par = p.canopy.radiative_transfer.par.abs
            f_abs_nir = p.canopy.radiative_transfer.nir.abs
            nir_d = p.canopy.radiative_transfer.nir_d
            par_d = p.canopy.radiative_transfer.par_d
            @test p.canopy.radiative_transfer.SW_n ==
                  @. f_abs_par * par_d + f_abs_nir * nir_d
            ϵ_canopy = p.canopy.radiative_transfer.ϵ
            T_canopy = FT.(T_atmos(t0))
            T_soil = FT.(soil_driver.T(t0))
            ϵ_soil = FT.(soil_driver.ϵ)
            LW_d = FT.(longwave_radiation(t0))
            LW_d_canopy = @. (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
            LW_u_soil = @. ϵ_soil * _σ * T_soil^4 + (1 - ϵ_soil) * LW_d_canopy
            @test p.canopy.radiative_transfer.LW_n == @. (
                ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 +
                ϵ_canopy * LW_u_soil
            )


            # Penman-monteith
            Δ = FT(100 * (0.444017302 + (290 - 273.15) * 0.0286064092))
            Rn = FT(shortwave_radiation(t0))
            G = FT(0.0)
            thermo_params = canopy.parameters.earth_param_set.thermo_params
            ts_in =
                Thermodynamics.PhaseEquil_pTq.(
                    thermo_params,
                    p.drivers.P,
                    p.drivers.T,
                    p.drivers.q,
                )
            ρa = Thermodynamics.air_density.(thermo_params, ts_in)
            cp = FT(
                cp_m(
                    thermo_params,
                    Thermodynamics.PhasePartition.(q_atmos(t0)),
                ),
            )

            es =
                Thermodynamics.saturation_vapor_pressure.(
                    Ref(thermo_params),
                    FT.(T_atmos(t0)),
                    Ref(Thermodynamics.Liquid()),
                )
            ea =
                Thermodynamics.partial_pressure_vapor.(
                    thermo_params,
                    FT(P_atmos(t0)),
                    Thermodynamics.PhasePartition.(FT.(q_atmos(t0))),
                )

            VPD = es .- ea

            conditions = turbulent_fluxes(atmos, canopy, Y, p, t0) #Per unit m^2 of leaf
            r_ae = Array(parent(conditions.r_ae))[1] # s/m
            ga = 1 / r_ae
            γ = FT(66)
            R = FT(LP.gas_constant(earth_param_set))
            gs = Array(
                parent(
                    ClimaLand.Canopy.upscale_leaf_conductance.(
                        p.canopy.conductance.gs,
                        LAI,
                        FT.(T_atmos(t0)),
                        R,
                        FT.(P_atmos(t0)),
                    ),
                ),
            )[1]
            Lv = FT(2453e6) #J/m^3

            ET = penman_monteith.(
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
                (Array(parent(turb_fluxes.transpiration .- ET))[1]) /
                Array(parent(turb_fluxes.transpiration))[1],
            ) < 0.5

            @test ClimaLand.surface_evaporative_scaling(canopy, Y, p) == FT(1.0)
            @test ClimaLand.surface_height(canopy, Y, p) == compartment_faces[1]
            T_sfc = FT.(T_atmos(t0))
            @test all(
                Array(
                    parent(ClimaLand.surface_temperature(canopy, Y, p, t0)),
                ) .== [T_sfc],
            )
            @test ClimaLand.surface_temperature(canopy, Y, p, t0) isa
                  ClimaCore.Fields.Field
            @test all(
                Array(
                    parent(
                        ClimaLand.Canopy.canopy_temperature(
                            canopy.energy,
                            canopy,
                            Y,
                            p,
                            t0,
                        ),
                    ),
                ) .== [T_sfc],
            )
            @test ClimaLand.Canopy.canopy_temperature(
                canopy.energy,
                canopy,
                Y,
                p,
                t0,
            ) isa ClimaCore.Fields.Field
        end
    end
end


@testset "Component prescribed fields" begin
    struct Default{FT} <: ClimaLand.Canopy.AbstractCanopyComponent{FT} end
    for FT in (Float32, Float64)
        # Plant Hydraulics
        LAI = FT(2)
        RAI = FT(1)
        SAI = FT(1)
        lai_fun = t -> LAI * sin(t * 2π / 365)
        ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
            TimeVaryingInput(lai_fun),
            SAI,
            RAI,
        )
        K_sat_plant = FT(1.8e-8) # m/s
        ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
        Weibull_param = FT(4) # unitless, Holtzman's original c param value
        a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
        plant_ν = FT(0.7) # m3/m3
        plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
        conductivity_model =
            PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
        retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
        root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0) # 1st element is the deepest root depth
        param_set = PlantHydraulics.PlantHydraulicsParameters(;
            ai_parameterization = ai_parameterization,
            ν = plant_ν,
            S_s = plant_S_s,
            rooting_depth = FT(0.5),
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
                    range(
                        start = 0.0,
                        step = Δz,
                        stop = Δz * (n_stem + n_leaf),
                    ),
                )
            )

        plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
            parameters = param_set,
            n_stem = n_stem,
            n_leaf = n_leaf,
            compartment_surfaces = compartment_faces,
            compartment_midpoints = compartment_centers,
        )

        t0 = FT(100)
        domain = Point(; z_sfc = FT(0.0))
        p = ClimaCore.fill(
            (;
                canopy = (;
                    hydraulics = (;
                        area_index = (
                            leaf = FT(0.0),
                            root = FT(0.0),
                            stem = FT(0.0),
                        )
                    )
                )
            ),
            domain.space.surface,
        )
        # Test that they are set properly
        set_canopy_prescribed_field!(plant_hydraulics, p, t0)
        @test all(
            Array(parent(p.canopy.hydraulics.area_index.leaf)) .==
            FT(LAI * sin(t0 * 2π / 365)),
        )
        @test all(
            Array(parent(p.canopy.hydraulics.area_index.stem)) .== FT(1.0),
        )
        @test all(
            Array(parent(p.canopy.hydraulics.area_index.root)) .== FT(1.0),
        )

        # Test that LAI is updated
        set_canopy_prescribed_field!(plant_hydraulics, p, FT(200))
        @test all(
            Array(parent(p.canopy.hydraulics.area_index.leaf)) .==
            ClimaLand.Canopy.PlantHydraulics.clip(
                FT(LAI * sin(200 * 2π / 365)),
                FT(0.05),
            ),
        )

        set_canopy_prescribed_field!(Default{FT}(), p, t0)
        set_canopy_prescribed_field!(Default{FT}(), p, t0)
        # Test that they are unchanged
        @test all(
            parent(p.canopy.hydraulics.area_index.leaf) .==
            ClimaLand.Canopy.PlantHydraulics.clip(
                FT(LAI * sin(200 * 2π / 365)),
                FT(0.05),
            ),
        )
        @test all(
            Array(parent(p.canopy.hydraulics.area_index.stem)) .== FT(1.0),
        )
        @test all(
            Array(parent(p.canopy.hydraulics.area_index.root)) .== FT(1.0),
        )
    end
end

@testset "PrescribedGroundConditions" begin
    for FT in (Float32, Float64)
        soil_driver = PrescribedGroundConditions(FT)
        @test ground_albedo_PAR(
            Val((:canopy,)),
            soil_driver,
            nothing,
            nothing,
            nothing,
        ) == FT(0.2)
        @test ground_albedo_NIR(
            Val((:canopy,)),
            soil_driver,
            nothing,
            nothing,
            nothing,
        ) == FT(0.4)
        @test soil_driver.root_depths ==
              SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)
        @test FT.(soil_driver.ψ(2.0)) == FT.(0.0)
        @test FT.(soil_driver.T(2.0)) == FT.(298.0)
    end
end

@testset "Canopy software pipes with energy model" begin
    for FT in (Float32, Float64)
        domain = ClimaLand.Domains.SphericalSurface(;
            radius = FT(100.0),
            nelements = 10,
            npolynomial = 1,
        )
        # create a field with both 1.0s and 0.0s
        mechanism_field = ClimaCore.Fields.Field(FT, domain.space.surface)
        ClimaCore.Fields.set!(
            x -> x.coordinates.lat > 0 ? 0.0 : 1.0,
            mechanism_field,
        )
        # create one case where parameters are spatially varying and one where not
        g1_cases = (FT(790), fill(FT(790), domain.space.surface))
        Vcmax25_cases = (FT(9e-5), fill(FT(9e-5), domain.space.surface))
        mechanism_cases = (FT(1), mechanism_field)
        rooting_cases = (FT(0.5), fill(FT(0.5), domain.space.surface))
        # test default values as field
        α_PAR_leaf_cases = (FT(0.1), fill(FT(0.1), domain.space.surface))
        α_NIR_leaf_cases = (FT(0.4), fill(FT(0.4), domain.space.surface))
        ld_cases = (FT(0.5), fill(FT(0.5), domain.space.surface))
        zipped_params = zip(
            g1_cases,
            Vcmax25_cases,
            mechanism_cases,
            rooting_cases,
            α_PAR_leaf_cases,
            α_NIR_leaf_cases,
            ld_cases,
        )
        for (g1, Vcmax25, is_c3, rooting_depth, α_PAR_leaf, α_NIR_leaf, ld) in
            zipped_params
            G_Function = ConstantGFunction(ld)
            RTparams =
                BeerLambertParameters(FT; α_PAR_leaf, α_NIR_leaf, G_Function)
            photosynthesis_params = FarquharParameters(FT, is_c3; Vcmax25)
            stomatal_g_params = MedlynConductanceParameters(FT; g1)

            stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
            photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
            rt_model = BeerLambertModel{FT}(RTparams)
            energy_model = BigLeafEnergyModel{FT}(BigLeafEnergyParameters{FT}())
            earth_param_set = LP.LandParameters(FT)
            LAI = FT(8.0) # m2 [leaf] m-2 [ground]
            z_0m = FT(2.0) # m, Roughness length for momentum - value from tall forest ChatGPT
            z_0b = FT(0.1) # m, Roughness length for scalars - value from tall forest ChatGPT
            h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
            shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
                z_0m,
                z_0b,
                earth_param_set,
            )
            lat = FT(0.0) # degree
            long = FT(-180) # degree

            function zenith_angle(
                t,
                start_date;
                latitude = lat,
                longitude = long,
                insol_params = earth_param_set.insol_params,
            )
                current_datetime = start_date + Dates.Second(round(t))
                d, δ, η_UTC =
                    FT.(
                        Insolation.helper_instantaneous_zenith_angle(
                            current_datetime,
                            start_date,
                            insol_params,
                        )
                    )
                return Insolation.instantaneous_zenith_angle(
                    d,
                    δ,
                    η_UTC,
                    longitude,
                    latitude,
                )[1]
            end

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
            start_date = DateTime(2005)
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
            )


            # Plant Hydraulics
            RAI = FT(1)
            SAI = FT(0)
            lai_fun = t -> LAI
            ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
                TimeVaryingInput(lai_fun),
                SAI,
                RAI,
            )
            K_sat_plant = FT(1.8e-8) # m/s
            ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
            Weibull_param = FT(4) # unitless, Holtzman's original c param value
            a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
            plant_ν = FT(0.7) # m3/m3
            plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
            conductivity_model =
                PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
            retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
            root_depths =
                SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0) # 1st element is the deepest root depth
            param_set = PlantHydraulics.PlantHydraulicsParameters(;
                ai_parameterization = ai_parameterization,
                ν = plant_ν,
                S_s = plant_S_s,
                rooting_depth = rooting_depth,
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
                        range(
                            start = 0.0,
                            step = Δz,
                            stop = Δz * (n_stem + n_leaf),
                        ),
                    )
                )

            ψ_soil0 = FT(0.0)
            T_soil0 = FT(290)
            soil_driver = PrescribedGroundConditions(
                root_depths,
                (t) -> ψ_soil0,
                (t) -> T_soil0,
                FT(0.2),
                FT(0.4),
                FT(0.98),
            )

            plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
                parameters = param_set,
                n_stem = n_stem,
                n_leaf = n_leaf,
                compartment_surfaces = compartment_faces,
                compartment_midpoints = compartment_centers,
            )
            autotrophic_parameters = AutotrophicRespirationParameters(FT)
            autotrophic_respiration_model =
                AutotrophicRespirationModel{FT}(autotrophic_parameters)

            canopy = ClimaLand.Canopy.CanopyModel{FT}(;
                parameters = shared_params,
                domain = domain,
                radiative_transfer = rt_model,
                photosynthesis = photosynthesis_model,
                conductance = stomatal_model,
                autotrophic_respiration = autotrophic_respiration_model,
                energy = energy_model,
                hydraulics = plant_hydraulics,
                boundary_conditions = Canopy.AtmosDrivenCanopyBC(
                    atmos,
                    radiation,
                    soil_driver,
                ),
            )
            # This test needs to match ClimaParams `canopy_emissivity`
            @test canopy.radiative_transfer.parameters.ϵ_canopy == FT(0.97)
            @test canopy.energy.parameters.ac_canopy == FT(2.0e3)
            Y, p, coords = ClimaLand.initialize(canopy)

            # Check that structure of Y is value (will error if not)
            @test !isnothing(zero(Y))
            @test propertynames(p) == (:canopy, :dss_buffer_2d, :drivers)
            for component in ClimaLand.Canopy.canopy_components(canopy)
                # Only hydraulics has a prognostic variable
                if component == :hydraulics
                    @test propertynames(getproperty(Y.canopy, component)) ==
                          ClimaLand.prognostic_vars(
                        getproperty(canopy, component),
                    )
                end
                @test propertynames(getproperty(p.canopy, component)) ==
                      ClimaLand.auxiliary_vars(getproperty(canopy, component))

                @test getproperty(auxiliary_types(canopy), component) ==
                      auxiliary_types(getproperty(canopy, component))
                @test getproperty(auxiliary_vars(canopy), component) ==
                      auxiliary_vars(getproperty(canopy, component))
                @test getproperty(prognostic_types(canopy), component) ==
                      prognostic_types(getproperty(canopy, component))
                @test getproperty(prognostic_types(canopy), component) ==
                      prognostic_types(getproperty(canopy, component))
            end
            Y.canopy.hydraulics .= plant_ν
            Y.canopy.energy.T = FT(289)

            set_initial_cache! = make_set_initial_cache(canopy)
            t0 = FT(0.0)
            set_initial_cache!(p, Y, t0)
            exp_tendency! = make_exp_tendency(canopy)
            imp_tendency! = ClimaLand.make_imp_tendency(canopy)
            jacobian! = ClimaLand.make_jacobian(canopy)
            # set up jacobian info
            jac_kwargs = (;
                jac_prototype = ImplicitEquationJacobian(Y),
                Wfact = jacobian!,
            )

            dY = similar(Y)
            exp_tendency!(dY, Y, p, t0)
            turb_fluxes = ClimaLand.turbulent_fluxes(atmos, canopy, Y, p, t0)

            @test p.canopy.hydraulics.fa.:1 == turb_fluxes.transpiration
            @test p.canopy.energy.turbulent_fluxes.lhf == turb_fluxes.lhf
            @test p.canopy.energy.turbulent_fluxes.shf == turb_fluxes.shf
            @test all(Array(parent(p.canopy.energy.fa_energy_roots)) .== FT(0))

            @test all(
                Array(
                    parent(ClimaLand.surface_temperature(canopy, Y, p, t0)),
                ) .== FT(289),
            )
            @test all(
                Array(
                    parent(
                        ClimaLand.Canopy.canopy_temperature(
                            canopy.energy,
                            canopy,
                            Y,
                            p,
                            t0,
                        ),
                    ),
                ) .== FT(289),
            )
        end
    end
end

@testset "Jacobian for Temperature" begin
    for FT in (Float32, Float64)
        domain = Point(; z_sfc = FT(0.0))

        g1 = FT(790)
        Vcmax25 = FT(9e-5)
        is_c3 = FT(1)
        RTparams = BeerLambertParameters(FT)
        photosynthesis_params = FarquharParameters(FT, is_c3; Vcmax25)
        stomatal_g_params = MedlynConductanceParameters(FT; g1)

        stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
        photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
        rt_model = BeerLambertModel{FT}(RTparams)
        energy_model = BigLeafEnergyModel{FT}(BigLeafEnergyParameters{FT}())
        earth_param_set = LP.LandParameters(FT)
        LAI = FT(8.0) # m2 [leaf] m-2 [ground]
        z_0m = FT(2.0) # m, Roughness length for momentum - value from tall forest ChatGPT
        z_0b = FT(0.1) # m, Roughness length for scalars - value from tall forest ChatGPT
        h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
        shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
            z_0m,
            z_0b,
            earth_param_set,
        )
        lat = FT(0.0) # degree
        long = FT(-180) # degree

        function zenith_angle(
            t,
            start_date;
            latitude = lat,
            longitude = long,
            insol_params = earth_param_set.insol_params,
        )
            current_datetime = start_date + Dates.Second(round(t))
            d, δ, η_UTC =
                FT.(
                    Insolation.helper_instantaneous_zenith_angle(
                        current_datetime,
                        start_date,
                        insol_params,
                    )
                )
            return Insolation.instantaneous_zenith_angle(
                d,
                δ,
                η_UTC,
                longitude,
                latitude,
            )[1]
        end

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
        start_date = DateTime(2005)
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
        )

        # Plant Hydraulics
        RAI = FT(1)
        SAI = FT(0)
        lai_fun = t -> LAI
        ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
            TimeVaryingInput(lai_fun),
            SAI,
            RAI,
        )
        K_sat_plant = FT(1.8e-8) # m/s
        ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
        Weibull_param = FT(4) # unitless, Holtzman's original c param value
        a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
        plant_ν = FT(0.7) # m3/m3
        plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
        conductivity_model =
            PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
        retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
        root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0) # 1st element is the deepest root depth
        rooting_depth = FT(0.5)
        param_set = PlantHydraulics.PlantHydraulicsParameters(;
            ai_parameterization = ai_parameterization,
            ν = plant_ν,
            S_s = plant_S_s,
            rooting_depth = rooting_depth,
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
                    range(
                        start = 0.0,
                        step = Δz,
                        stop = Δz * (n_stem + n_leaf),
                    ),
                )
            )

        ψ_soil0 = FT(0.0)
        T_soil0 = FT(290)
        soil_driver = PrescribedGroundConditions(
            root_depths,
            (t) -> ψ_soil0,
            (t) -> T_soil0,
            FT(0.2),
            FT(0.4),
            FT(0.98),
        )

        plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
            parameters = param_set,
            n_stem = n_stem,
            n_leaf = n_leaf,
            compartment_surfaces = compartment_faces,
            compartment_midpoints = compartment_centers,
        )
        autotrophic_parameters = AutotrophicRespirationParameters(FT)
        autotrophic_respiration_model =
            AutotrophicRespirationModel{FT}(autotrophic_parameters)

        canopy = ClimaLand.Canopy.CanopyModel{FT}(;
            parameters = shared_params,
            domain = domain,
            radiative_transfer = rt_model,
            photosynthesis = photosynthesis_model,
            conductance = stomatal_model,
            autotrophic_respiration = autotrophic_respiration_model,
            energy = energy_model,
            hydraulics = plant_hydraulics,
            boundary_conditions = Canopy.AtmosDrivenCanopyBC(
                atmos,
                radiation,
                soil_driver,
            ),
        )

        Y, p, coords = ClimaLand.initialize(canopy)

        Y.canopy.hydraulics .= plant_ν
        Y.canopy.energy.T = FT(289)

        set_initial_cache! = make_set_initial_cache(canopy)
        t0 = FT(0.0)
        imp_tendency! = ClimaLand.make_imp_tendency(canopy)
        jacobian! = ClimaLand.make_jacobian(canopy)

        set_initial_cache!(p, Y, t0)
        T_sfc = ClimaLand.surface_temperature(canopy, Y, p, t0)
        ρ_sfc = ClimaLand.surface_air_density(
            canopy.boundary_conditions.atmos,
            canopy,
            Y,
            p,
            t0,
            T_sfc,
        )
        thermo_params = canopy.parameters.earth_param_set.thermo_params
        q_sfc =
            Thermodynamics.q_vap_saturation_generic.(
                thermo_params,
                T_sfc,
                ρ_sfc,
                Thermodynamics.Liquid(),
            )
        dY = similar(Y)
        imp_tendency!(dY, Y, p, t0)
        jac = ImplicitEquationJacobian(Y)
        jacobian!(jac, Y, p, FT(1), t0)
        jac_value =
            parent(jac.matrix[@name(canopy.energy.T), @name(canopy.energy.T)])
        ΔT = FT(0.01)

        Y_2 = deepcopy(Y)
        Y_2.canopy.energy.T = FT(289 + ΔT)
        p_2 = deepcopy(p)
        set_initial_cache!(p_2, Y_2, t0)
        T_sfc2 = ClimaLand.surface_temperature(canopy, Y_2, p_2, t0)
        ρ_sfc2 = ClimaLand.surface_air_density(
            canopy.boundary_conditions.atmos,
            canopy,
            Y_2,
            p_2,
            t0,
            T_sfc2,
        )
        q_sfc2 =
            Thermodynamics.q_vap_saturation_generic.(
                thermo_params,
                T_sfc2,
                ρ_sfc2,
                Thermodynamics.Liquid(),
            )
        dY_2 = similar(Y_2)
        imp_tendency!(dY_2, Y_2, p_2, t0)

        finitediff_LW =
            (
                p_2.canopy.radiative_transfer.LW_n .-
                p.canopy.radiative_transfer.LW_n
            ) ./ ΔT
        estimated_LW = p.canopy.energy.∂LW_n∂Tc
        @test parent(abs.(finitediff_LW .- estimated_LW) ./ finitediff_LW)[1] <
              0.01

        finitediff_SHF =
            (
                p_2.canopy.energy.turbulent_fluxes.shf .-
                p.canopy.energy.turbulent_fluxes.shf
            ) ./ ΔT
        estimated_SHF = p.canopy.energy.turbulent_fluxes.∂SHF∂Tc
        @test parent(abs.(finitediff_SHF .- estimated_SHF) ./ finitediff_SHF)[1] <
              0.15

        finitediff_LHF =
            (
                p_2.canopy.energy.turbulent_fluxes.lhf .-
                p.canopy.energy.turbulent_fluxes.lhf
            ) ./ ΔT
        estimated_LHF =
            p.canopy.energy.turbulent_fluxes.∂LHF∂qc .* p.canopy.energy.∂qc∂Tc
        @test parent(abs.(finitediff_LHF .- estimated_LHF) ./ finitediff_LHF)[1] <
              0.3

        # Error in `q` derivative is large
        finitediff_q = (q_sfc2 .- q_sfc) ./ ΔT
        estimated_q = p.canopy.energy.∂qc∂Tc
        @test parent(abs.(finitediff_q .- estimated_q) ./ finitediff_q)[1] <
              0.25

        # Im not sure why this is not smaller! There must be an error in ∂LHF∂qc also.
        estimated_LHF_with_correct_q =
            p.canopy.energy.turbulent_fluxes.∂LHF∂qc .* finitediff_q
        @test parent(
            abs.(finitediff_LHF .- estimated_LHF_with_correct_q) ./
            finitediff_LHF,
        )[1] < 0.5

        # Recall jac = ∂Ṫ∂T - 1 [dtγ = 1]
        ∂Ṫ∂T = jac_value .+ 1
        @test (abs.(
            parent((dY_2.canopy.energy.T .- dY.canopy.energy.T) ./ ΔT) - ∂Ṫ∂T
        ) / ∂Ṫ∂T)[1] < 0.25 # Error propagates here from ∂LHF∂T
    end
end


@testset "Zero LAI;" begin
    for FT in (Float32, Float64)
        domain = ClimaLand.Domains.SphericalSurface(;
            radius = FT(100.0),
            nelements = 10,
            npolynomial = 1,
        )
        # create a field with both 1.0s and 0.0s
        mechanism_field = ClimaCore.Fields.Field(FT, domain.space.surface)
        ClimaCore.Fields.set!(
            x -> x.coordinates.lat > 0 ? 0.0 : 1.0,
            mechanism_field,
        )
        # create one case where parameters are spatially varying and one where not
        g1_cases = (FT(790), fill(FT(790), domain.space.surface))
        Vcmax25_cases = (FT(9e-5), fill(FT(9e-5), domain.space.surface))
        mechanism_cases = (FT(1), mechanism_field)
        rooting_cases = (FT(0.5), fill(FT(0.5), domain.space.surface))
        α_PAR_leaf_cases = (FT(0.1), fill(FT(0.1), domain.space.surface))
        τ_PAR_leaf_cases = (FT(0.05), fill(FT(0.05), domain.space.surface))
        α_NIR_leaf_cases = (FT(0.45), fill(FT(0.45), domain.space.surface))
        τ_NIR_leaf_cases = (FT(0.25), fill(FT(0.25), domain.space.surface))
        χl_cases = (FT(0.1), fill(FT(0.1), domain.space.surface))
        zipped_params = zip(
            g1_cases,
            Vcmax25_cases,
            mechanism_cases,
            rooting_cases,
            α_PAR_leaf_cases,
            τ_PAR_leaf_cases,
            α_NIR_leaf_cases,
            τ_NIR_leaf_cases,
            χl_cases,
        )
        for (
            g1,
            Vcmax25,
            is_c3,
            rooting_depth,
            α_PAR_leaf,
            τ_PAR_leaf,
            α_NIR_leaf,
            τ_NIR_leaf,
            χl,
        ) in zipped_params
            BeerLambertparams = BeerLambertParameters(FT)
            # TwoStreamModel parameters
            Ω = FT(0.69)
            G_Function = CLMGFunction(χl)
            λ_γ_PAR = FT(5e-7)
            ϵ_canopy = FT(0.97)
            BeerLambertparams =
                BeerLambertParameters(FT; α_PAR_leaf, α_NIR_leaf, λ_γ_PAR)
            TwoStreamparams = TwoStreamParameters(
                FT;
                Ω,
                α_PAR_leaf,
                τ_PAR_leaf,
                α_NIR_leaf,
                τ_NIR_leaf,
                G_Function,
            )
            photosynthesis_params = FarquharParameters(FT, is_c3; Vcmax25)
            stomatal_g_params = MedlynConductanceParameters(FT; g1)

            stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
            photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
            rt_models = (
                BeerLambertModel{FT}(BeerLambertparams),
                TwoStreamModel{FT}(TwoStreamparams),
            )
            energy_model = BigLeafEnergyModel{FT}(BigLeafEnergyParameters{FT}())
            earth_param_set = LP.LandParameters(FT)
            LAI = FT(0.0) # m2 [leaf] m-2 [ground]
            z_0m = FT(2.0) # m, Roughness length for momentum - value from tall forest ChatGPT
            z_0b = FT(0.1) # m, Roughness length for scalars - value from tall forest ChatGPT
            h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
            shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
                z_0m,
                z_0b,
                earth_param_set,
            )
            lat = FT(0.0) # degree
            long = FT(-180) # degree

            function zenith_angle(
                t,
                start_date;
                latitude = lat,
                longitude = long,
                insol_params = earth_param_set.insol_params,
            )
                current_datetime = start_date + Dates.Second(round(t))
                d, δ, η_UTC =
                    FT.(
                        Insolation.helper_instantaneous_zenith_angle(
                            current_datetime,
                            start_date,
                            insol_params,
                        )
                    )
                return Insolation.instantaneous_zenith_angle(
                    d,
                    δ,
                    η_UTC,
                    longitude,
                    latitude,
                )[1]
            end

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
            start_date = DateTime(2005)
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
            )

            # Plant Hydraulics
            RAI = FT(0)
            SAI = FT(0)
            lai_fun = t -> LAI
            ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
                TimeVaryingInput(lai_fun),
                SAI,
                RAI,
            )
            K_sat_plant = FT(1.8e-8) # m/s
            ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
            Weibull_param = FT(4) # unitless, Holtzman's original c param value
            a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
            plant_ν = FT(0.7) # m3/m3
            plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
            conductivity_model =
                PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
            retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
            root_depths =
                SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0) # 1st element is the deepest root depth
            param_set = PlantHydraulics.PlantHydraulicsParameters(;
                ai_parameterization = ai_parameterization,
                ν = plant_ν,
                S_s = plant_S_s,
                rooting_depth = rooting_depth,
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
                        range(
                            start = 0.0,
                            step = Δz,
                            stop = Δz * (n_stem + n_leaf),
                        ),
                    )
                )

            ψ_soil0 = FT(0.0)
            T_soil0 = FT(290)
            soil_driver = PrescribedGroundConditions(
                root_depths,
                (t) -> ψ_soil0,
                (t) -> T_soil0,
                FT(0.2),
                FT(0.4),
                FT(0.98),
            )

            plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
                parameters = param_set,
                n_stem = n_stem,
                n_leaf = n_leaf,
                compartment_surfaces = compartment_faces,
                compartment_midpoints = compartment_centers,
            )
            autotrophic_parameters = AutotrophicRespirationParameters(FT)
            autotrophic_respiration_model =
                AutotrophicRespirationModel{FT}(autotrophic_parameters)
            for rt_model in rt_models
                canopy = ClimaLand.Canopy.CanopyModel{FT}(;
                    parameters = shared_params,
                    domain = domain,
                    radiative_transfer = rt_model,
                    photosynthesis = photosynthesis_model,
                    conductance = stomatal_model,
                    autotrophic_respiration = autotrophic_respiration_model,
                    energy = energy_model,
                    hydraulics = plant_hydraulics,
                    boundary_conditions = Canopy.AtmosDrivenCanopyBC(
                        atmos,
                        radiation,
                        soil_driver,
                    ),
                )

                Y, p, coords = ClimaLand.initialize(canopy)
                Y.canopy.hydraulics .= plant_ν
                Y.canopy.energy.T = FT(289)
                set_initial_cache! = make_set_initial_cache(canopy)
                t0 = FT(0.0)
                set_initial_cache!(p, Y, t0)

                @test all(parent(p.canopy.hydraulics.fa.:1) .== FT(0))
                @test all(
                    parent(p.canopy.energy.turbulent_fluxes.lhf) .== FT(0),
                )
                @test all(
                    parent(p.canopy.energy.turbulent_fluxes.shf) .== FT(0),
                )
                @test all(parent(p.canopy.energy.fa_energy_roots) .== FT(0))
                @test all(parent(p.canopy.hydraulics.fa_roots) .== FT(0))
                @test all(
                    parent(p.canopy.energy.turbulent_fluxes.transpiration) .==
                    FT(0),
                )
                @test all(parent(p.canopy.radiative_transfer.LW_n) .== FT(0))
                @test all(parent(p.canopy.radiative_transfer.SW_n) .== FT(0))
                @test all(parent(p.canopy.radiative_transfer.par.abs) .== FT(0))
                @test all(parent(p.canopy.radiative_transfer.nir.abs) .== FT(0))
                @test all(parent(p.canopy.autotrophic_respiration.Ra) .== FT(0))
            end
        end
    end
end

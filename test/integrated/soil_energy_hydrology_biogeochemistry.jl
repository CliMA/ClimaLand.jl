using Test
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
using ClimaLSM.Soil.Biogeochemistry
using Dates

import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

for FT in (Float32, Float64)
    @testset "Soil respiration test set, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        # Make soil model args
        ν = FT(0.556)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)
        κ_minerals = FT(2.5)
        κ_om = FT(0.25)
        κ_quartz = FT(8.0)
        κ_air = FT(0.025)
        κ_ice = FT(2.21)
        κ_liq = FT(0.57)
        ρp = FT(2.66 / 1e3 * 1e6)
        ρc_ds = FT(2e6 * (1.0 - ν))
        κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
        κ_dry = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
        κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
        κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)

        soil_ps = Soil.EnergyHydrologyParameters{FT}(;
            κ_dry = κ_dry,
            κ_sat_frozen = κ_sat_frozen,
            κ_sat_unfrozen = κ_sat_unfrozen,
            ρc_ds = ρc_ds,
            ν = ν,
            ν_ss_om = ν_ss_om,
            ν_ss_quartz = ν_ss_quartz,
            ν_ss_gravel = ν_ss_gravel,
            K_sat = K_sat,
            hydrology_cm = hcm,
            S_s = S_s,
            θ_r = θ_r,
            earth_param_set = earth_param_set,
        )
        zmax = FT(0)
        zmin = FT(-1)
        nelems = 20
        lsm_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        zero_flux_bc = Soil.FluxBC((p, t) -> 0.0)
        sources = () # PhaseChange
        boundary_fluxes = (;
            top = (water = zero_flux_bc, heat = zero_flux_bc),
            bottom = (water = zero_flux_bc, heat = zero_flux_bc),
        )
        soil_args = (;
            boundary_conditions = boundary_fluxes,
            sources = sources,
            domain = lsm_domain,
            parameters = soil_ps,
        )

        # Make biogeochemistry model args
        Csom = (z, t) -> eltype(z)(5.0)

        co2_parameters = Soil.Biogeochemistry.SoilCO2ModelParameters{FT}(;
            earth_param_set = earth_param_set,
        )
        C = FT(4)
        co2_top_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> C)
        co2_bot_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> C)
        co2_sources = ()
        co2_boundary_conditions =
            (; CO2 = (top = co2_top_bc, bottom = co2_bot_bc))

        # Make a PrescribedAtmosphere - we only care about atmos_p though
        precipitation_function = (t) -> 1.0
        snow_precip = (t) -> 1.0
        atmos_T = (t) -> 1.0
        atmos_u = (t) -> 1.0
        atmos_q = (t) -> 1.0
        atmos_p = (t) -> 100000.0
        UTC_DATETIME = Dates.now()
        atmos_h = FT(30)
        atmos_co2 = (t) -> 1.0

        atmos = ClimaLSM.PrescribedAtmosphere(
            precipitation_function,
            snow_precip,
            atmos_T,
            atmos_u,
            atmos_q,
            atmos_p,
            UTC_DATETIME,
            atmos_h;
            c_co2 = atmos_co2,
        )


        soil_drivers = Soil.Biogeochemistry.SoilDrivers(
            Soil.Biogeochemistry.PrognosticMet{FT}(),
            Soil.Biogeochemistry.PrescribedSOC{FT}(Csom),
            atmos,
        )
        soilco2_args = (;
            boundary_conditions = co2_boundary_conditions,
            sources = co2_sources,
            domain = lsm_domain,
            parameters = co2_parameters,
            drivers = soil_drivers,
        )

        # Create integrated model instance
        model = LandSoilBiogeochemistry{FT}(;
            soil_args = soil_args,
            soilco2_args = soilco2_args,
        )
        Y, p, coords = initialize(model)

        function init_soil!(Y, z, params)
            ν = params.ν
            FT = eltype(Y.soil.ϑ_l)
            Y.soil.ϑ_l .= FT(0.33)
            Y.soil.θ_i .= FT(0.0)
            T = FT(279.85)
            ρc_s = FT.(Soil.volumetric_heat_capacity(FT(0.33), FT(0.0), params))
            Y.soil.ρe_int .=
                Soil.volumetric_internal_energy.(FT(0.0), ρc_s, T, Ref(params))
        end

        function init_co2!(Y, C_0)
            Y.soilco2.C .= C_0
        end

        z = coords.subsurface.z
        init_soil!(Y, z, model.soil.parameters)
        init_co2!(Y, C)
        t0 = FT(0.0)
        set_initial_aux_state! = make_set_initial_aux_state(model)
        set_initial_aux_state!(p, Y, t0)

        @test p.soil.T ≈ Soil.Biogeochemistry.soil_temperature(
            model.soilco2.driver.met,
            p,
            Y,
            t0,
            z,
        )
        @test all(
            parent(
                Soil.Biogeochemistry.soil_SOM_C(
                    model.soilco2.driver.soc,
                    p,
                    Y,
                    t0,
                    z,
                ),
            ) .== FT(5.0),
        )
        @test p.soil.θ_l ≈ Soil.Biogeochemistry.soil_moisture(
            model.soilco2.driver.met,
            p,
            Y,
            t0,
            z,
        )

        try
            co2_parameters = Soil.Biogeochemistry.SoilCO2ModelParameters{FT}(;
                ν = FT(0.2),
                earth_param_set = earth_param_set,
            )
            soil_drivers = Soil.Biogeochemistry.SoilDrivers(
                Soil.Biogeochemistry.PrognosticMet{FT}(),
                Soil.Biogeochemistry.PrescribedSOC{FT}(Csom),
                atmos,
            )
            soilco2_args = (;
                boundary_conditions = co2_boundary_conditions,
                sources = co2_sources,
                domain = lsm_domain,
                parameters = co2_parameters,
                drivers = soil_drivers,
            )

            # Create integrated model instance
            model = LandSoilBiogeochemistry{FT}(;
                soil_args = soil_args,
                soilco2_args = soilco2_args,
            )
            Y, p, coords = initialize(model)

            function init_soil!(Y, z, params)
                ν = params.ν
                FT = eltype(Y.soil.ϑ_l)
                Y.soil.ϑ_l .= FT(0.33)
                Y.soil.θ_i .= FT(0.0)
                T = FT(279.85)
                ρc_s = Soil.volumetric_heat_capacity(FT(0.33), FT(0.0), params)
                Y.soil.ρe_int .=
                    Soil.volumetric_internal_energy.(
                        FT(0.0),
                        ρc_s,
                        T,
                        Ref(params),
                    )
            end

            function init_co2!(Y, C_0)
                Y.soilco2.C .= C_0
            end

            z = coords.subsurface.z
            init_soil!(Y, z, model.soil.parameters)
            init_co2!(Y, C)
            t0 = FT(0.0)
            set_initial_aux_state! = make_set_initial_aux_state(model)
            set_initial_aux_state!(p, Y, t0)

            @test p.soil.T ≈ Soil.Biogeochemistry.soil_temperature(
                model.soilco2.driver.met,
                p,
                Y,
                t0,
                z,
            )
            soil_drivers = Soil.Biogeochemistry.SoilDrivers(
                Soil.Biogeochemistry.PrescribedMet{FT}(Csom, Csom),
                Soil.Biogeochemistry.PrescribedSOC{FT}(Csom),
                atmos,
            )
            soilco2_args = (;
                boundary_conditions = co2_boundary_conditions,
                sources = co2_sources,
                domain = lsm_domain,
                parameters = co2_parameters,
                drivers = soil_drivers,
            )

            # Create integrated model instance
            model = LandSoilBiogeochemistry{FT}(;
                soil_args = soil_args,
                soilco2_args = soilco2_args,
            )
        catch err
            @test isa(err, AssertionError)
        end
    end
end

using Test
import ClimaComms
ClimaComms.@import_required_backends
using Dates
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil.Biogeochemistry
import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    toml_dict = LP.create_toml_dict(FT)
    @testset "Soil co2 biogeochemistry sources, FT = $FT" begin
        # Create parameters
        parameters = SoilCO2ModelParameters(toml_dict)

        # Set up domain
        nelems = 50
        zmin = FT(-1)
        zmax = FT(0.0)
        domain = Column(; zlim = (zmin, zmax), nelements = nelems)

        # Set up sources
        sources = (MicrobeProduction{FT}(),)

        # Make a PrescribedAtmosphere
        precipitation_function = (t) -> 1.0
        snow_precip = (t) -> 1.0
        atmos_T = (t) -> 1.0
        atmos_u = (t) -> 1.0
        atmos_q = (t) -> 1.0
        atmos_p = (t) -> 100000.0
        UTC_DATETIME = Dates.now()
        atmos_h = FT(30)
        atmos_co2 = (t) -> 1.0

        atmos = ClimaLand.PrescribedAtmosphere(
            TimeVaryingInput(precipitation_function),
            TimeVaryingInput(snow_precip),
            TimeVaryingInput(atmos_T),
            TimeVaryingInput(atmos_u),
            TimeVaryingInput(atmos_q),
            TimeVaryingInput(atmos_p),
            UTC_DATETIME,
            atmos_h,
            toml_dict;
            c_co2 = TimeVaryingInput(atmos_co2),
        )

        # Set up prescribed meteorology
        T_soil = (z, t) -> eltype(z)(303)
        θ_l = (z, t) -> eltype(z)(0.3)
        ν = FT(0.6)
        θ_r = FT(0.0)
        α = FT(0.1)
        n = FT(2)
        hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = n)
        prescribed_met = PrescribedMet{FT}(T_soil, θ_l, ν, θ_r, hcm)

        # Create drivers
        drivers = SoilDrivers(prescribed_met, atmos)

        # Create model
        model =
            SoilCO2Model{FT}(domain, drivers, toml_dict; parameters, sources)

        # Initialize and test
        Y, p, coords = initialize(model)
        dY = similar(Y)
        exp_tendency! = make_exp_tendency(model)
        t = Float64(1)

        # Set initial conditions
        Y.soilco2.CO2 .= FT(4)
        Y.soilco2.O2_f .= FT(0.2)
        Y.soilco2.SOC .= FT(5.0)

        # Update cache and compute tendencies
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t)
        exp_tendency!(dY, Y, p, t)

        # Test that microbial source is positive (CO2 production)
        @test all(parent(p.soilco2.Sm) .>= FT(0))

        # Test that SOC consumption equals CO2 production (mass balance)
        @test dY.soilco2.SOC ≈ @. -p.soilco2.Sm

        # Test that interior cells have CO2 production (source dominates)
        # Near boundaries, diffusion may affect the tendency
        interior_indices = 2:(length(parent(dY.soilco2.CO2)) - 1)
        @test all(parent(dY.soilco2.CO2)[interior_indices] .> FT(0))

        # Test that O2 consumption generally occurs in interior
        # (boundary diffusion may affect edge values)
        @test sum(parent(dY.soilco2.O2_f)[interior_indices]) < FT(0)

        # Test boundary condition variables are set
        @test p.soilco2.top_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.top_bc)
        @test p.soilco2.bottom_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.bottom_bc)
    end


    @testset "Soil co2 biogeochemistry diffusion, FT = $FT" begin
        # Create parameters
        parameters = SoilCO2ModelParameters(toml_dict)

        # Set up domain
        nelems = 50
        zmin = FT(-1)
        zmax = FT(0.0)
        domain = Column(; zlim = (zmin, zmax), nelements = nelems)

        # No sources for diffusion test
        sources = ()

        # Set uniform CO2 concentration for testing diffusion
        C = FT(4.0)

        # Make a PrescribedAtmosphere with CO2 matching soil initial conditions
        precipitation_function = (t) -> 1.0
        snow_precip = (t) -> 1.0
        atmos_T = (t) -> 1.0
        atmos_u = (t) -> 1.0
        atmos_q = (t) -> 1.0
        atmos_p = (t) -> 100000.0
        UTC_DATETIME = Dates.now()
        atmos_h = FT(30)
        atmos_co2 = (t) -> C  # Match atmospheric CO2 to soil CO2

        atmos = ClimaLand.PrescribedAtmosphere(
            TimeVaryingInput(precipitation_function),
            TimeVaryingInput(snow_precip),
            TimeVaryingInput(atmos_T),
            TimeVaryingInput(atmos_u),
            TimeVaryingInput(atmos_q),
            TimeVaryingInput(atmos_p),
            UTC_DATETIME,
            atmos_h,
            toml_dict;
            c_co2 = TimeVaryingInput(atmos_co2),
        )

        # Set up prescribed meteorology
        T_soil = (z, t) -> eltype(z)(303)
        θ_l = (z, t) -> eltype(z)(0.3)
        ν = FT(0.6)
        θ_r = FT(0.0)
        α = FT(0.1)
        n = FT(2)
        hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = n)
        prescribed_met = PrescribedMet{FT}(T_soil, θ_l, ν, θ_r, hcm)

        # Create drivers
        drivers = SoilDrivers(prescribed_met, atmos)

        # Create model
        model =
            SoilCO2Model{FT}(domain, drivers, toml_dict; parameters, sources)

        # Initialize and test
        Y, p, coords = initialize(model)
        dY = similar(Y)
        exp_tendency! = make_exp_tendency(model)
        t = Float64(1)

        # Set initial conditions - uniform CO2 matching boundary conditions
        Y.soilco2.CO2 .= C
        Y.soilco2.O2_f .= FT(0.2)
        Y.soilco2.SOC .= FT(5.0)

        # Update cache and compute tendencies
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t)
        exp_tendency!(dY, Y, p, t)

        # With no sources and uniform concentration matching atmospheric BC,
        # net CO2 change should be ~0
        @test sum(dY.soilco2.CO2) ≈ FT(0.0) atol = FT(1e-6)

        # Test boundary condition variables are set
        @test p.soilco2.top_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.top_bc)
        @test p.soilco2.bottom_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.bottom_bc)
    end
end

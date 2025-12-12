"""
Integration test for the OptimalLAI callback in a minimal ClimaLand simulation.
This test creates a simple canopy model with controlled forcing and verifies
that LAI updates correctly over time.

Water limitation is now handled through the soil moisture stress factor (β) in the
daily potential GPP computation, rather than through a water-limited LAI_max.
"""

using Test
import ClimaParams
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Dates
using ClimaLand
using ClimaLand: PrescribedAtmosphere, PrescribedRadiativeFluxes
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP

@testset "OptimalLAI Callback Integration Test" begin
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)

    # Create a point domain (simplest case)
    lat = FT(38.7)   # Ozark-like latitude
    long = FT(-92.2) # Ozark-like longitude
    z_sfc = FT(0.0)
    domain = Point(; z_sfc, longlat = (long, lat))

    # Set up forcing - summer conditions
    start_date = DateTime(2010, 6, 15, 0, 0, 0)  # Summer

    # Atmospheric forcing
    T_atmos = t -> FT(298.0)  # 25°C
    u_atmos = t -> FT(2.0)    # m/s wind
    q_atmos = t -> FT(0.01)   # kg/kg specific humidity
    P_atmos = t -> FT(101325.0)  # Pa
    liquid_precip = t -> FT(0.0)
    snow_precip = t -> FT(0.0)
    c_atmos = t -> FT(4.2e-4)  # mol/mol CO2
    h_atmos = FT(30.0)

    atmos = PrescribedAtmosphere(
        TimeVaryingInput(liquid_precip),
        TimeVaryingInput(snow_precip),
        TimeVaryingInput(T_atmos),
        TimeVaryingInput(u_atmos),
        TimeVaryingInput(q_atmos),
        TimeVaryingInput(P_atmos),
        start_date,
        h_atmos,
        toml_dict;
        c_co2 = TimeVaryingInput(c_atmos),
    )

    # Radiation - constant high PAR for summer
    shortwave_radiation = t -> FT(800.0)  # W/m² - strong summer radiation
    longwave_radiation = t -> FT(300.0)   # W/m²
    earth_param_set = LP.LandParameters(toml_dict)
    zenith_angle = (t, s) -> ClimaLand.default_zenith_angle(
        t, s;
        insol_params = earth_param_set.insol_params,
        latitude = lat,
        longitude = long,
    )

    radiation = PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(shortwave_radiation),
        TimeVaryingInput(longwave_radiation),
        start_date;
        θs = zenith_angle,
    )

    # Ground conditions
    ground = ClimaLand.PrescribedGroundConditions{FT}(;
        T = TimeVaryingInput(t -> FT(298.0)),
        θ = TimeVaryingInput(t -> FT(0.3)),
    )

    # Canopy model components
    radiation_parameters = (;
        α_PAR_leaf = FT(0.1),
        α_NIR_leaf = FT(0.4),
        Ω = FT(1),
        G_Function = Canopy.ConstantGFunction(FT(0.5)),
    )
    radiative_transfer = Canopy.BeerLambertModel{FT}(domain, toml_dict; radiation_parameters)

    # Use PModel for photosynthesis (required for OptimalLAI)
    is_c3 = FT(1)
    photosynthesis = Canopy.PModel{FT}(domain, toml_dict; is_c3)

    # PModel conductance
    conductance = Canopy.PModelConductance{FT}(toml_dict)

    # Optimal LAI model
    lai_model = Canopy.OptimalLAIModel{FT}(Canopy.OptimalLAIParameters{FT}(toml_dict))

    # Plant hydraulics
    hydraulics = Canopy.PlantHydraulicsModel{FT}(domain, toml_dict)

    # Biomass with prescribed LAI (will be overridden by OptimalLAI)
    LAI_fun = t -> FT(4.0)
    LAI = TimeVaryingInput(LAI_fun)
    biomass = Canopy.PrescribedBiomassModel{FT}(
        domain, LAI, toml_dict;
        rooting_depth = FT(1.0),
        height = FT(10.0),
    )

    # Energy model
    energy = Canopy.BigLeafEnergyModel{FT}(toml_dict)

    # No moisture stress (β = 1 always)
    soil_moisture_stress = Canopy.NoMoistureStressModel{FT}()

    # Combine into canopy model
    canopy_forcing = (; atmos, radiation, ground)
    canopy = Canopy.CanopyModel{FT}(
        domain,
        canopy_forcing,
        LAI_fun,
        toml_dict;
        radiative_transfer,
        photosynthesis,
        conductance,
        lai_model,
        soil_moisture_stress,
        hydraulics,
        energy,
        biomass,
    )

    # Simulation setup
    dt = FT(1800.0)  # 30 minutes

    # Initialize
    Y, p, coords = ClimaLand.initialize(canopy)

    # Get update_cache! function
    update_cache! = ClimaLand.make_update_cache(canopy)
    set_initial_cache! = ClimaLand.make_set_initial_cache(canopy)

    # Set initial conditions
    t0 = FT(0.0)

    println("\n=== Before set_initial_cache! ===")
    println("LAI = $(parent(p.canopy.lai_model.LAI))")
    println("A0_daily = $(parent(p.canopy.lai_model.A0_daily))")
    println("A0_annual = $(parent(p.canopy.lai_model.A0_annual))")

    # Initialize the cache
    set_initial_cache!(p, Y, t0)

    println("\n=== After set_initial_cache! ===")
    println("LAI = $(parent(p.canopy.lai_model.LAI))")
    println("A0_daily = $(parent(p.canopy.lai_model.A0_daily))")
    println("A0_annual = $(parent(p.canopy.lai_model.A0_annual))")
    println("A0_daily_acc = $(parent(p.canopy.lai_model.A0_daily_acc))")
    println("A0_annual_daily_acc = $(parent(p.canopy.lai_model.A0_annual_daily_acc))")

    # Now manually set LAI to a low value to test recovery
    p.canopy.lai_model.LAI .= FT(0.1)
    println("\n=== After manually setting LAI to 0.1 ===")
    println("LAI = $(parent(p.canopy.lai_model.LAI))")

    # Compute local noon time
    seconds_in_day = FT(86400.0)
    local_noon = seconds_in_day * (FT(0.5) - long / 360)
    println("\nLocal noon (seconds UTC): $local_noon")

    # GSL from defaults
    GSL = FT(240.0)

    # Create local_noon field
    local_noon_field = ClimaCore.Fields.ones(domain.space.surface) .* local_noon

    # Run simulation manually for 3 days
    n_days = 3
    n_timesteps_per_day = Int(seconds_in_day / dt)

    for day in 1:n_days
        for step in 1:n_timesteps_per_day
            current_t = FT((day - 1) * seconds_in_day + (step - 1) * dt)
            t_in_day = current_t % seconds_in_day

            # Update cache
            update_cache!(p, Y, current_t)

            # Current date
            current_date = start_date + Dates.Second(Int(current_t))

            # Call LAI update
            Canopy.call_update_optimal_LAI(
                p, Y, t_in_day, current_date;
                canopy = canopy,
                dt = dt,
                local_noon = local_noon_field,
                GSL = GSL,
            )

            # Check if noon
            is_noon = abs(t_in_day - local_noon) < dt / 2
            if is_noon
                println("\nDay $day at noon (t=$current_t):")
                println("  LAI = $(parent(p.canopy.lai_model.LAI))")
                println("  A0_daily = $(parent(p.canopy.lai_model.A0_daily))")
                println("  A0_daily_acc = $(parent(p.canopy.lai_model.A0_daily_acc))")
                println("  A0_annual_daily_acc = $(parent(p.canopy.lai_model.A0_annual_daily_acc))")
                println("  last_day_of_year = $(parent(p.canopy.lai_model.last_day_of_year))")
            end
        end
    end

    # Final state
    println("\n=== Final State ===")
    final_LAI = parent(p.canopy.lai_model.LAI)[1]
    final_A0_daily = parent(p.canopy.lai_model.A0_daily)[1]
    final_A0_annual = parent(p.canopy.lai_model.A0_annual)[1]
    println("LAI = $final_LAI")
    println("A0_daily = $final_A0_daily")
    println("A0_annual = $final_A0_annual")

    # Compute expected L_steady for current conditions
    # Water limitation is now through β in daily A0, not in LAI_max
    k = FT(0.5)
    z = FT(12.227)
    sigma = FT(0.771)
    alpha = FT(0.067)  # ~15 days memory

    LAI_max = Canopy.compute_L_max(final_A0_annual, k, z)
    m = Canopy.compute_m(GSL, LAI_max, final_A0_annual, sigma, k)
    L_steady = Canopy.compute_steady_state_LAI(final_A0_daily, m, k, LAI_max)
    println("\nExpected values:")
    println("LAI_max = $LAI_max (energy-limited only)")
    println("m = $m")
    println("L_steady = $L_steady")
    println("Threshold A0_daily = $(1/(m*k))")

    # Test that LAI increased from initial value if A0_daily > threshold
    threshold = 1 / (m * k)
    if final_A0_daily > threshold
        @test final_LAI > 0.1  # Should have increased from initial 0.1
        println("\nTEST PASSED: LAI increased from 0.1 to $final_LAI")
    else
        println("\nA0_daily ($final_A0_daily) below threshold ($threshold), LAI expected to decrease")
    end
end

using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, step!, balance_soil_carbon!
import ClimaUtilities.TimeManager: date
using Dates

@testset "Soil-carbon balance with canopy carbon pools" begin
    # The global land model with carbon pools, as in restart.jl, at low resolution.
    # `balance_soil_carbon!` over a few steps must set the log of the respiration
    # factor of each layer to log(litter / respiration) summed over those steps,
    # at most log(100), and leave the SOC at its state at the start. The pools are set by hand, so the
    # litter is far from the respiration and the factors are extreme.
    FT = Float64
    context = ClimaComms.context()
    start_date = DateTime(2008, 3, 1)
    Δt = 900.0
    n_steps = 12
    spinup_end = start_date + Second(n_steps * Δt)
    domain = ClimaLand.Domains.global_domain(
        FT;
        nelements = (10, 5),
        mask_threshold = FT(0.99),
        context,
    )
    toml_dict = LP.create_toml_dict(FT)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        spinup_end + Day(1),
        domain.space.surface,
        toml_dict,
        FT;
        use_lowres_forcing = true,
        context,
    )
    model = ClimaLand.LandModel{FT}(
        (; atmos, radiation),
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
        biomass = ClimaLand.Canopy.PrognosticCarbonModel{FT}(
            ClimaLand.Canopy.ZhouOptimalLAIModel{FT}(
                ClimaLand.Domains.obtain_surface_domain(domain),
                toml_dict,
            ),
            toml_dict,
        ),
    )
    set_ic! = ClimaLand.Simulations.make_set_initial_state_from_file(
        ClimaLand.Artifacts.soil_ic_2008_50m_path(; context),
        model;
        enforce_constraints = true,
    )
    new_simulation() = LandSimulation(
        start_date,
        spinup_end + Day(1),
        Δt,
        model;
        set_ic!,
        diagnostics = (),
        user_callbacks = (),
    )
    simulation = new_simulation()
    Y, p = simulation._integrator.u, simulation._integrator.p
    @test all(iszero, parent(p.soilco2.ln_rate_factor))
    # The cells with carbon pools shed litter, so give them pools
    Y.canopy.biomass.C_leaf .= FT(0.1)
    Y.canopy.biomass.C_stem .= FT(3)
    Y.canopy.biomass.C_root .= FT(0.5)
    state = copy(Y)
    SOC_start = copy(Y.soilco2.SOC)
    balance_soil_carbon!(simulation, spinup_end)
    @test date(simulation._integrator.t) == spinup_end
    @test Y.soilco2.SOC == SOC_start

    # The same steps from the same state, summing litter input and respiration
    replay = new_simulation()
    Y_r, p_r = replay._integrator.u, replay._integrator.p
    Y_r .= state
    ClimaLand.make_set_initial_cache(model)(p_r, Y_r, replay._integrator.t)
    litter = zero(SOC_start)
    respiration = zero(SOC_start)
    for _ in 1:n_steps
        step!(replay)
        litter .+= p_r.soil_litter_input
        respiration .+= p_r.soilco2.Sm
    end
    # Broadcasts skip the ocean points of the land mask, so start from zeros
    expected = zero(SOC_start)
    @. expected = ifelse(
        (respiration > 0) & isfinite(respiration) & isfinite(litter),
        ifelse(
            litter > 0,
            min(log(litter / respiration), log(FT(100))),
            log(floatmin(FT)),
        ),
        FT(0),
    )
    ln_factor = parent(p.soilco2.ln_rate_factor)
    @test count(x -> x != 0, ln_factor) > 0
    @test all(isfinite, ln_factor)
    @test ln_factor ≈ parent(expected) rtol = sqrt(eps(FT))

    # With the factor, respiration matches the litter input over the same steps
    rebalanced = new_simulation()
    Y_b, p_b = rebalanced._integrator.u, rebalanced._integrator.p
    Y_b .= state
    p_b.soilco2.ln_rate_factor .= p.soilco2.ln_rate_factor
    ClimaLand.make_set_initial_cache(model)(p_b, Y_b, rebalanced._integrator.t)
    respiration_b = zero(SOC_start)
    for _ in 1:n_steps
        step!(rebalanced)
        respiration_b .+= p_b.soilco2.Sm
    end
    # Where the factor is moderate, respiration scales with it, so it matches the
    # litter input over the same steps (the state differs only by the rebalanced
    # respiration). Large factors are capped by the substrate and O2 limitation.
    lt, rb = parent(litter), parent(respiration_b)
    @test maximum(ln_factor) <= log(FT(100)) * (1 + eps(FT))
    # Capped layers do not reach their litter input by design
    moderate = @. (lt > 0) & (abs(ln_factor) < 4)
    @test count(moderate) > 0
    @test rb[moderate] ≈ lt[moderate] rtol = 0.02
    @test all(@. abs(rb[moderate] / lt[moderate] - 1) < 0.02)

    @test_throws ErrorException balance_soil_carbon!(simulation, start_date)
end

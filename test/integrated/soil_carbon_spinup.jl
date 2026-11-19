using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, step!
import ClimaLand.Soil.Biogeochemistry as BGC
import ClimaUtilities.TimeManager: date
using Dates

@testset "Soil carbon spin-up of the land model" begin
    # A coarse version of the long run with prognostic LAI and carbon pools
    FT = Float64
    context = ClimaComms.context()
    start_date = DateTime(2008, 3, 1)
    Δt = 900.0
    stop_date = start_date + Day(1)
    domain = ClimaLand.Domains.global_domain(
        FT;
        nelements = (10, 5),
        mask_threshold = FT(0.99),
        context,
    )
    toml_dict = LP.create_toml_dict(FT)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
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
    simulation =
        LandSimulation(start_date, stop_date, Δt, model; diagnostics = nothing)
    Y, p = simulation._integrator.u, simulation._integrator.p
    # Pools from their equilibrium over a few hours, so that they shed litter
    ClimaLand.Simulations.spin_up_carbon_pools!(
        simulation,
        start_date + Hour(6),
    )
    step!(simulation)

    # The sums of one sample give back the litter at the equilibrium they set
    soilco2 = model.soilco2
    sums = ClimaLand.Simulations.SoilCarbonSums(Y)
    ClimaLand.Simulations.add_soil_carbon_sums!(
        sums,
        soilco2,
        Y,
        p,
        simulation._integrator.t,
    )
    @test sums.litter == p.soil_litter_input
    SOC = copy(Y.soilco2.SOC)
    SOC .= -1 # layers without an equilibrium keep it
    ClimaLand.Simulations.soil_carbon_equilibrium!(
        SOC,
        sums.litter,
        sums.respiration,
    )
    @test all(isfinite, parent(SOC))
    params = soilco2.parameters
    (; z, z_sfc) = soilco2.domain.fields
    θ_l = BGC.soil_moisture(
        soilco2.drivers.met,
        p,
        Y,
        simulation._integrator.t,
        z,
    )
    # Broadcasts skip the ocean columns of the masked domain, which stay 0 here
    space = axes(SOC)
    misfit = ClimaCore.Fields.zeros(space)
    solved = ClimaCore.Fields.zeros(space)
    @. solved = ifelse(SOC > 0, FT(1), FT(0))
    @. misfit = ifelse(
        solved > 0,
        abs(
            BGC.microbe_source(
                p.soilco2.T,
                θ_l,
                SOC,
                sums.o2_availability,
                params,
            ) * BGC.decomposition_depth_factor(z_sfc - z, params) /
            sums.litter - 1,
        ),
        FT(0),
    )
    @test sum(parent(solved)) > 0
    @test maximum(parent(misfit)) < FT(0.02)

    # The spin-up sets the soil carbon and keeps the simulation running
    ClimaLand.Simulations.spin_up_carbon_pools!(
        simulation,
        start_date + Hour(10);
        soil_carbon = true,
    )
    @test date(simulation._integrator.t) == start_date + Hour(10)
    @test all(isfinite, parent(Y.soilco2.SOC))
    @test all(parent(Y.soilco2.SOC) .>= 0)
    step!(simulation)
    @test all(isfinite, parent(Y.soilco2.SOC))
end

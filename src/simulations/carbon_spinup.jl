"""
    spin_up_carbon_pools!(landsim::LandSimulation, spinup_end::DateTime)

Spins up the carbon pools of the `PrognosticCarbonModel` of the canopy of `landsim`:
advances `landsim` to `spinup_end`, averaging the drivers of the pools over its steps,
then sets the leaf, wood and fine-root pools to their equilibrium under these mean
drivers (`Canopy.equilibrium_carbon_pools`): the GPP, the leaf respiration, the Q10
factor of the maintenance respiration at the canopy temperature, the air temperature,
which also becomes the mean annual temperature `T_annual`, and the tree share.

The pools follow the GPP and do not feed back on it, so the drivers of the spin-up are
those of the simulation, which should span whole years. Without it, wood and fine roots
take decades and years to fill the empty pools of the initial conditions. The sugar pool
keeps its state, which follows the fluxes of the last days: after a large change of the
other pools, as from empty ones, it settles over weeks, so a second spin-up over the
following year leaves a state without transients. The steps of the spin-up are part of
the simulation and its diagnostics; a restart does not repeat it.
"""
function spin_up_carbon_pools!(landsim::LandSimulation, spinup_end::DateTime)
    model = landsim.model
    canopy = model isa ClimaLand.Canopy.CanopyModel ? model : model.canopy
    biomass = canopy.biomass
    biomass isa ClimaLand.Canopy.PrognosticCarbonModel ||
        error("The canopy of the simulation has no carbon pools")
    integrator = landsim._integrator
    Y, p = integrator.u, integrator.p
    date(integrator.t) < spinup_end ||
        error("The spin-up must end after the current date of the simulation")
    surface_space = axes(Y.canopy.biomass.C_leaf)
    sums = (;
        GPP = ClimaCore.Fields.zeros(surface_space),
        Rd = ClimaCore.Fields.zeros(surface_space),
        f_T = ClimaCore.Fields.zeros(surface_space),
        T_air = ClimaCore.Fields.zeros(surface_space),
        tree_share = ClimaCore.Fields.zeros(surface_space),
    )
    parameters = biomass.parameters
    (; Q10, T_ref) = parameters
    nsteps = 0
    while date(integrator.t) < spinup_end
        step!(landsim)
        GPP = ClimaLand.Canopy.get_GPP(p, canopy.photosynthesis)
        Rd = ClimaLand.Canopy.get_Rd_canopy(p, canopy.photosynthesis)
        T_canopy =
            ClimaLand.Canopy.canopy_temperature(canopy.energy, canopy, Y, p)
        tree_share = ClimaLand.Canopy.carbon_tree_share(biomass, p)
        @. sums.GPP += GPP
        @. sums.Rd += Rd
        @. sums.f_T += Q10^((T_canopy - T_ref) / 10)
        @. sums.T_air += p.drivers.T
        @. sums.tree_share += tree_share
        nsteps += 1
    end

    FT = eltype(parameters)
    n = FT(nsteps)
    pools = ClimaCore.Fields.Field(
        NamedTuple{(:C_sugar, :C_leaf, :C_stem, :C_root), NTuple{4, FT}},
        surface_space,
    )
    @. pools = ClimaLand.Canopy.equilibrium_carbon_pools(
        parameters,
        sums.GPP / n,
        sums.Rd / n,
        sums.f_T / n,
        sums.T_air / n,
        sums.tree_share / n,
    )
    b = Y.canopy.biomass
    @. b.C_leaf = finite_or_zero(pools.C_leaf)
    @. b.C_stem = finite_or_zero(pools.C_stem)
    @. b.C_root = finite_or_zero(pools.C_root)
    @. b.T_annual = sums.T_air / n
    return nothing
end

finite_or_zero(x) = ifelse(isfinite(x), x, zero(x))

"""
    balance_soil_carbon!(landsim::LandSimulation, spinup_end::DateTime)

Holds the soil organic carbon (SOC) of `landsim` at its current state, as an
observed stock, and sets the decomposition rate of each soil layer so that the
SOC is in steady state, as in data-constrained models such as CARDAMOM: advances
`landsim` to `spinup_end`, averaging the litter input `S_litter` and the
microbial respiration `Sm` of each layer over its steps, then multiplies the
respiration rate by `S_litter / Sm` (stored as its logarithm in
`p.soilco2.ln_rate_factor`), and puts the SOC back to its state at the start.

With the new rates, the mean respiration of each layer over the spin-up equals
its mean litter input, so the SOC has no trend under the drivers of the
spin-up; it still varies with the climate, the litter and the moisture of the
other years. A layer without litter input keeps its SOC (factor `floatmin`, in effect 0); a layer
without respiration keeps its rate. The spin-up should span whole years and
follow `spin_up_carbon_pools!`, so the litter comes from pools in equilibrium.
The factor lives in the cache, not in the state: a simulation started from the
state must set `p.soilco2.ln_rate_factor` again (see `experiments/trendy`).
"""
function balance_soil_carbon!(landsim::LandSimulation, spinup_end::DateTime)
    model = landsim.model
    hasproperty(model, :soilco2) && hasproperty(model, :canopy) ||
        error("The simulation has no soil CO2 model and canopy")
    model.canopy.biomass isa ClimaLand.Canopy.PrognosticCarbonModel ||
        error("The canopy of the simulation has no carbon pools")
    integrator = landsim._integrator
    Y, p = integrator.u, integrator.p
    date(integrator.t) < spinup_end ||
        error("The spin-up must end after the current date of the simulation")
    SOC_start = copy(Y.soilco2.SOC)
    litter = zero(SOC_start)
    respiration = zero(SOC_start)
    while date(integrator.t) < spinup_end
        step!(landsim)
        @. litter += p.soil_litter_input
        @. respiration += p.soilco2.Sm
    end
    FT = eltype(SOC_start)
    ln_factor = p.soilco2.ln_rate_factor
    @. ln_factor = ifelse(
        (respiration > 0) & isfinite(respiration) & isfinite(litter),
        ifelse(
            litter > 0,
            ln_factor + log(litter / respiration),
            log(floatmin(FT)),
        ),
        ln_factor,
    )
    Y.soilco2.SOC .= SOC_start
    return nothing
end

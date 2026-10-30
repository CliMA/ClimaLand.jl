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

"""
    spin_up_carbon_pools!(landsim::LandSimulation, spinup_end::DateTime; soil_carbon = false)

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

With `soil_carbon = true`, the soil organic carbon (SOC) of a `LandModel`, fed by the
litter of the pools, is also set to its equilibrium. Every fourth step, each soil layer
sums its litter input and, for each candidate SOC value of `SOIL_CARBON_GRID`, the
microbial respiration that its temperature, moisture and O₂ would give at that SOC. At
the end, each layer's SOC is set to the value at which its mean respiration equals its
mean litter input (`soil_carbon_equilibrium!`). SOC does not feed back on the vegetation
or the soil physics, so one period is enough. The exceptions are the soil O₂ and the
litter of the pools set at the end, which adjust afterwards. Over centuries, the
respiration of deep layers, slowed by `decomposition_depth_factor`, would bring them to
the same state.
"""
function spin_up_carbon_pools!(
    landsim::LandSimulation,
    spinup_end::DateTime;
    soil_carbon = false,
)
    model = landsim.model
    canopy = model isa ClimaLand.Canopy.CanopyModel ? model : model.canopy
    biomass = canopy.biomass
    biomass isa ClimaLand.Canopy.PrognosticCarbonModel ||
        error("The canopy of the simulation has no carbon pools")
    soil_carbon &&
        !(model isa ClimaLand.LandModel && hasproperty(model, :soilco2)) &&
        error("The soil carbon spin-up needs a LandModel with soil CO2")
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
    soil_sums = soil_carbon ? SoilCarbonSums(Y) : nothing
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
        soil_carbon &&
            nsteps % 4 == 0 &&
            add_soil_carbon_sums!(soil_sums, model.soilco2, Y, p, integrator.t)
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
    soil_carbon && soil_carbon_equilibrium!(
        Y.soilco2.SOC,
        soil_sums.litter,
        soil_sums.respiration,
    )
    return nothing
end

finite_or_zero(x) = ifelse(isfinite(x), x, zero(x))

"""
    SOIL_CARBON_GRID

The candidate SOC values (kg C m⁻³) of the soil carbon spin-up, log-spaced.
"""
const SOIL_CARBON_GRID = exp10.(range(-3, 3; length = 49))

"""
    SOIL_CARBON_MAX

The largest SOC (kg C m⁻³) the soil carbon spin-up sets; a layer whose equilibrium
would exceed it keeps its SOC.
"""
const SOIL_CARBON_MAX = 100.0

"""
    SoilCarbonSums(Y)

The sums of the soil carbon spin-up over the soil layers of `Y`: the litter input, and a
tuple of the microbial respiration at each SOC value of `SOIL_CARBON_GRID`; and a scratch
field for the O₂ availability.
"""
struct SoilCarbonSums{F, R}
    litter::F
    respiration::R
    o2_availability::F
end

function SoilCarbonSums(Y)
    space = axes(Y.soilco2.SOC)
    return SoilCarbonSums(
        ClimaCore.Fields.zeros(space),
        ntuple(_ -> ClimaCore.Fields.zeros(space), length(SOIL_CARBON_GRID)),
        ClimaCore.Fields.zeros(space),
    )
end

"""
    add_soil_carbon_sums!(sums::SoilCarbonSums, soilco2, Y, p, t)

Adds the litter input of `p` and, for each SOC value of `SOIL_CARBON_GRID`, the
microbial respiration of `soilco2` at that SOC under the soil temperature, moisture and
O₂ of `Y` and `p`, as in its `update_aux!`.
"""
function add_soil_carbon_sums!(sums::SoilCarbonSums, soilco2, Y, p, t)
    BGC = ClimaLand.Soil.Biogeochemistry
    params = soilco2.parameters
    met = soilco2.drivers.met
    (; z, z_sfc) = soilco2.domain.fields
    T_soil = p.soilco2.T
    θ_l = BGC.soil_moisture(met, p, Y, t, z)
    θ_i = BGC.soil_ice(met, p, Y, t, z)
    P_sfc = p.drivers.P
    O2_avail = sums.o2_availability
    @. O2_avail = BGC.o2_availability(
        BGC.o2_fraction_from_concentration(
            max(Y.soilco2.O2, 0) / p.soilco2.θ_eff_o2,
            T_soil,
            P_sfc,
            params,
        ),
        BGC.volumetric_air_content(θ_l + θ_i, met.ν),
        params.D_oa,
    )
    FT = eltype(Y.soilco2.SOC)
    for (respiration, c) in zip(sums.respiration, SOIL_CARBON_GRID)
        @. respiration +=
            BGC.microbe_source(T_soil, θ_l, FT(c), O2_avail, params) *
            BGC.decomposition_depth_factor(z_sfc - z, params)
    end
    @. sums.litter += p.soil_litter_input
    return nothing
end

"""
    soil_carbon_equilibrium!(SOC, litter, respiration; grid = SOIL_CARBON_GRID, soc_max = SOIL_CARBON_MAX)

Sets `SOC` in each layer to the value at which the respiration equals the litter input,
from their sums over the same samples: `respiration[k]` sums the respiration at
SOC `grid[k]`. The respiration is interpolated linearly between the candidates in log
SOC, and taken linear in SOC below the first one. Layers without litter get no SOC.
Layers whose litter exceeds the respiration at `soc_max`, or at the last candidate,
keep their SOC.
"""
function soil_carbon_equilibrium!(
    SOC,
    litter,
    respiration;
    grid = SOIL_CARBON_GRID,
    soc_max = SOIL_CARBON_MAX,
)
    FT = eltype(SOC)
    soc = similar(SOC)
    R1 = respiration[1]
    c1 = FT(grid[1])
    @. soc = ifelse(R1 > 0, c1 * litter / R1, FT(Inf))
    for k in 1:(length(grid) - 1)
        lo, hi = log(FT(grid[k])), log(FT(grid[k + 1]))
        Rk, Rk1 = respiration[k], respiration[k + 1]
        @. soc = ifelse(
            (litter >= Rk) & (litter < Rk1),
            exp(lo + (litter - Rk) / (Rk1 - Rk) * (hi - lo)),
            soc,
        )
    end
    RK = respiration[end]
    @. SOC = ifelse(
        litter > 0,
        ifelse((litter < RK) & (soc <= FT(soc_max)), soc, SOC),
        FT(0),
    )
    return nothing
end

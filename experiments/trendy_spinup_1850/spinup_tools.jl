import ClimaComms
import ClimaCore
import ClimaCore: InputOutput
import ClimaLand
import ClimaLand.Canopy
import ClimaLand.Soil.Biogeochemistry as BGC

"""Candidate SOC values (kg C m⁻³) of the respiration curves, log-spaced."""
const SOC_GRID = exp10.(range(-3, 3; length = 49))
soc_key(k) = Symbol("c", lpad(k, 2, '0'))

"""
    SpinupAccumulator(simulation; soc = true, substep = 4)

Time sums over a spin-up cycle, kept in a `FieldVector` so that they can be saved and
resumed across jobs:

- `pools`: the drivers of `Canopy.equilibrium_carbon_pools`, as in
  `ClimaLand.Simulations.spin_up_carbon_pools!` (GPP, leaf respiration, Q10 factor of the
  canopy temperature, air temperature, tree share), summed every step;
- `soc` (with `soc = true`): in every soil layer, the litter input, the O2 availability
  and, for each candidate SOC value `SOC_GRID[k]`, the microbial respiration
  `microbe_source` would have at that SOC under the step's temperature, moisture and O2.
  These are summed every `substep` steps, weighted by the time they stand for.
- `alts` (with `soc = true` and a non-empty `alts`): the same curves for alternative
  DAMM parameters, one set per element `(; ea, kmo2, doa)` of `alts`: the activation
  energy `ea` (J/mol, or `NaN` for the model's), and multipliers of `kM_o2` and `D_oa`
  (the O2 availability is recomputed with `D_oa × doa` from the state's O2). Changes of
  `V_ref` and of the SOC scale `kM/(p_sx D_liq)` need no extra set: they act exactly on
  any curve.

The time means give the equilibrium vegetation pools and, per layer, the SOC at which the
mean respiration equals the mean litter input (`solve_soc`).
"""
mutable struct SpinupAccumulator{V, S}
    acc::V
    pool_seconds::Float64
    soc_seconds::Float64
    nstep::Int
    substep::Int
    soc::Bool
    alts::Vector{NamedTuple{(:ea, :kmo2, :doa), NTuple{3, Float64}}}
    alt_params::Vector{Any}
    scratch::S
end

"""The DAMM parameters `params` with the alternative `alt` of a `SpinupAccumulator`."""
function alternative_parameters(params, alt)
    T = typeof(params)
    values = map(fieldnames(T)) do name
        v = getfield(params, name)
        name == :Ea_sx && !isnan(alt.ea) ? oftype(v, alt.ea) :
        name == :kM_o2 ? v * alt.kmo2 : name == :D_oa ? v * alt.doa : v
    end
    return T(values...)
end

"""Parses `SOC_ALTS`, e.g. "ea=60000,kmo2=0.3;doa=3", into the `alts` of a `SpinupAccumulator`."""
function parse_alts(spec)
    alts = NamedTuple{(:ea, :kmo2, :doa), NTuple{3, Float64}}[]
    for item in filter(!isempty, split(spec, ';'))
        kv = Dict(Symbol(strip(k)) => parse(Float64, v) for (k, v) in split.(split(item, ','), '='))
        push!(alts, (; ea = get(kv, :ea, NaN), kmo2 = get(kv, :kmo2, 1.0), doa = get(kv, :doa, 1.0)))
    end
    return alts
end

alts_string(alts) = join(["ea=$(a.ea),kmo2=$(a.kmo2),doa=$(a.doa)" for a in alts], ";")

function SpinupAccumulator(simulation; soc = true, substep = 4, alts = NamedTuple{(:ea, :kmo2, :doa), NTuple{3, Float64}}[])
    Y = simulation._integrator.u
    surface = axes(Y.canopy.biomass.C_leaf)
    subsurface = axes(Y.soilco2.SOC)
    zs() = ClimaCore.Fields.zeros(surface)
    pools = ClimaCore.Fields.FieldVector(;
        GPP = zs(),
        Rd = zs(),
        f_T = zs(),
        T_air = zs(),
        tree_share = zs(),
    )
    acc = if soc
        new_curves() = ClimaCore.Fields.FieldVector(;
            NamedTuple{ntuple(soc_key, length(SOC_GRID))}(
                ntuple(_ -> ClimaCore.Fields.zeros(subsurface), length(SOC_GRID)),
            )...,
        )
        base = (;
            pools,
            litter = ClimaCore.Fields.zeros(subsurface),
            curves = new_curves(),
            o2 = ClimaCore.Fields.zeros(subsurface),
            o2_sum = ClimaCore.Fields.zeros(subsurface),
        )
        isempty(alts) ? ClimaCore.Fields.FieldVector(; base...) :
        ClimaCore.Fields.FieldVector(;
            base...,
            alts = ClimaCore.Fields.FieldVector(;
                NamedTuple{ntuple(j -> Symbol("a", j), length(alts))}(ntuple(_ -> new_curves(), length(alts)))...,
            ),
        )
    else
        ClimaCore.Fields.FieldVector(; pools)
    end
    scratch = soc && !isempty(alts) ? ClimaCore.Fields.zeros(subsurface) : nothing
    return SpinupAccumulator(acc, 0.0, 0.0, 0, substep, soc, collect(alts), Any[], scratch)
end

function accumulate_spinup!(sa::SpinupAccumulator, simulation, Δt)
    land = simulation.model
    canopy = land.canopy
    biomass = canopy.biomass
    integrator = simulation._integrator
    Y, p, t = integrator.u, integrator.p, integrator.t
    (; Q10, T_ref) = biomass.parameters
    s = sa.acc.pools
    GPP = Canopy.get_GPP(p, canopy.photosynthesis)
    Rd = Canopy.get_Rd_canopy(p, canopy.photosynthesis)
    T_canopy = Canopy.canopy_temperature(canopy.energy, canopy, Y, p)
    tree_share = Canopy.carbon_tree_share(biomass, p)
    @. s.GPP += GPP * Δt
    @. s.Rd += Rd * Δt
    @. s.f_T += Q10^((T_canopy - T_ref) / 10) * Δt
    @. s.T_air += p.drivers.T * Δt
    @. s.tree_share += tree_share * Δt
    sa.pool_seconds += Δt
    sa.nstep += 1
    if sa.soc && sa.nstep % sa.substep == 0
        w = sa.substep * Δt
        soilco2 = land.soilco2
        params = soilco2.parameters
        met = soilco2.drivers.met
        z = soilco2.domain.fields.z
        T_soil = p.soilco2.T
        θ_l = BGC.soil_moisture(met, p, Y, t, z)
        θ_i = BGC.soil_ice(met, p, Y, t, z)
        P_sfc = p.drivers.P
        ν = met.ν
        D_oa = params.D_oa
        O2_avail = sa.acc.o2
        @. O2_avail = BGC.o2_availability(
            BGC.o2_fraction_from_concentration(
                max(Y.soilco2.O2, 0) / p.soilco2.θ_eff_o2,
                T_soil,
                P_sfc,
                params,
            ),
            BGC.volumetric_air_content(θ_l + θ_i, ν),
            D_oa,
        )
        for (k, c) in enumerate(SOC_GRID)
            curve = getproperty(sa.acc.curves, soc_key(k))
            @. curve += BGC.microbe_source(T_soil, θ_l, c, O2_avail, params) * w
        end
        @. sa.acc.o2_sum += O2_avail * w
        isempty(sa.alt_params) && append!(sa.alt_params, [alternative_parameters(params, a) for a in sa.alts])
        for (j, (alt, alt_params)) in enumerate(zip(sa.alts, sa.alt_params))
            O2_alt = sa.scratch
            @. O2_alt = BGC.o2_availability(
                BGC.o2_fraction_from_concentration(
                    max(Y.soilco2.O2, 0) / p.soilco2.θ_eff_o2,
                    T_soil,
                    P_sfc,
                    params,
                ),
                BGC.volumetric_air_content(θ_l + θ_i, ν),
                D_oa * alt.doa,
            )
            set = getproperty(sa.acc.alts, Symbol("a", j))
            for (k, c) in enumerate(SOC_GRID)
                curve = getproperty(set, soc_key(k))
                @. curve += BGC.microbe_source(T_soil, θ_l, c, O2_alt, alt_params) * w
            end
        end
        @. sa.acc.litter += p.soil_litter_input * w
        sa.soc_seconds += w
    end
    return nothing
end

function save_spinup(sa::SpinupAccumulator, path, context)
    writer = InputOutput.HDF5Writer(path, context)
    InputOutput.write_attributes!(
        writer,
        "/",
        Dict(
            "pool_seconds" => sa.pool_seconds,
            "soc_seconds" => sa.soc_seconds,
            "nstep" => sa.nstep,
            "alts" => alts_string(sa.alts),
        ),
    )
    InputOutput.write!(writer, sa.acc, "acc")
    Base.close(writer)
end

function load_spinup!(sa::SpinupAccumulator, path, context)
    read_into!(sa.acc, path, "acc")
    sa.pool_seconds = read_root_attribute(path, "pool_seconds")
    sa.soc_seconds = read_root_attribute(path, "soc_seconds")
    sa.nstep = read_root_attribute(path, "nstep")
    return sa
end

finite_or_zero(x) = ifelse(isfinite(x), x, zero(x))

"""
    apply_pool_equilibrium!(sa, simulation)

Sets the leaf, wood and fine-root pools and `T_annual` to their equilibrium under the
mean drivers of the cycle, as `spin_up_carbon_pools!` does for its own period.
"""
function apply_pool_equilibrium!(sa::SpinupAccumulator, simulation)
    canopy = simulation.model.canopy
    parameters = canopy.biomass.parameters
    FT = eltype(parameters)
    Y = simulation._integrator.u
    s = sa.acc.pools
    T = FT(sa.pool_seconds)
    surface = axes(Y.canopy.biomass.C_leaf)
    pools = ClimaCore.Fields.Field(
        NamedTuple{(:C_sugar, :C_leaf, :C_stem, :C_root), NTuple{4, FT}},
        surface,
    )
    @. pools = Canopy.equilibrium_carbon_pools(
        parameters,
        s.GPP / T,
        s.Rd / T,
        s.f_T / T,
        s.T_air / T,
        s.tree_share / T,
    )
    b = Y.canopy.biomass
    @. b.C_leaf = finite_or_zero(pools.C_leaf)
    @. b.C_stem = finite_or_zero(pools.C_stem)
    @. b.C_root = finite_or_zero(pools.C_root)
    @. b.T_annual = s.T_air / T
    return nothing
end

"""
    solve_soc(sa, simulation; cap = 100.0)

Per soil layer, the SOC (kg C m⁻³) at which the mean respiration of the cycle equals its
mean litter input, with log respiration interpolated linearly in log SOC between the
candidate values (exact in the linear regime). Below the first
candidate the respiration is nearly linear in SOC, so the SOC scales with the litter.
Where the litter exceeds the respiration at the largest candidate, no equilibrium exists
in range and the SOC is set to `cap`; the returned mask marks those points.
"""
function solve_soc(sa::SpinupAccumulator, simulation; cap = 100.0)
    Y = simulation._integrator.u
    FT = eltype(Y.soilco2.SOC)
    T = FT(sa.soc_seconds)
    L = @. sa.acc.litter / T
    R(k) = getproperty(sa.acc.curves, soc_key(k))
    K = length(SOC_GRID)
    soc = similar(Y.soilco2.SOC)
    c1 = FT(SOC_GRID[1])
    R1, RK = R(1), R(K)
    @. soc = ifelse(R1 > 0, c1 * L / (R1 / T), FT(0))
    for k in 1:(K - 1)
        lo, hi = log(FT(SOC_GRID[k])), log(FT(SOC_GRID[k + 1]))
        Rk, Rk1 = R(k), R(k + 1)
        @. soc = ifelse(
            (L * T >= Rk) & (L * T < Rk1),
            ifelse(
                Rk > 0,
                exp(lo + log(L * T / Rk) / log(Rk1 / Rk) * (hi - lo)),
                exp(lo + (L * T - Rk) / (Rk1 - Rk) * (hi - lo)),
            ),
            soc,
        )
    end
    capped = @. FT(L * T >= RK)
    @. soc = ifelse(capped > 0, FT(cap), soc)
    @. soc = ifelse(L > 0, finite_or_zero(soc), FT(0))
    return soc, capped
end

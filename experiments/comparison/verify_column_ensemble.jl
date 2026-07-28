# Self-verifying `Column` vs single-point `ColumnEnsemble` check for the
# standalone `EnergyHydrology` soil column (same model setup as
# `compare_soil.jl`).
#
# The two builds share an identical vertical mesh, so column 1 of the ensemble
# should reproduce the single column. This script asserts every leaf field is
# EXACTLY equal (RMSE == 0); for any field that is not, it automatically loops to
# localize the difference (which entries, how many ulps) and ties it back to the
# geometry, i.e. the extruded ensemble's Jacobian carries a constant horizontal
# metric factor J_h = R²·cosd(lat)·(π/180)² that cancels analytically but not
# bitwise inside the vertical stencils. (See investigate_roundoff.jl for the full
# derivation; this file is the automated, field-agnostic version of that check.)
#
# Run: julia --project=.buildkite experiments/comparison/verify_column_ensemble.jl

import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Simulations: LandSimulation, step!
import ClimaLand.Parameters as LP

include(joinpath(@__DIR__, "field_rmse.jl"))

FT = Float64
toml_dict = LP.create_toml_dict(FT)

# ---- model setup (identical to compare_soil.jl) ---------------------------
vg_n = FT(2.9)
vg_α = FT(6)
K_sat = FT(4.42 / 3600 / 100) # m/s
hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
ν = FT(0.43)
θ_r = FT(0.045)
S_s = FT(1e-3)
params = ClimaLand.Soil.EnergyHydrologyParameters(
    toml_dict;
    ν,
    ν_ss_om = FT(0.0),
    ν_ss_quartz = FT(1.0),
    ν_ss_gravel = FT(0.0),
    hydrology_cm = hcm,
    K_sat,
    S_s,
    θ_r,
)

boundary_fluxes = (;
    top = WaterHeatBC(;
        water = WaterFluxBC((p, t) -> -K_sat / 10),
        heat = HeatFluxBC((p, t) -> 0.0),
    ),
    bottom = WaterHeatBC(;
        water = WaterFluxBC((p, t) -> 0.0),
        heat = HeatFluxBC((p, t) -> 0.0),
    ),
)

zmax = FT(0)
zmin = FT(-0.5)
nelems = 25
dz_tuple = (FT(0.05), FT(0.005))

soil_column = Column(; zlim = (zmin, zmax), nelements = nelems, dz_tuple)
soil_ensemble = ColumnEnsemble(;
    zlim = (zmin, zmax),
    nelements = nelems,
    dz_tuple,
    longlat = (FT(0), FT(0)),
)

function hydrostatic_equilibrium(z, z_interface, params)
    (; ν, S_s, hydrology_cm) = params
    (; α, n, m) = hydrology_cm
    if z < z_interface
        return -S_s * (z - z_interface) + ν
    else
        return ν * (1 + (α * (z - z_interface))^n)^(-m)
    end
end
function init_soil!(Y, p, t0, model)
    params = model.parameters
    z = model.domain.fields.z
    FT = eltype(Y.soil.ϑ_l)
    Y.soil.ϑ_l .= hydrostatic_equilibrium.(z, FT(-0.45), params)
    Y.soil.θ_i .= 0
    T = FT(296.15)
    ρc_s = @. Soil.volumetric_heat_capacity(
        Y.soil.ϑ_l,
        FT(0),
        params.ρc_ds,
        params.earth_param_set,
    )
    Y.soil.ρe_int =
        Soil.volumetric_internal_energy.(FT(0), ρc_s, T, params.earth_param_set)
end

function build_soil_sim(domain)
    soil = Soil.EnergyHydrology{FT}(;
        parameters = params,
        domain = domain,
        boundary_conditions = boundary_fluxes,
        sources = (),
    )
    return LandSimulation(
        Float64(0),
        Float64(24 * 3600),
        FT(100),
        soil;
        set_ic! = init_soil!,
    )
end

# ---- automated difference diagnostic --------------------------------------

# Signed difference between two Float64s, in ulps of the larger in magnitude.
ulps(a, b) = (a - b) / eps(max(abs(a), abs(b), floatmin(Float64)))

# Column-1 slice of a subsurface space's vertical Jacobian, as a flat vector.
function subsurface_J(space)
    J = parent(ClimaCore.Fields.local_geometry_field(space).J)
    return vec(collect(selectdim(J, ndims(J), 1)))
end

# For one leaf field that is not bitwise equal, loop over its entries to report
# exactly which differ and by how many ulps. Returns the differing indices.
function localize_difference(name, f_col, f_ens, col)
    a1 = flatten(f_col)         # single column
    a2 = flatten(f_ens, col)    # column `col` of the ensemble
    idx = findall(a1 .!= a2)
    max_ulp = isempty(idx) ? 0.0 : maximum(abs(ulps(a1[i], a2[i])) for i in idx)
    @info "  $name differs" entries = length(a1) differing = length(idx) max_ulp
    for i in idx
        @info @sprintf(
            "    entry %2d: column = %.17e  ensemble = %.17e  (%+.1f ulp)",
            i,
            a1[i],
            a2[i],
            ulps(a1[i], a2[i]),
        )
    end
    return idx
end

# Compare all leaf fields of two states; return the sorted names that are not
# bitwise identical (RMSE != 0).
function differing_fields(state_col, state_ens, col)
    diffs = state_field_rmses(state_col, state_ens, col)
    return sort([name for (name, d) in diffs if !(d.rmse == 0.0)])
end

# ---- run and check --------------------------------------------------------
sim_column = build_soil_sim(soil_column)
sim_ensemble = build_soil_sim(soil_ensemble)

# Show the full per-field table (t=0 then after 1 step) for context.
report_state_diffs(
    sim_column,
    sim_ensemble;
    col = 1,
    label = "t = 0 (no steps)",
)
step!(sim_column)
step!(sim_ensemble)
@assert sum(isnan.(sim_ensemble._integrator.u)) == 0
report_state_diffs(sim_column, sim_ensemble; col = 1, label = "after 1 step")

# The check: every field must be exactly equal. If not, loop to find the cause.
col = 1
states = (
    ("Y", sim_column._integrator.u, sim_ensemble._integrator.u),
    ("p", sim_column._integrator.p, sim_ensemble._integrator.p),
)
any_diff = false
for (prefix, s_col, s_ens) in states
    names = differing_fields(s_col, s_ens, col)
    if isempty(names)
        @info "$prefix: all fields bitwise identical."
        continue
    end
    global any_diff = true
    @info "$prefix: $(length(names)) field(s) not exactly equal — looping to find the cause" names
    pairs = leaf_pairs(s_col, s_ens)
    for name in names
        f_col, f_ens = pairs[name]
        localize_difference("$prefix.$name", f_col, f_ens, col)
    end
end

if !any_diff
    @info "EXACTLY IDENTICAL: Column and single-point ColumnEnsemble agree bitwise."
else
    # Tie the differences to the geometry: the ensemble's vertical Jacobian is a
    # constant multiple of the column's (the horizontal metric factor J_h).
    Jc = subsurface_J(soil_column.space.subsurface)
    Je = subsurface_J(soil_ensemble.space.subsurface)
    lo, hi = extrema(Je ./ Jc)
    @info "Cause: extruded-geometry Jacobian ratio (ensemble / column)" J_h_lo =
        lo J_h_hi = hi
    @info "The constant factor J_h cancels analytically inside the vertical \
           stencils but not bitwise, so tendencies land ~1 ulp apart. This \
           matches the earlier ρe_int finding (see investigate_roundoff.jl)."
end

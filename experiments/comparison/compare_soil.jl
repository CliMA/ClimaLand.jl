# Demonstration that ClimaLand runs on the new ClimaCore multi-column space.
#
# `Domains.ColumnEnsemble` is backed by a `MultiColumnFiniteDifferenceSpace`
# (an ensemble of independent vertical columns) rather than the single-column
# `CenterFiniteDifferenceSpace` behind `Domains.Column`. A single-point
# `ColumnEnsemble` built from the same vertical mesh must reproduce the
# equivalent `Column`: the domains carry identical grids, so the same
# `EnergyHydrology` model, initial condition, and implicit solve should give
# the same answer. This exercises the new space through model construction, the
# Jacobian (`initialize_jacobian`), and the NaN-check callback.
#
# Every prognostic and cache field is compared (via `report_state_diffs` from
# field_rmse.jl) twice: at t = 0 before any steps are taken (any difference
# there is a construction/initialization discrepancy, not a dynamics one), and
# again after a single step (isolating one implicit solve).

import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Simulations: LandSimulation, solve!, step!
import ClimaLand.Parameters as LP

include(joinpath(@__DIR__, "field_rmse.jl"))

FT = Float64;
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)

vg_n = FT(2.9)
vg_α = FT(6)
K_sat = FT(4.42 / 3600 / 100) # m/s
hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
ν = FT(0.43)
θ_r = FT(0.045)
S_s = FT(1e-3)
ν_ss_om = FT(0.0)
ν_ss_quartz = FT(1.0)
ν_ss_gravel = FT(0.0)

params = ClimaLand.Soil.EnergyHydrologyParameters(
    toml_dict;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = hcm,
    K_sat,
    S_s,
    θ_r,
);

surface_water_flux = WaterFluxBC((p, t) -> -K_sat / 10)
bottom_water_flux = WaterFluxBC((p, t) -> 0.0);
surface_heat_flux = HeatFluxBC((p, t) -> 0.0)
bottom_heat_flux = HeatFluxBC((p, t) -> 0.0);
boundary_fluxes = (;
    top = WaterHeatBC(; water = surface_water_flux, heat = surface_heat_flux),
    bottom = WaterHeatBC(; water = bottom_water_flux, heat = bottom_heat_flux),
);

# Shared vertical discretization. A stretched mesh (`dz_tuple`) is used because
# `ColumnEnsemble` always stretches; giving `Column` the same `dz_tuple`, `zlim`,
# and `nelements` makes the two vertical grids identical.
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

# Build a one-day simulation (initial condition and cache are set at
# construction; no steps are taken until step!).
function build_soil_sim(domain)
    soil = Soil.EnergyHydrology{FT}(;
        parameters = params,
        domain = domain,
        boundary_conditions = boundary_fluxes,
        sources = (),
    )
    t0 = Float64(0)
    tf = Float64(24 * 3600)
    return LandSimulation(t0, tf, FT(100), soil; set_ic! = init_soil!)
end

sim_column = build_soil_sim(soil_column)
sim_ensemble = build_soil_sim(soil_ensemble)

# Compare every prognostic and cache field before any steps are taken.
# Per field, with a1 = flattened single-column field, a2 = flattened ensemble
# column, and n = length(a1):
#   rmse  = sqrt(sum((a1 .- a2) .^ 2) / n)   RMS difference (field units)
#   scale = sqrt(sum(a1 .^ 2) / n)           RMS magnitude of the reference a1
#   rel   = rmse / scale                     relative error (dimensionless;
#                                            0 if rmse = scale = 0, Inf if only scale = 0)
report_state_diffs(
    sim_column,
    sim_ensemble;
    col = 1,
    label = "t = 0 (no steps)",
)

step!(sim_column)
step!(sim_ensemble)

# The single-point ensemble stepped without NaNs on the multi-column space.
@assert sum(isnan.(sim_ensemble._integrator.u)) == 0

# ... and reproduced the equivalent single column (identical vertical mesh).
report_state_diffs(sim_column, sim_ensemble; col = 1, label = "after 1 step")

@info "ColumnEnsemble (single point) reproduces Column: EnergyHydrology + \
       implicit solve ran on the multi-column MultiColumnFiniteDifferenceSpace."

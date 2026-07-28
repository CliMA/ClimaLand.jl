# Compare the single-column `Column` build against a single-point `ColumnEnsemble`
# (multi-column `MultiColumnFiniteDifferenceSpace`) for the standalone soil
# biogeochemistry experiment. The model setup is copied from `experiment.jl`: a
# `LandSoilBiogeochemistry` of
#   - soil: `EnergyHydrology` with a prescribed surface water flux, free drainage
#     at the bottom, and phase change,
#   - soil CO2 biogeochemistry (`SoilCO2Model`, prognostic microbes),
# driven by a `PrescribedAtmosphere` (only `atmos_p` matters here).
#
# Both builds share identical vertical mesh, parameters, forcing, and initial
# condition, so column 1 of the ensemble must reproduce the single column. Every
# prognostic and cache field is compared (via `report_state_diffs` from
# field_rmse.jl) at t = 0 before any steps, then again after a single step.
#
# Run: julia --project=.buildkite experiments/comparison/compare_biogeochemistry.jl

import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
import ClimaLand.Simulations: LandSimulation, step!
using Dates

import ClimaLand.Parameters as LP
import ClimaParams

include(joinpath(@__DIR__, "field_rmse.jl"))

const FT = Float64
toml_dict = LP.create_toml_dict(FT)

t0 = Float64(0)
tf = Float64(10000)
dt = Float64(10)

# Soil hydrology/energy parameters (shared across both builds).
ν = FT(0.556)
K_sat = FT(0.0443 / 3600 / 100) # m/s
S_s = FT(1e-3) # inverse meters
vg_n = FT(2.0)
vg_α = FT(2.6) # inverse meters
θ_r = FT(0.1)
ν_ss_om = FT(0.0)
ν_ss_quartz = FT(1.0)
ν_ss_gravel = FT(0.0)
soil_parameters = Soil.EnergyHydrologyParameters(
    toml_dict;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n),
    K_sat,
    S_s,
    θ_r,
)

zmax = FT(0)
zmin = FT(-1)
nelems = 20

# Prescribed atmosphere; identical for both builds (only atmos_p is used).
precipitation_function = (t) -> 1.0
snow_precip = (t) -> 1.0
atmos_T = (t) -> 1.0
atmos_u = (t) -> 1.0
atmos_q = (t) -> 1.0
atmos_p = (t) -> 100000.0
atmos_h = FT(30)
atmos_co2 = (t) -> 1.0
atmos = ClimaLand.PrescribedAtmosphere(
    TimeVaryingInput(precipitation_function),
    TimeVaryingInput(snow_precip),
    TimeVaryingInput(atmos_T),
    TimeVaryingInput(atmos_u),
    TimeVaryingInput(atmos_q),
    TimeVaryingInput(atmos_p),
    Dates.DateTime(2005),
    atmos_h,
    toml_dict;
    c_co2 = TimeVaryingInput(atmos_co2),
)

function init_soil!(Y, z, params)
    FT = eltype(Y.soil.ϑ_l)
    Y.soil.ϑ_l .= FT(0.33)
    Y.soil.θ_i .= FT(0.1)
    T = FT(279.85)
    ρc_s = Soil.volumetric_heat_capacity(
        FT(0.33),
        FT(0.1),
        params.ρc_ds,
        params.earth_param_set,
    )
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            FT(0.0),
            ρc_s,
            T,
            params.earth_param_set,
        )
end

function init_co2!(Y, z)
    function CO2_profile(z::FT) where {FT}
        C = FT(0.0)
        return FT(C)
    end
    function O2_profile(z::FT) where {FT}
        O2 = FT(0.21)
        return FT(O2)
    end
    function SOC_profile(z::FT) where {FT}
        Csom = FT(5.0)
        return FT(Csom)
    end
    Y.soilco2.CO2 .= CO2_profile.(z)
    Y.soilco2.O2 .= O2_profile.(z)
    Y.soilco2.SOC .= SOC_profile.(z)
end

function set_ic!(Y, p, t0, model)
    z = model.soil.domain.fields.z
    init_co2!(Y, z)
    init_soil!(Y, z, model.soil.parameters)
end

# Build the LandSoilBiogeochemistry model + simulation on the given domain,
# parameterized by `domain` so the setup is identical on a Column and a
# single-point ColumnEnsemble.
function build_biogeochem_sim(domain)
    top_flux_bc_w = Soil.WaterFluxBC((p, t) -> -0.00001)
    bot_flux_bc_w = Soil.FreeDrainage()
    top_flux_bc_h = Soil.HeatFluxBC((p, t) -> 0.0)
    bot_flux_bc_h = Soil.HeatFluxBC((p, t) -> 0.0)
    sources = (PhaseChange{FT}(),)
    boundary_conditions = (;
        top = WaterHeatBC(; water = top_flux_bc_w, heat = bot_flux_bc_h),
        bottom = WaterHeatBC(; water = bot_flux_bc_w, heat = bot_flux_bc_h),
    )
    soil = EnergyHydrology{FT}(;
        boundary_conditions,
        sources,
        domain,
        parameters = soil_parameters,
    )

    drivers =
        Soil.Biogeochemistry.SoilDrivers(PrognosticMet(soil_parameters), atmos)
    soilco2 = SoilCO2Model{FT}(domain, drivers, toml_dict)

    model = LandSoilBiogeochemistry{FT}(soil, soilco2)

    ode_algo = CTS.ExplicitAlgorithm(CTS.RK4())
    return LandSimulation(
        t0,
        tf,
        dt,
        model;
        set_ic!,
        timestepper = ode_algo,
        updateat = dt,
        diagnostics = nothing,
    )
end

column_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
ensemble_domain = ColumnEnsemble(;
    zlim = (zmin, zmax),
    nelements = nelems,
    longlat = (FT(0), FT(0)),
)

sim_column = build_biogeochem_sim(column_domain)
sim_ensemble = build_biogeochem_sim(ensemble_domain)

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

for i in 1:500
    step!(sim_column)
    step!(sim_ensemble)
end

@assert sum(isnan.(sim_ensemble._integrator.u)) == 0

report_state_diffs(sim_column, sim_ensemble; col = 1, label = "after 500 step")

@info "Standalone soil biogeochemistry LandSoilBiogeochemistry ran on the \
       multi-column MultiColumnFiniteDifferenceSpace and reproduced the single \
       Column."

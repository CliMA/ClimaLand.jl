# Compare the single-column `Column` build against a single-point `ColumnEnsemble`
# (multi-column `MultiColumnFiniteDifferenceSpace`) for the standalone Richards
# equation experiment. The model setup is copied verbatim from
# `richards_comparison.jl` (the Bonan clay and sand parameter sets); the only
# thing that differs between the two builds is the domain, so column 1 of the
# ensemble must reproduce the single column.
#
# Every prognostic and cache field is compared (via `report_state_diffs` from
# field_rmse.jl) at t = 0 before any steps, then again after a single step.
#
# Run: julia --project=.buildkite experiments/comparison/compare_richards.jl

import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
using ClimaLand.Soil
import ClimaLand.Simulations: LandSimulation, step!

import ClimaLand.Parameters as LP

include(joinpath(@__DIR__, "field_rmse.jl"))

const FT = Float64
toml_dict = LP.create_toml_dict(FT)

t0 = Float64(0)

# Compare a Column build against a single-point ColumnEnsemble build of the same
# model. Per field, with a1 = flattened single-column field, a2 = flattened
# ensemble column, and n = length(a1):
#   rmse  = sqrt(sum((a1 .- a2) .^ 2) / n)   RMS difference (field units)
#   scale = sqrt(sum(a1 .^ 2) / n)           RMS magnitude of the reference a1
#   rel   = rmse / scale                     relative error (dimensionless;
#                                            0 if rmse = scale = 0, Inf if only scale = 0)
function compare_spaces(name, sim_column, sim_ensemble)
    report_state_diffs(
        sim_column,
        sim_ensemble;
        col = 1,
        label = "$name: t = 0 (no steps)",
    )
    step!(sim_column)
    step!(sim_ensemble)
    @assert sum(isnan.(sim_ensemble._integrator.u)) == 0
    report_state_diffs(
        sim_column,
        sim_ensemble;
        col = 1,
        label = "$name: after 1 step",
    )
end

# ---- clay (from the "Richards comparison to Bonan; clay" testset) ----------
let
    dt = Float64(1e3)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(1.43)
    vg_α = FT(2.6) # inverse meters
    hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
    θ_r = FT(0.124)
    zmax = FT(0)
    zmin = FT(-1.5)
    nelems = 150

    top_state_bc = MoistureStateBC((p, t) -> ν - 1e-3)
    bot_flux_bc = FreeDrainage()
    sources = ()
    boundary_states = (; top = top_state_bc, bottom = bot_flux_bc)
    params = Soil.RichardsParameters(;
        ν = ν,
        hydrology_cm = hcm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )

    function set_ic!(Y, p, t0, model)
        Y.soil.ϑ_l .= FT(0.24)
    end

    stepper = CTS.ARS111()
    norm_condition = CTS.MaximumError(FT(1e-8))
    conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 50,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
            convergence_checker = conv_checker,
        ),
    )

    build(domain) = LandSimulation(
        t0,
        2 * dt,
        dt,
        Soil.RichardsModel{FT}(;
            parameters = params,
            domain = domain,
            boundary_conditions = boundary_states,
            sources = sources,
        );
        set_ic!,
        timestepper = ode_algo,
        updateat = dt,
        diagnostics = nothing,
    )

    compare_spaces(
        "clay",
        build(Column(; zlim = (zmin, zmax), nelements = nelems)),
        build(
            ColumnEnsemble(;
                zlim = (zmin, zmax),
                nelements = nelems,
                longlat = (FT(0), FT(0)),
            ),
        ),
    )
end

# ---- sand (from the "Richards comparison to Bonan; sand" testset) ----------
let
    dt = Float64(1)
    ν = FT(0.287)
    K_sat = FT(34 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(3.96)
    vg_α = FT(2.7) # inverse meters
    vg_m = FT(1)
    S_c = ν + 1
    hcm = vanGenuchten(vg_α, vg_n, vg_m, S_c)
    θ_r = FT(0.075)
    zmax = FT(0)
    zmin = FT(-1.5)
    nelems = 150

    top_state_bc = MoistureStateBC((p, t) -> 0.267)
    bot_flux_bc = FreeDrainage()
    sources = ()
    boundary_states = (; top = top_state_bc, bottom = bot_flux_bc)
    params = Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)

    function set_ic!(Y, p, t0, model)
        Y.soil.ϑ_l .= FT(0.1)
    end

    stepper = CTS.ARS111()
    norm_condition = CTS.MaximumError(FT(1e-8))
    conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 50,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
            convergence_checker = conv_checker,
        ),
    )

    build(domain) = LandSimulation(
        t0,
        2 * dt,
        dt,
        Soil.RichardsModel{FT}(;
            parameters = params,
            domain = domain,
            boundary_conditions = boundary_states,
            sources = sources,
        );
        set_ic!,
        timestepper = ode_algo,
        updateat = dt,
        diagnostics = nothing,
    )

    compare_spaces(
        "sand",
        build(Column(; zlim = (zmin, zmax), nelements = nelems)),
        build(
            ColumnEnsemble(;
                zlim = (zmin, zmax),
                nelements = nelems,
                longlat = (FT(0), FT(0)),
            ),
        ),
    )
end

@info "Standalone Richards RichardsModel ran on the multi-column \
       MultiColumnFiniteDifferenceSpace and reproduced the single Column for \
       both the Bonan clay and sand parameter sets."

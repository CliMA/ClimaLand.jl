# ══════════════════════════════════════════════════════════════════════════════════
# DK-Sor Site-Level Calibration  —  Production Version  (parallel, SLURM)
#
# Identical EKP setup to calibrate_DK-Sor_tutorial.jl but parallelises the
# ensemble forward-model evaluations across Julia workers using Distributed.pmap.
#
# Usage (submit via companion SLURM script):
#   sbatch experiments/integrated/fluxnet/slurm_calibrate_dksor.sh
#
# Worker count is read from $SLURM_CPUS_PER_TASK (set in the SLURM script).
# Each worker independently builds and integrates a DK-Sor model for the
# assigned parameter vector and minibatch years.
#
# Per-iteration wall time estimate:
#   ~500 s/model-year × 2 years × 1 (parallel over 23 workers) ≈ 17 min/iter
#   10 iterations + ~10 min compilation  ≈  3 hours total
# ══════════════════════════════════════════════════════════════════════════════════

# ── Distributed setup (must come before @everywhere imports) ──────────────────────
using Distributed

# Respect SLURM CPU allocation; fall back to all available hardware threads - 1.
const N_WORKERS = max(1,
    parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", string(Sys.CPU_THREADS))) - 1)

addprocs(N_WORKERS;
    exeflags = "--project=$(Base.active_project())",
    topology  = :master_worker,
)

@info "Started $(nworkers()) Julia workers (requested $(N_WORKERS))"

# ── ClimaComms backend must be loaded in two separate @everywhere calls because
#    Julia lowers (macro-expands) an entire @everywhere begin...end block before
#    evaluating any import statements.  If @import_required_backends is in the
#    same block as `import ClimaComms`, workers fail with
#    "UndefVarError: ClimaComms not defined".
@everywhere import ClimaComms
@everywhere ClimaComms.@import_required_backends

# ══════════════════════════════════════════════════════════════════════════════════
#  Everything inside @everywhere is loaded on all workers AND the master process.
# ══════════════════════════════════════════════════════════════════════════════════
@everywhere begin

using ClimaCore
import ClimaParams as CP
using Dates
using ClimaDiagnostics
using ClimaUtilities
using Logging

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.TimeManager: date
using NCDatasets

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using LinearAlgebra
using Statistics
import Random
using JLD2
using CairoMakie

# ── Configuration (replicated on every worker) ────────────────────────────────────
const FT          = Float64
const SITE_ID_VAL = :DK_Sor
const DT          = Float64(450)
const CLIMALAND_DIR = pkgdir(ClimaLand)
const MET_NC_PATH   = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
const FLUX_NC_PATH  = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")

const _SITE_LOC    = FluxnetSimulations.get_location(FT, Val(SITE_ID_VAL))
const _SITE_HEIGHT = FluxnetSimulations.get_fluxtower_height(FT, Val(SITE_ID_VAL))
const _SITE_PARAMS = FluxnetSimulations.get_parameters(FT, Val(SITE_ID_VAL))

const time_offset = _SITE_LOC.time_offset
const lat         = _SITE_LOC.lat
const long        = _SITE_LOC.long
const atmos_h     = _SITE_HEIGHT.atmos_h

const CALIB_YEARS    = collect(2003:2013)
const MINIBATCH_SIZE = 2
const N_ITERATIONS   = 10
const RNG_SEED       = 42

const MIN_VALID_FRAC = 0.5
const SIGMA2_LHF     = 400.0
const SIGMA2_SHF     = 625.0
const SIGMA2_NEE     = 0.25
const SIGMA2_MISS    = 1.0e6

const PARAM_NAMES = [
    "pmodel_cstar",
    "pmodel_β_c3",
    "pmodel_α",
    "moisture_stress_c",
    "soilCO2_reference_rate",
    "soilCO2_activation_energy",
    "michaelis_constant",
    "O2_michaelis_constant",
    "root_leaf_nitrogen_ratio",
    "relative_contribution_factor",
    "emissivity_bare_soil",
]
const N_PARAMS = length(PARAM_NAMES)
const N_ENS    = 30   # EKI ensemble size (same order as global calibration)

# ── Helper: write parameter override TOML ────────────────────────────────────────
function write_override_toml(path, param_names, param_values)
    open(path, "w") do io
        for (name, val) in zip(param_names, Float64.(param_values))
            println(io, "[\"$(name)\"]")
            println(io, "value = $(val)")
            println(io, "type  = \"float\"")
            println(io)
        end
    end
end

# ── Helper: extract monthly means from DictWriter diagnostics ─────────────────────
function extract_monthly_means(simulation, start_date)
    writer    = simulation.diagnostics[1].output_writer
    diag_keys = ["lhf_30m_average", "shf_30m_average", "nee_30m_average"]
    convs     = [1.0, 1.0, 12.0 * 86400.0]   # NEE: mol/m²/s → gC/m²/d
    result    = Float64[]
    for (key, conv) in zip(diag_keys, convs)
        times, data = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, key)
        dates = [date(t) for t in times]
        for month in 1:12
            mask = Dates.month.(dates) .== month
            push!(result, mean(data[mask]) * conv)
        end
    end
    return result
end

# ── Core model function (runs on workers) ────────────────────────────────────────
function run_model_year(params_vec, year)
    (; soil_ν, soil_K_sat, soil_S_s, soil_vg_n, soil_vg_α, θ_r,
       ν_ss_quartz, ν_ss_om, ν_ss_gravel,
       z_0m_soil, z_0b_soil, soil_α_PAR, soil_α_NIR,
       Ω, χl, α_PAR_leaf, λ_γ_PAR, τ_PAR_leaf,
       α_NIR_leaf, τ_NIR_leaf, ϵ_canopy, ac_canopy,
       Drel, g0, SAI, f_root_to_shoot,
       conductivity_model, retention_model,
       plant_ν, plant_S_s, rooting_depth, h_canopy) = _SITE_PARAMS

    start_date = DateTime(year,     1, 1)
    stop_date  = DateTime(year + 1, 1, 1)

    # Guard: TransformUnscented sigma points near prior bounds can map to NaN
    # in constrained space. Return NaN output so EKP update ignores this member.
    if any(isnan, params_vec)
        @warn "run_model_year($year): NaN in params_vec, skipping"
        return fill(NaN, 36)
    end

    tmpfile = tempname() * ".toml"
    write_override_toml(tmpfile, PARAM_NAMES, params_vec)

    try
        local_toml = LP.create_toml_dict(FT; override_files = [tmpfile])
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

        (; dz_tuple, nelements, zmin, zmax) =
            FluxnetSimulations.get_domain_info(FT, Val(SITE_ID_VAL))
        land_domain    = Column(; zlim = (zmin, zmax), nelements, dz_tuple,
                                  longlat = (long, lat))
        surface_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

        (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
            MET_NC_PATH, lat, long, time_offset, atmos_h,
            start_date, local_toml, FT,
        )

        LAI, maxLAI = NCDataset(MET_NC_PATH, "r") do ds
            time_vals = ds["time"][:]
            lai_data  = Float64.(coalesce.(ds["LAI_alternative"][1, 1, :], NaN))
            lai_secs  = Float64[
                Second(t - Hour(time_offset) - start_date).value
                for t in time_vals
            ]
            valid = .!isnan.(lai_data)
            TimeVaryingInput(lai_secs[valid], lai_data[valid]), maximum(lai_data[valid])
        end
        RAI = maxLAI * f_root_to_shoot

        forcing     = (; atmos, radiation)
        soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
                          PAR_albedo = soil_α_PAR, NIR_albedo = soil_α_NIR)
        retention_parameters = (;
            ν = soil_ν, θ_r, K_sat = soil_K_sat,
            hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
        )
        composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

        # emissivity not passed explicitly → read from local_toml
        soil = Soil.EnergyHydrology{FT}(
            land_domain, forcing, local_toml;
            prognostic_land_components,
            additional_sources  = (ClimaLand.RootExtraction{FT}(),),
            albedo              = soil_albedo,
            runoff              = ClimaLand.Soil.Runoff.SurfaceRunoff(),
            retention_parameters, composition_parameters,
            S_s = soil_S_s, z_0m = z_0m_soil, z_0b = z_0b_soil,
        )

        co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
        drivers  = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
        soilco2  = Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, drivers, local_toml)

        radiation_parameters = (;
            Ω, G_Function = CLMGFunction(χl),
            α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
        )
        radiative_transfer   = Canopy.TwoStreamModel{FT}(
            surface_domain, local_toml; radiation_parameters, ϵ_canopy)
        photosynthesis       = Canopy.PModel{FT}(surface_domain, local_toml)
        conductance          = Canopy.PModelConductance{FT}(local_toml; Drel)
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, local_toml; soil_params = (; ν = soil_ν, θ_r))
        hydraulics = Canopy.PlantHydraulicsModel{FT}(
            surface_domain, local_toml;
            ν = plant_ν, S_s = plant_S_s, conductivity_model, retention_model)
        biomass = Canopy.PrescribedBiomassModel{FT}(;
            LAI, SAI, RAI, rooting_depth, height = h_canopy)
        energy  = Canopy.BigLeafEnergyModel{FT}(local_toml; ac_canopy)
        ground  = ClimaLand.PrognosticGroundConditions{FT}()

        canopy = Canopy.CanopyModel{FT}(
            surface_domain, (; atmos, radiation, ground), LAI, local_toml;
            prognostic_land_components,
            radiative_transfer, photosynthesis, conductance,
            soil_moisture_stress, hydraulics, energy, biomass,
        )
        snow = Snow.SnowModel(FT, surface_domain, forcing, local_toml, DT;
            prognostic_land_components)
        land = LandModel{FT}(canopy, snow, soil, soilco2, nothing)

        T_init_K = NCDataset(MET_NC_PATH, "r") do ds
            Float64(coalesce(ds["Tair"][1, 1, 1], 283.15))
        end
        function set_ic!(Y, p, t, model)
            FT_l = eltype(Y.soil.ρe_int)
            θ_l_0 = soil.parameters.θ_r .+
                    (soil.parameters.ν .- soil.parameters.θ_r) .* FT_l(0.95)
            Y.soil.ϑ_l .= θ_l_0
            Y.soil.θ_i .= FT_l(0)
            ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l, Y.soil.θ_i,
                soil.parameters.ρc_ds, soil.parameters.earth_param_set)
            Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
                Y.soil.θ_i, ρc_s, FT_l(T_init_K),
                soil.parameters.earth_param_set)
            Y.snow.S .= FT_l(0); Y.snow.S_l .= FT_l(0); Y.snow.U .= FT_l(0)
            Y.canopy.energy.T .= FT_l(T_init_K)
            ψ_leaf_0 = FT_l(-2e5 / 9800)
            hyd = model.canopy.hydraulics
            S_l_ini = ClimaLand.Canopy.inverse_water_retention_curve.(
                hyd.parameters.retention_model, ψ_leaf_0,
                hyd.parameters.ν, hyd.parameters.S_s)
            Y.canopy.hydraulics.ϑ_l .= ClimaLand.Canopy.augmented_liquid_fraction.(
                hyd.parameters.ν, S_l_ini)
            Y.soilco2.CO2 .= FT_l(6e-5)
            Y.soilco2.O2  .= FT_l(0.08)
            Y.soilco2.SOC .= FT_l(5.0)
        end

        diags = ClimaLand.default_diagnostics(
            land, start_date;
            output_writer    = ClimaDiagnostics.Writers.DictWriter(),
            output_vars      = ["lhf", "shf", "nee"],
            reduction_period = :halfhourly,
        )
        simulation = LandSimulation(
            start_date, stop_date, DT, land;
            set_ic!, updateat = DT, diagnostics = diags,
        )
        Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
            solve!(simulation)
        end
        return extract_monthly_means(simulation, start_date)

    catch e
        @warn "run_model_year($year) failed on worker $(myid()): $e"
        return fill(NaN, 36)
    finally
        rm(tmpfile; force = true)
    end
end

# G function: run all minibatch years, return concatenated monthly means
function G(params_vec, years)
    return vcat([run_model_year(params_vec, yr) for yr in years]...)
end

end  # @everywhere

# ══════════════════════════════════════════════════════════════════════════════════
#  Master-only code: EKP setup, checkpoint, calibration loop, plots
# ══════════════════════════════════════════════════════════════════════════════════

const OUTDIR          = joinpath(CLIMALAND_DIR,
    "experiments/integrated/fluxnet/DK-Sor/out/calibration")
const CHECKPOINT_FILE = joinpath(OUTDIR, "ekp_checkpoint.jld2")
mkpath(OUTDIR)

function get_calibration_prior()
    priors = [
        PD.constrained_gaussian("pmodel_cstar",                 0.448094,   0.05,    0.2,    0.7),
        PD.constrained_gaussian("pmodel_β_c3",                  68.2139,    40.0,   10.0,  300.0),
        PD.constrained_gaussian("pmodel_α",                     0.964115,   0.02,   0.85,   0.999),
        PD.constrained_gaussian("moisture_stress_c",            0.377838,   0.15,   0.05,   1.0),
        PD.constrained_gaussian("soilCO2_reference_rate",       2.17564e-7, 2.0e-7, 1.0e-9, 2.0e-6),
        PD.constrained_gaussian("soilCO2_activation_energy",    37357.0,    25000., 5000., 150000.),
        PD.constrained_gaussian("michaelis_constant",           0.459997,   0.25,   0.001,  2.0),
        PD.constrained_gaussian("O2_michaelis_constant",        0.00255006, 0.002,  1.0e-5, 0.1),
        PD.constrained_gaussian("root_leaf_nitrogen_ratio",     0.0625785,  0.04,   0.01,   3.0),
        PD.constrained_gaussian("relative_contribution_factor", 0.844417,   0.15,   0.0,    1.5),
        PD.constrained_gaussian("emissivity_bare_soil",         0.96,       0.03,   0.6,    1.0),
    ]
    return PD.combine_distributions(priors)
end

function load_obs_per_year(flux_nc_path, calib_years)
    obs_vec = EKP.Observation[]
    NCDataset(flux_nc_path, "r") do ds
        flux_times = DateTime.(ds["time"][:])
        lhf_all    = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
        shf_all    = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
        nee_all    = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
        for yr in calib_years
            yr_mask  = Dates.year.(flux_times) .== yr
            lhf_yr   = lhf_all[yr_mask]
            shf_yr   = shf_all[yr_mask]
            nee_yr   = nee_all[yr_mask]
            times_yr = flux_times[yr_mask]
            mean_vec  = Float64[]
            noise_vec = Float64[]
            for (data_yr, σ2) in [(lhf_yr, SIGMA2_LHF),
                                  (shf_yr, SIGMA2_SHF),
                                  (nee_yr, SIGMA2_NEE)]
                for month in 1:12
                    m_mask = Dates.month.(times_yr) .== month
                    vals   = data_yr[m_mask]
                    valid  = filter(v -> !isnan(v) && v < 1e19, vals)
                    frac   = length(valid) / max(sum(m_mask), 1)
                    if frac >= MIN_VALID_FRAC
                        push!(mean_vec,  mean(valid))
                        push!(noise_vec, σ2)
                    else
                        push!(mean_vec,  0.0)
                        push!(noise_vec, SIGMA2_MISS)
                    end
                end
            end
            push!(obs_vec,
                EKP.Observation(mean_vec, Diagonal(noise_vec), "DK_Sor_$(yr)"))
        end
    end
    return obs_vec
end

function save_checkpoint(ekp, iteration)
    JLD2.jldsave(CHECKPOINT_FILE; ekp, iteration)
    @info "Checkpoint saved at iteration $iteration"
end

function load_checkpoint()
    isfile(CHECKPOINT_FILE) || return nothing, 0
    d = JLD2.load(CHECKPOINT_FILE)
    @info "Resuming from checkpoint: iteration $(d["iteration"])"
    return d["ekp"], d["iteration"]
end

function plot_parameter_evolution(ekp, prior)
    fig = Figure(size = ((N_PARAMS + 1) * 320, 380))
    for i in 1:N_PARAMS
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, N_PARAMS + 1], ekp)
    outpath = joinpath(OUTDIR, "calibration_parameter_evolution.png")
    save(outpath, fig)
    @info "Parameter evolution → $(outpath)"
end

function plot_daily_comparison(ϕ_final; year = 2009)
    @info "Running posterior model for $year (daily comparison)…"
    monthly_means = run_model_year(ϕ_final, year)
    lhf_model_monthly = monthly_means[1:12]
    shf_model_monthly = monthly_means[13:24]
    nee_model_monthly = monthly_means[25:36]

    lhf_obs_daily = Float64[]
    shf_obs_daily = Float64[]
    nee_obs_daily = Float64[]
    obs_dates     = Date[]

    NCDataset(FLUX_NC_PATH, "r") do ds
        flux_times = DateTime.(ds["time"][:])
        lhf_all    = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
        shf_all    = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
        nee_all    = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
        mask = Dates.year.(flux_times) .== year
        append!(obs_dates,     [Date(t) for t in flux_times[mask]])
        append!(lhf_obs_daily, lhf_all[mask])
        append!(shf_obs_daily, shf_all[mask])
        append!(nee_obs_daily, nee_all[mask])
    end
    for arr in (lhf_obs_daily, shf_obs_daily, nee_obs_daily)
        arr[arr .>= 1e19] .= NaN
    end

    month_centers = [Date(year, m, 15) for m in 1:12]
    fig = Figure(size = (1400, 1000))
    vars = [
        (lhf_obs_daily, lhf_model_monthly, "LHF", "W/m²"),
        (shf_obs_daily, shf_model_monthly, "SHF", "W/m²"),
        (nee_obs_daily, nee_model_monthly, "NEE", "gC/m²/d"),
    ]
    for (row, (obs_daily, model_monthly, label, unit)) in enumerate(vars)
        ax = Axis(fig[row, 1];
            title  = "$(label) — $(year)  (posterior)",
            xlabel = "Date",
            ylabel = "$(label) ($(unit))",
        )
        valid = .!isnan.(obs_daily)
        scatter!(ax, obs_dates[valid], obs_daily[valid];
            color = :red, markersize = 4, label = "Obs (daily)")
        lines!(ax, month_centers, model_monthly;
            color = :blue, linewidth = 2, label = "Model monthly mean")
        axislegend(ax; position = :lt)
    end
    outpath = joinpath(OUTDIR, "calibration_daily_comparison_$(year).png")
    save(outpath, fig)
    @info "Daily comparison → $(outpath)"
end

# ── Main ─────────────────────────────────────────────────────────────────────────
function main()
    rng   = Random.MersenneTwister(RNG_SEED)
    prior = get_calibration_prior()

    @info "DK-Sor parallel calibration | workers=$(nworkers()) | N_ens=$(N_ENS)"

    obs_per_year = load_obs_per_year(FLUX_NC_PATH, CALIB_YEARS)

    ekp, start_iter = load_checkpoint()
    if isnothing(ekp)
        minibatcher = EKP.RandomFixedSizeMinibatcher(MINIBATCH_SIZE)
        obs_series  = EKP.ObservationSeries(
            obs_per_year, minibatcher, string.(CALIB_YEARS))
        # Use Inversion() (standard EKI) — TransformUnscented sigma points near
        # prior bounds produce NaN in constrained space → invalid TOML.
        # Inversion() requires an explicit initial ensemble sampled from the prior.
        initial_ensemble = EKP.construct_initial_ensemble(rng, prior, N_ENS)
        ekp = EKP.EnsembleKalmanProcess(
            initial_ensemble,
            obs_series,
            EKP.Inversion();
            verbose   = true,
            rng,
            scheduler = EKP.DataMisfitController(terminate_at = 100),
        )
    end

    N_ens = EKP.get_N_ens(ekp)

    for i in (start_iter + 1):N_ITERATIONS
        @info "══ Iteration $(i) / $(N_ITERATIONS) ══"
        mb       = EKP.get_current_minibatch(ekp)
        mb_years = CALIB_YEARS[mb]
        @info "  Minibatch years: $(mb_years)"

        params_i = EKP.get_ϕ_final(prior, ekp)   # N_PARAMS × N_ens

        # pmap distributes members across workers; each worker runs its years serially.
        t_wall = @elapsed begin
            G_results = pmap(1:N_ens) do j
                G(params_i[:, j], mb_years)
            end
            G_ens = hcat(G_results...)
        end
        @info "  Parallel forward models done in $(round(t_wall/60, digits=1)) min"

        EKP.update_ensemble!(ekp, G_ens)
        save_checkpoint(ekp, i)

        ϕ_mean = EKP.get_ϕ_mean_final(prior, ekp)
        @info "  Posterior means:"
        for (name, val) in zip(PARAM_NAMES, ϕ_mean)
            @info "    $(rpad(name, 35)) = $(round(val, sigdigits = 4))"
        end
    end

    ϕ_final = EKP.get_ϕ_mean_final(prior, ekp)
    @info "══ Calibration complete ══"
    for (name, val) in zip(PARAM_NAMES, ϕ_final)
        @info "  $(rpad(name, 35)) = $(val)"
    end

    plot_parameter_evolution(ekp, prior)
    plot_daily_comparison(ϕ_final; year = 2009)
end

main()

"""
Per-site CalLMIP Phase 1b deliverable runs (prior/posterior, Cal + Val sites).

CONSISTENCY REQUIREMENT: uses the exact model of the calibration engine
(default LandModel, MODIS LAI, 1-site ColumnEnsemble via the same forcing
helpers), so the calibrated parameters are applied to the model they were
calibrated on.

Protocol handling per site (native met window y0–y1):
  - spin-up: min(5, window−1) years cycling the site's own met from default
    ICs, then the FULL prognostic state is handed to the production run
    (parent-copy per leaf; D5 documented deviation from 1850 equilibrium).
  - production run: (y0-01-02 → y1+1-01-02) UTC. Day 1 of the output axis
    (Jan 1, y0) is fill (the UTC start must lie inside every offset's native
    record); the axis tail is padded by year-recycling so the final daily
    window closes and Dec 31, y1 is real.
  - daily CalLMIP diagnostics via DictWriter → Phase 1a-format JLD2 →
    write_callmip_nc_phase1b.

Library: `run_site_deliverable(site_id; params_vec, stage)`.
CLI:  julia --project=.buildkite deliverable_runs_phase1b.jl <site_id> <Prior|Posterior>
      (Posterior reads output_calibration/posterior_ensemble.jld2 ϕ_mean)
"""

using Dates
using Statistics
using JLD2
import ClimaCore
import ClimaDiagnostics
import ClimaUtilities.TimeManager: date

include(joinpath(@__DIR__, "forward_model_phase1b.jl"))
include(joinpath(@__DIR__, "write_callmip_netcdf_phase1b.jl"))

const CALLMIP_SURFACE_VARS_1B = [
    "nee",
    "lhf",
    "shf",
    "gpp",
    "er",
    "hr",
    "trans",
    "ct",
    "lai",
    "cveg",
    "soilrn",
    "soillhf",
    "soilshf",
]
const CALLMIP_COLUMN_VARS_1B = ["swc", "tsoil", "soc"]
const DELIV_OUTDIR =
    get(ENV, "CALLMIP1B_DELIV_DIR", joinpath(@__DIR__, "output_deliverables"))

copy_state!(dst::ClimaCore.Fields.Field, src::ClimaCore.Fields.Field) =
    (parent(dst) .= parent(src))
function copy_state!(dst, src)
    for p in propertynames(dst)
        copy_state!(getproperty(dst, p), getproperty(src, p))
    end
end

function _build_sim(
    site,
    start_date,
    stop_date,
    toml_dict;
    diagnostics,
    set_ic!,
)
    domain = ColumnEnsemble(;
        zlim = FT1B.((-2, 0)),
        nelements = 10,
        longlat = [FT1B.((site.long, site.lat))],
    )
    surface_space = domain.space.surface
    forcing = prescribed_forcing_callmip(
        [site],
        surface_space,
        start_date,
        toml_dict,
        FT1B,
    )
    LAI = prescribed_lai_callmip([site], surface_space, start_date)
    land = Base.invokelatest(
        LandModel{FT1B},
        forcing,
        LAI,
        toml_dict,
        domain,
        DT1B,
    )
    diags =
        isnothing(diagnostics) ? () :
        vcat(
            [
                ClimaLand.default_diagnostics(
                    land,
                    start_date,
                    "";
                    output_writer = diagnostics,
                    output_vars = vars,
                    reduction_period = :daily,
                ) for vars in (CALLMIP_SURFACE_VARS_1B, CALLMIP_COLUMN_VARS_1B)
            ]...,
        )
    kwargs = isnothing(set_ic!) ? (;) : (; set_ic!)
    sim = Base.invokelatest(
        LandSimulation,
        start_date,
        stop_date,
        DT1B,
        land;
        updateat = Second(DT1B),
        diagnostics = diags,
        kwargs...,
    )
    return sim
end

"""
    run_site_deliverable(site_id; params_vec = Float64[], stage = "Prior",
                         spinup_years = 5, outdir = DELIV_OUTDIR)

Spin up, run the native window, extract daily CalLMIP diagnostics, and write
both the JLD2 and the final NetCDF. Returns the NetCDF path.
"""
function run_site_deliverable(
    site_id;
    params_vec = Float64[],
    stage = "Prior",
    spinup_years = 5,
    outdir = DELIV_OUTDIR,
)
    mkpath(outdir)
    info = site_info(site_id)
    y0, y1 = info.met_years
    nspin = min(spinup_years, y1 - y0)
    toml_dict, tmpfile = toml_with_overrides(params_vec)
    try
        raw = load_site(site_id)
        # pad only the tail sliver so the final daily window can close
        axis = collect(first(raw.utc_dates):Minute(30):DateTime(y1 + 1, 1, 2))
        site = pad_to_axis(raw, axis)
        start_date = DateTime(y0, 1, 2)

        @info "[$site_id $stage] spin-up $(nspin) yr"
        sim_spin = _build_sim(
            site,
            start_date,
            DateTime(y0 + nspin, 1, 1),
            toml_dict;
            diagnostics = nothing,
            set_ic! = nothing,
        )
        solve!(sim_spin)
        Y_spun = sim_spin._integrator.u

        @info "[$site_id $stage] production run $(y0)–$(y1)"
        writer = ClimaDiagnostics.Writers.DictWriter()
        sim = _build_sim(
            site,
            start_date,
            DateTime(y1 + 1, 1, 1),
            toml_dict;
            diagnostics = writer,
            set_ic! = (Y, p, t, model) -> copy_state!(Y, Y_spun),
        )
        solve!(sim)

        # extraction (single column) in the Phase 1a JLD2 layout
        function series(var)
            key = "$(var)_1d_average"
            times = sort!(collect(keys(writer[key])))
            vals = [vec(Array(parent(writer[key][t])))[1] for t in times]
            return Date.(date.(times)), vals
        end
        dates, _ = series(CALLMIP_SURFACE_VARS_1B[1])
        surface_data = Dict{String, Vector{Float64}}(
            v => series(v)[2] for v in CALLMIP_SURFACE_VARS_1B
        )
        column_data = Dict{String, Matrix{Float64}}()
        z_soil = Float64[]
        for v in CALLMIP_COLUMN_VARS_1B
            key = "$(v)_1d_average"
            times = sort!(collect(keys(writer[key])))
            cols = [vec(Array(parent(writer[key][t]))) for t in times]
            column_data[v] = reduce(hcat, cols)
            if isempty(z_soil)
                z_soil = vec(
                    Array(
                        parent(
                            ClimaCore.Fields.coordinate_field(
                                axes(writer[key][times[1]]),
                            ).z,
                        ),
                    ),
                )
            end
        end

        diag_path = joinpath(outdir, "$(site_id)_$(stage).jld2")
        jldsave(diag_path; dates, surface_data, column_data, z_soil)

        cal_val = site_id in CALIBRATION_SITE_IDS ? "Calibration" : "Validation"
        tag = cal_val == "Calibration" ? "Cal" : "Val"
        nc_path = joinpath(
            outdir,
            "ClimaLand.v1_Phase1b_Scen1_$(site_id)_$(tag)_$(stage).nc",
        )
        write_callmip_nc_phase1b(
            diag_path,
            nc_path;
            site_id,
            lat = info.lat,
            long = info.long,
            met_years = info.met_years,
            stage,
            cal_val,
        )
        @info "[$site_id $stage] wrote $nc_path"
        return nc_path
    finally
        isnothing(tmpfile) || rm(tmpfile; force = true)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    site_id = ARGS[1]
    stage = length(ARGS) >= 2 ? ARGS[2] : "Prior"
    params = Float64[]
    if stage == "Posterior"
        post = JLD2.load(
            joinpath(@__DIR__, "output_calibration", "posterior_ensemble.jld2"),
        )
        params = Float64.(post["ϕ_mean"])
        @assert post["param_names"] == PARAM_NAMES
    end
    run_site_deliverable(site_id; params_vec = params, stage)
end

"""
ClimaCalibrate forward model for CalLMIP Phase 1a full diagnostic simulations.

Runs a single-column ClimaLand simulation for the full 1997-01-01 – 2014-12-31
period (with 2-year spinup starting 1995-01-01) and saves all CalLMIP-required
diagnostic variables to JLD2.

Called by run_callmip_simulations.jl (not by EKI).

Output variables saved:
  Surface (scalar daily):
    nee, lhf, shf, gpp, er, trans, ct, lai, cveg
  Column (profile daily → integrated in write_callmip_netcdf.jl):
    swc, tsoil, soc

sign convention: nee positive = source (mol CO2/m²/s).
write_callmip_netcdf.jl flips sign and converts to kgC/m²/s for `nep`.
"""

import ClimaCalibrate
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date

import EnsembleKalmanProcesses as EKP
import JLD2
using Dates
using NCDatasets
using Statistics

# Include shared helpers (forcing loader, model builder, IC setter)
include(joinpath(@__DIR__, "model_interface.jl"))

const _CALLMIP_SURFACE_VARS = [
    "nee", "lhf", "shf", "gpp", "er", "trans", "ct", "lai", "cveg",
]
const _CALLMIP_COLUMN_VARS = ["swc", "tsoil", "soc"]
const OUTPUT_START_DATE = Date(1997, 1, 1)

# ── Forward model ─────────────────────────────────────────────────────────────

function ClimaCalibrate.forward_model(iteration, member)
    FT = Float64
    climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
    met_nc_path   = joinpath(climaland_dir, "DK_Sor",
                             "DK-Sor_1997-2014_FLUXNET2015_Met.nc")

    # Full simulation: 2-year spinup (1995-1996) + full 1997-2014
    sim_start = DateTime(1996, 1, 1)
    sim_stop  = DateTime(2015, 1, 1)
    @info "CalLMIP simulation: member $member  ($sim_start → $sim_stop)"

    calibrate_params_path =
        ClimaCalibrate.parameter_path(CALLMIP_OUTPUT_DIR, iteration, member)
    toml_dict = LP.create_toml_dict(FT; override_files = [calibrate_params_path])

    land, forcing, ν, θ_r = build_dk_sor_model(
        FT, sim_start, sim_stop, toml_dict, met_nc_path,
    )
    set_ic! = make_dk_sor_ic(forcing.atmos, ν, θ_r)

    # Diagnostics: :short preset + explicit column vars
    output_writer = ClimaDiagnostics.Writers.DictWriter()

    possible = ClimaLand.Diagnostics.get_possible_diagnostics(land)
    req_col  = filter(v -> v in possible, _CALLMIP_COLUMN_VARS)

    diags_short = ClimaLand.default_diagnostics(
        land, sim_start, "";   # outdir unused (DictWriter handles output)
        output_writer,
        output_vars      = :short,
        reduction_period = :daily,
    )
    diags_col = isempty(req_col) ? [] :
        ClimaLand.default_diagnostics(
            land, sim_start, "";   # outdir unused (DictWriter handles output)
            output_writer,
            output_vars      = req_col,
            reduction_period = :daily,
        )
    diags = vcat(diags_short, diags_col)

    simulation = LandSimulation(
        sim_start, sim_stop, Float64(DT), land;
        set_ic!, updateat = Second(DT), user_callbacks = (), diagnostics = diags,
    )
    solve!(simulation)

    member_path = ClimaCalibrate.path_to_ensemble_member(
        CALLMIP_OUTPUT_DIR, iteration, member,
    )
    save_callmip_diagnostics(simulation, member_path)
    return nothing
end

# ── Save ──────────────────────────────────────────────────────────────────────

function save_callmip_diagnostics(simulation, member_path)
    # Surface diagnostics
    surface_data = Dict{String, Any}()
    ref_dates    = nothing
    for var in _CALLMIP_SURFACE_VARS
        try
            dates, vals = extract_daily_diag(simulation, "$(var)_1d_average")
            surface_data[var] = Float64.(vals)
            ref_dates === nothing && (ref_dates = dates)
        catch e
            @warn "Surface diagnostic '$var' unavailable" exception = e
            surface_data[var] = Float64[]
        end
    end

    # Column diagnostics (matrix: n_z × n_days)
    column_data = Dict{String, Any}()
    for var in _CALLMIP_COLUMN_VARS
        try
            dates, mat = extract_daily_column_diag(simulation, "$(var)_1d_average")
            column_data[var] = mat
        catch e
            @warn "Column diagnostic '$var' unavailable" exception = e
            column_data[var] = Matrix{Float64}(undef, 0, 0)
        end
    end

    # Soil z-coordinates
    z_soil = try
        Float64.(vec(parent(
            ClimaCore.Fields.coordinate_field(
                axes(simulation.model.soil.domain.fields.z),
            ).z,
        )))
    catch
        Float64[]
    end

    # Discard spinup — keep 1997-01-01 onwards
    keep = isnothing(ref_dates) ? Bool[] : ref_dates .>= OUTPUT_START_DATE
    kept_dates = isnothing(ref_dates) ? Date[] : ref_dates[keep]
    for (k, v) in surface_data
        if v isa Vector && !isempty(v) && length(v) == length(ref_dates)
            surface_data[k] = v[keep]
        end
    end
    for (k, v) in column_data
        if v isa Matrix && !isempty(v) && size(v, 2) == length(ref_dates)
            column_data[k] = v[:, keep]
        end
    end

    JLD2.jldsave(
        joinpath(member_path, "callmip_diagnostics.jld2");
        dates        = kept_dates,
        surface_data = surface_data,
        column_data  = column_data,
        z_soil       = z_soil,
    )
    @info "Saved callmip_diagnostics.jld2 ($(length(kept_dates)) days)"
end

function extract_daily_column_diag(simulation, diag_name)
    writer = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            writer = d.output_writer
            break
        end
    end
    isnothing(writer) && error("Column diagnostic '$diag_name' not found.")
    (times, raw_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, diag_name)
    dates  = Date.(times isa Vector{DateTime} ? times : date.(times))
    n_days = length(dates)
    first_val = raw_data[1]
    if isa(first_val, Number)
        return dates, reshape(Float64.(raw_data), 1, n_days)
    else
        n_z = length(parent(first_val))
        mat = Matrix{Float64}(undef, n_z, n_days)
        for t in 1:n_days
            mat[:, t] = Float64.(vec(parent(raw_data[t])))
        end
        return dates, mat
    end
end

function ClimaCalibrate.observation_map(iteration)
    return zeros(1, 1)   # placeholder — CalLMIP runs don't use EKP
end

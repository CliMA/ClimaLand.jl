"""
Forward map G(θ) for the CalLMIP Phase 1b multisite calibration.

One member evaluation = one 21-column `ColumnEnsemble` run of the default
`LandModel` (decision D1: out-of-the-box global-map parameters at each tower;
default initial conditions) with the calibrated parameters injected as global
scalars via a ClimaParams TOML override (Phase 1a mechanism). The run covers
(first window year − 1) → end of window; the leading year is water/temperature
spin-up (D5) and is excluded from extraction.

G layout matches observations_phase1b.jl exactly: year-major, then site-major,
then [NEE×12 (gC m⁻² d⁻¹), LHF×12, SHF×12 (W m⁻²)] — monthly all-day means
(Phase 1a GMATCH="allday" convention).

The 21 padded sites (full union axis) are loaded once per process and cached;
per-window forcing is an in-memory crop.
"""

using Dates
using Statistics
import ClimaCore
import ClimaDiagnostics
import ClimaUtilities.TimeManager: date

include(joinpath(@__DIR__, "forcing_phase1b.jl"))

const FT1B = Float64
const DT1B = Float64(450)
const PARAM_NAMES = [
    "pmodel_cstar",
    "pmodel_β_c3",
    "pmodel_β_c4",
    "pmodel_α",
    "moisture_stress_c",
    "soilCO2_reference_rate",
    "soilCO2_activation_energy",
    "michaelis_constant",
    "O2_michaelis_constant",
    "autotrophic_respiration_Rd_ref",
    "relative_contribution_factor",
]
const N_SITES = length(CALIBRATION_SITE_IDS)
const N_PER_SITE_YEAR = 36

const _SITES_FULL = Ref{Any}(nothing)
"""Padded 21-site set on the full union axis, loaded once per process."""
function sites_full()
    isnothing(_SITES_FULL[]) && (_SITES_FULL[] = load_calibration_sites())
    return _SITES_FULL[]
end

function toml_with_overrides(params_vec)
    isempty(params_vec) && return LP.create_toml_dict(FT1B), nothing
    tmpfile = tempname() * ".toml"
    open(tmpfile, "w") do io
        for (nm, v) in zip(PARAM_NAMES, params_vec)
            println(io, "[\"$nm\"]")
            println(io, "value = $(Float64(v))")
            println(io, "type  = \"float\"")
            println(io)
        end
    end
    return LP.create_toml_dict(FT1B; override_files = [tmpfile]), tmpfile
end

"""
    run_member(params_vec, years) -> Vector{Float64}

Evaluate G for one parameter vector over the (consecutive) `years` window.
Returns NaNs of the right length on failure. `params_vec = Float64[]` runs the
default (prior) parameters.
"""
function run_member(params_vec, years)
    n_out = length(years) * N_SITES * N_PER_SITE_YEAR
    toml_dict, tmpfile = toml_with_overrides(params_vec)
    try
        start_date = DateTime(years[1] - 1, 1, 1)   # leading spin-up year
        stop_date = DateTime(years[end] + 1, 1, 1)
        window = (start_date - Day(1), stop_date + Day(1))
        sites = crop_sites(sites_full(), window)

        land_domain = ColumnEnsemble(;
            zlim = FT1B.((-2, 0)),
            nelements = 10,
            longlat = [FT1B.((s.long, s.lat)) for s in sites],
        )
        surface_space = land_domain.space.surface
        forcing = prescribed_forcing_callmip(
            sites,
            surface_space,
            start_date,
            toml_dict,
            FT1B,
        )
        LAI = prescribed_lai_callmip(sites, surface_space, start_date)
        land = Base.invokelatest(
            LandModel{FT1B},
            forcing,
            LAI,
            toml_dict,
            land_domain,
            DT1B;
            prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        )
        @assert !isnothing(land.soilco2) "SoilCO2 component is required (DAMM params + cSoil)"
        writer = ClimaDiagnostics.Writers.DictWriter()
        diags = ClimaLand.default_diagnostics(
            land,
            start_date,
            "";
            output_writer = writer,
            output_vars = ["nee", "lhf", "shf"],
            reduction_period = :daily,
        )
        simulation = Base.invokelatest(
            LandSimulation,
            start_date,
            stop_date,
            DT1B,
            land;
            updateat = Second(DT1B),
            diagnostics = diags,
        )
        solve!(simulation)

        # per-column daily series → monthly all-day means per site
        function monthly_by_site(diag_name, scale)
            times = sort!(collect(keys(writer[diag_name])))
            ddates = Date.(date.(times))
            data = Matrix{Float64}(undef, length(times), N_SITES)
            for (k, t) in enumerate(times)
                data[k, :] = vec(Array(parent(writer[diag_name][t])))[1:N_SITES]
            end
            out = Dict{Tuple{Int, Int}, Vector{Float64}}()
            for y in years, m in 1:12
                sel =
                    (Dates.year.(ddates) .== y) .& (Dates.month.(ddates) .== m)
                out[(y, m)] =
                    any(sel) ? vec(mean(data[sel, :]; dims = 1)) .* scale :
                    fill(NaN, N_SITES)
            end
            return out
        end
        nee = monthly_by_site("nee_1d_average", 12.0 * 86400.0)  # → gC/m²/d
        lhf = monthly_by_site("lhf_1d_average", 1.0)
        shf = monthly_by_site("shf_1d_average", 1.0)

        g = Float64[]
        for y in years, j in 1:N_SITES
            for tab in (nee, lhf, shf), m in 1:12
                push!(g, tab[(y, m)][j])
            end
        end
        @assert length(g) == n_out
        return g
    catch e
        @warn "run_member failed: $(typeof(e)): $e"
        return fill(NaN, n_out)
    finally
        isnothing(tmpfile) || rm(tmpfile; force = true)
    end
end

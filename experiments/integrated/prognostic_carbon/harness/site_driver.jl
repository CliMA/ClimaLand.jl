## Single-column ERA5 driver for the prognostic-carbon test battery.
##
## One column, parametrized entirely through ENV so a single driver can run any
## row of `test_sites.csv`. Ported from `origin/ar/derecho_loop`'s
## `experiments/integrated/era5/single_site.jl`, dropping the `ONLINE_*`
## switches that became permanent behaviour on this branch.
##
## LAI_MODE selects which biomass model sits underneath the carbon model:
## `prescribed` (MODIS, the default for stages 1-3, which isolates carbon-model
## error from LAI-model error) or `prognostic` (ZhouOptimalLAIModel, stage 4 on).
## The carbon model must work under both.
##
## Writes `carbon_metrics.txt` to SITE_OUTDIR: post-spinup annual means of GPP,
## LAI, Ra and Rh. That is the stage-0 baseline every later stage is checked
## against — rule 1 (no change to GPP or LAI) is verified against these numbers.

import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar, evaluate!
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

using Dates
using Statistics

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)

site_lon = parse(FT, get(ENV, "SITE_LON", "5.0"))
site_lat = parse(FT, get(ENV, "SITE_LAT", "25.0"))
# Domains.Column builds the column as a box with a half-width that scales with
# the coordinate magnitude, so a coordinate of exactly 0 gives a degenerate
# zero-width box and trips an assertion in Domains.Plane. Nudge by ~100 m.
site_lon = site_lon == 0 ? FT(1e-3) : site_lon
site_lat = site_lat == 0 ? FT(1e-3) : site_lat
site_name = get(ENV, "SITE_NAME", "site")
start_date = DateTime(get(ENV, "START", "2000-09-01"))
stop_date = DateTime(get(ENV, "STOP", "2002-09-01"))
Δt = parse(FT, get(ENV, "DT", "450.0"))
lai_mode = get(ENV, "LAI_MODE", "prescribed")
lai_mode in ("prescribed", "prognostic") ||
    error("LAI_MODE must be \"prescribed\" or \"prognostic\", got $(lai_mode)")

root_path = get(ENV, "SITE_OUTDIR", "prognostic_carbon_$(site_name)")
mkpath(root_path)

# Optional parameter overrides, applied only when the ENV var is non-empty.
# Restricted to the `optimal_lai_*` parameters that still exist on this branch;
# the f0/vpd_gs/GSL online switches of the derecho_loop driver are gone, their
# behaviour now being unconditional.
override_specs = (
    ("optimal_lai_z", get(ENV, "OPT_Z", "")),
    ("optimal_lai_z_c4", get(ENV, "OPT_Z_C4", "")),
    ("optimal_lai_sigma", get(ENV, "OPT_SIGMA", "")),
    ("optimal_lai_sigma_c4", get(ENV, "OPT_SIGMA_C4", "")),
    ("optimal_lai_alpha", get(ENV, "OPT_ALPHA", "")),
    ("optimal_lai_f0", get(ENV, "F0", "")),
    ("optimal_lai_z_a0", get(ENV, "Z_A0", "")),
    ("optimal_lai_online_c3c4", get(ENV, "ONLINE_C3C4", "")),
)
active_overrides = filter(s -> !isempty(s[2]), collect(override_specs))
override_files = String[]
if !isempty(active_overrides)
    override_path = joinpath(root_path, "override_params.toml")
    open(override_path, "w") do io
        for (nm, ev) in active_overrides
            println(io, "[\"$(nm)\"]")
            println(io, "value = $(parse(FT, ev))")
            println(io, "type = \"float\"")
        end
    end
    push!(override_files, override_path)
end

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
    lai_mode,
) where {FT}
    surface_space = domain.space.surface
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
        use_lowres_forcing = true,
    )
    forcing = (; atmos, radiation)
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    # Passing LAI selects PrescribedBiomassModel; omitting it makes the
    # constructor build the prognostic ZhouOptimalLAIModel instead.
    land = if lai_mode == "prescribed"
        LAI = ClimaLand.Canopy.prescribed_lai_modis(
            surface_space,
            start_date,
            stop_date,
        )
        LandModel{FT}(
            forcing,
            LAI,
            toml_dict,
            domain,
            Δt;
            prognostic_land_components,
        )
    else
        LandModel{FT}(forcing, toml_dict, domain, Δt; prognostic_land_components)
    end
    return land
end

longlat = (site_lon, site_lat)
zlim = FT.((-15, 0))
nelements = 15
dz_tuple = FT.((3, 0.05))
domain = ClimaLand.Domains.Column(; zlim, longlat, nelements, dz_tuple);
surface_space = domain.space.surface
toml_dict = LP.create_toml_dict(FT; override_files)

model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict, lai_mode);

# In-memory diagnostics so the baseline metrics can be pulled straight out of
# the writer without a NetCDF round-trip.
diag_writer = ClimaDiagnostics.Writers.DictWriter()
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date;
    output_writer = diag_writer,
    reduction_period = :daily,
    output_vars = ["gpp", "lai", "ra", "hr", "tair", "swc"],
);
simulation = LandSimulation(start_date, stop_date, Δt, model; diagnostics);

@info "Prognostic-carbon single-column battery site"
@info "Site: $site_name  (lon=$site_lon, lat=$site_lat)  LAI: $lai_mode"
@info "Timestep: $Δt s   Window: $start_date -> $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation);

# Baseline metrics. The running means inside the LAI and climate machinery need
# a year to spin up, so the first SPINUP_YEARS are excluded from the averages.
spinup_date = start_date + Year(parse(Int, get(ENV, "SPINUP_YEARS", "1")))
const M_C = 0.012011  # kg C per mol
# mol CO2 m^-2 s^-1 -> g C m^-2 day^-1
to_gC_per_day(x) = x * M_C * 86400 * 1000

diag_names = [d.output_short_name for d in diagnostics]
series_for(short) = begin
    idx = findfirst(n -> startswith(lowercase(n), short * "_"), diag_names)
    idx === nothing && return nothing
    ClimaLand.Diagnostics.diagnostic_as_vectors(diag_writer, diag_names[idx])
end

open(joinpath(root_path, "carbon_metrics.txt"), "w") do io
    println(io, "site $site_name")
    println(io, "lon $site_lon")
    println(io, "lat $site_lat")
    println(io, "lai_mode $lai_mode")
    println(io, "start $start_date")
    println(io, "stop $stop_date")
    println(io, "dt $Δt")
    for (short, label, conv) in (
        ("gpp", "GPP_gC_m2_day", to_gC_per_day),
        ("ra", "Ra_gC_m2_day", to_gC_per_day),
        ("hr", "Rh_gC_m2_day", to_gC_per_day),
        ("lai", "LAI_m2_m2", identity),
        ("tair", "Tair_K", identity),
        ("swc", "SWC_m3_m3", identity),
    )
        s = series_for(short)
        if s === nothing
            println(io, "$label missing")
            continue
        end
        times, values = s
        dates = date.(times)
        keep = dates .>= spinup_date
        # A run shorter than the spinup still yields a number, so smoke tests
        # produce comparable output instead of an empty metrics file.
        any(keep) || (keep = dates .>= start_date + Day(20))
        println(io, "$label $(conv(mean(values[keep])))")
        println(io, "$(label)_n $(count(keep))")
    end
end

for line in eachline(joinpath(root_path, "carbon_metrics.txt"))
    println("METRIC $site_name $line")
end

# Plots are cheap and land in the site's scratch dir; a failure here must not
# fail a site that simulated fine.
try
    LandSimVis.make_timeseries(root_path, diagnostics, start_date; spinup_date)
catch e
    @warn "timeseries plotting failed (simulation itself succeeded)" exception =
        (e, catch_backtrace())
end

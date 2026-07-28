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
import Thermodynamics
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

ClimaComms.@import_required_backends
const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
# Coordinates and run window are parametrized via ENV so one driver can run any
# row of the test battery (experiments/integrated/generic_site/test_sites.csv).
# Defaults reproduce the validated Sahara point (lon=5, lat=25). lon is degrees
# East; if a negative lon fails the ERA5 regridder, set SITE_LON = lon + 360.
site_lon = parse(FT, get(ENV, "SITE_LON", "5.0"))
site_lat = parse(FT, get(ENV, "SITE_LAT", "25.0"))
# ClimaLand.Domains.Column builds the single-column box with a half-width that
# scales with the coordinate magnitude (xlim=(long,long), ylim=(lat,lat)), so a
# coordinate of exactly 0 gives a degenerate zero-width box and trips an assertion
# in Domains.Plane. Nudge exact-zero coordinates by ~100 m to keep the box valid.
site_lon = site_lon == 0 ? FT(1e-3) : site_lon
site_lat = site_lat == 0 ? FT(1e-3) : site_lat
site_name = get(ENV, "SITE_NAME", "site")
start_date = DateTime(get(ENV, "START", "2000-09-01"))
stop_date = DateTime(get(ENV, "STOP", "2001-09-01"))
Δt = parse(FT, get(ENV, "DT", "450.0"))

root_path = get(ENV, "SITE_OUTDIR", "snowy_land_pmodel_$(site_name)")
mkpath(root_path)

# Toggle whether soil-moisture stress βm enters the potential GPP A0 (see
# optimal_lai_beta_in_A0). BETA_IN_A0=0 (default) makes A0 β-free (a true
# potential GPP); 1 is Zhou's convention. Applied as a parameter override.
beta_in_A0 = parse(FT, get(ENV, "BETA_IN_A0", "0"))
# ONLINE_F0=1 computes f0 online from the aridity index AI=<PET>/<P>; 0 (default)
# keeps the static artifact f0. (F0 is the inert legacy param knob; see the log.)
online_f0 = parse(FT, get(ENV, "ONLINE_F0", "1"))
# F0_SCALE multiplies f0 (calibration knob for the semi-arid LAI over-prediction).
f0_scale = parse(FT, get(ENV, "F0_SCALE", "0.8"))
# ONLINE_VPD_GS=1 computes vpd_gs online as an A0-weighted running-mean VPD.
online_vpd_gs = parse(FT, get(ENV, "ONLINE_VPD_GS", "0"))
# ONLINE_GSL=1 computes GSL online as a trailing-year running count of growing days.
online_gsl = parse(FT, get(ENV, "ONLINE_GSL", "0"))
# ONLINE_C3C4=1 computes fractional_c3 online from the C3/C4 competition.
online_c3c4 = parse(FT, get(ENV, "ONLINE_C3C4", "0"))
f0_env = get(ENV, "F0", "")
f0_setting = isempty(f0_env) ? "default" : f0_env
override_path = joinpath(root_path, "override_params.toml")
open(override_path, "w") do io
    println(io, "[\"optimal_lai_beta_in_A0\"]")
    println(io, "value = $(beta_in_A0)")
    println(io, "type = \"float\"")
    println(io, "[\"optimal_lai_online_f0\"]")
    println(io, "value = $(online_f0)")
    println(io, "type = \"float\"")
    println(io, "[\"optimal_lai_f0_scale\"]")
    println(io, "value = $(f0_scale)")
    println(io, "type = \"float\"")
    println(io, "[\"optimal_lai_online_vpd_gs\"]")
    println(io, "value = $(online_vpd_gs)")
    println(io, "type = \"float\"")
    println(io, "[\"optimal_lai_online_gsl\"]")
    println(io, "value = $(online_gsl)")
    println(io, "type = \"float\"")
    println(io, "[\"optimal_lai_online_c3c4\"]")
    println(io, "value = $(online_c3c4)")
    println(io, "type = \"float\"")
    if !isempty(f0_env)
        println(io, "[\"optimal_lai_f0\"]")
        println(io, "value = $(parse(FT, f0_env))")
        println(io, "type = \"float\"")
    end
    # Optional RMSE-tuning overrides (empty ENV = keep the calibrated default).
    for (nm, ev) in (
        ("optimal_lai_z", get(ENV, "OPT_Z", "")),
        ("optimal_lai_z_c4", get(ENV, "OPT_Z_C4", "")),
        ("optimal_lai_sigma", get(ENV, "OPT_SIGMA", "")),
        ("optimal_lai_sigma_c4", get(ENV, "OPT_SIGMA_C4", "")),
        ("optimal_lai_alpha", get(ENV, "OPT_ALPHA", "")),
        ("optimal_lai_z_a0", get(ENV, "Z_A0", "")),
    )
        if !isempty(ev)
            println(io, "[\"$(nm)\"]")
            println(io, "value = $(parse(FT, ev))")
            println(io, "type = \"float\"")
        end
    end
end

# The goal of this driver is to reproduce MODIS LAI with the prognostic
# optimal-LAI model (Zhou et al. 2025, `ZhouOptimalLAIModel`). We therefore run
# the model with PROGNOSTIC LAI (not prescribed) and compare the modeled LAI to
# the MODIS observations at the same coordinate (RMSE + overlay plot below).
function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
) where {FT}
    surface_space = domain.space.surface
    # Forcing data - high resolution
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

    # Omitting the prescribed-LAI argument makes the LandModel constructor build
    # the prognostic ZhouOptimalLAIModel biomass (the model under test). The
    # default canopy path also threads the soil's own retention ν into the
    # moisture-stress model, so no explicit soil_params wiring is needed here.
    land = LandModel{FT}(
        forcing,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
    )
    return land, atmos
end

longlat = (site_lon, site_lat)
zlim = FT.((-15, 0))
nelements = 15
dz_tuple = FT.((3, 0.05))
domain = ClimaLand.Domains.Column(; zlim, longlat, nelements, dz_tuple);
surface_space = domain.space.surface
toml_dict = LP.create_toml_dict(FT; override_files = [override_path])

model, atmos = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict);
# In-memory diagnostics (DictWriter) so we can pull the modeled LAI timeseries
# directly for the MODIS comparison. `lai` is the modeled leaf area index.
diag_writer = ClimaDiagnostics.Writers.DictWriter()
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date;
    output_writer = diag_writer,
    reduction_period = :daily,
    output_vars = [
        "lai",
        "a0a",
        "a0c3",
        "a0c4",
        "pra",
        "tsoil",
        "swc",
        "hr",
        "gpp",
        "tair",
        "fc3",
    ],
);
simulation = LandSimulation(start_date, stop_date, Δt, model; diagnostics);

@info "Run: prognostic-LAI single-column P-model (reproduce MODIS LAI)"
@info "Site: $site_name  (lon=$site_lon, lat=$site_lat)"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation);

# Compare modeled LAI to MODIS LAI at this coordinate: RMSE + overlay plot.
# Wrapped defensively so a comparison hiccup never fails the simulation itself.
# Spinup: the optimal-LAI RunningMeans (τ=1yr) need a full year to spin up before the
# LAI is scored; A0_annual also self-corrects off its stale artifact value, which can
# take longer (esp. with beta_in_A0=0). Default 1 year; SPINUP_YEARS overrides. Run a
# window at least one year longer than the spinup; short runs fall back below.
spinup_date = start_date + Year(parse(Int, get(ENV, "SPINUP_YEARS", "1")))
try
    diag_names = [d.output_short_name for d in diagnostics]
    lai_idx = findfirst(n -> occursin("lai", lowercase(n)), diag_names)
    lai_idx === nothing && error("no lai diagnostic found in $(diag_names)")

    model_time, model_lai = ClimaLand.Diagnostics.diagnostic_as_vectors(
        diag_writer,
        diag_names[lai_idx],
    )
    model_dates = date.(model_time)

    # MODIS LAI observations at this column, sampled at the model output times.
    lai_obs = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )
    tmp = ClimaCore.Fields.zeros(surface_space)
    modis_lai = map(model_time) do t
        evaluate!(tmp, lai_obs, t)
        first(Array(parent(tmp)))
    end

    # RMSE over the post-spinup window (fall back to a short spinup if the run is
    # under a year, so 1-year test runs still produce a metric).
    keep = model_dates .>= spinup_date
    any(keep) || (keep = model_dates .>= start_date + Day(20))
    rmse = sqrt(mean(abs2, model_lai[keep] .- modis_lai[keep]))
    bias = mean(model_lai[keep] .- modis_lai[keep])
    open(joinpath(root_path, "lai_metrics.txt"), "w") do io
        println(io, "site $site_name")
        println(io, "beta_in_A0 $beta_in_A0")
        println(io, "online_f0 $online_f0")
        println(io, "f0_scale $f0_scale")
        println(io, "online_vpd_gs $online_vpd_gs")
        println(io, "online_gsl $online_gsl")
        println(io, "online_c3c4 $online_c3c4")
        println(io, "opt_z $(get(ENV, "OPT_Z", "default"))")
        println(io, "opt_z_c4 $(get(ENV, "OPT_Z_C4", "default"))")
        println(io, "opt_sigma $(get(ENV, "OPT_SIGMA", "default"))")
        println(io, "opt_sigma_c4 $(get(ENV, "OPT_SIGMA_C4", "default"))")
        println(io, "opt_alpha $(get(ENV, "OPT_ALPHA", "default"))")
        println(io, "f0 $f0_setting")
        println(io, "LAI_RMSE $rmse")
        println(io, "LAI_BIAS $bias")
        println(io, "n $(count(keep))")
        # Energy- vs water-limitation diagnostic for LAI_max (Zhou Eq. 11):
        # fAPAR_energy = 1 - z/(k*A0_annual). Near 1 ⇒ A0 large ⇒ NOT energy-
        # limited ⇒ LAI_max is set by the water term f0*P/A0 (so f0/vpd_gs is the
        # lever); low/negative ⇒ energy/A0-limited.
        a0a_idx = findfirst(n -> occursin("a0a", lowercase(n)), diag_names)
        if a0a_idx !== nothing
            _, a0a_series = ClimaLand.Diagnostics.diagnostic_as_vectors(
                diag_writer,
                diag_names[a0a_idx],
            )
            a0a_mean = mean(a0a_series[keep])
            z = FT(toml_dict["optimal_lai_z"])
            k = FT(toml_dict["optimal_lai_k"])
            fapar_energy = 1 - z / (k * max(a0a_mean, eps(FT)))
            println(io, "A0_annual_mean $a0a_mean")
            println(io, "fAPAR_energy $fapar_energy")
        end
        pra_idx = findfirst(n -> occursin("pra", lowercase(n)), diag_names)
        if pra_idx !== nothing
            _, pra_series = ClimaLand.Diagnostics.diagnostic_as_vectors(
                diag_writer,
                diag_names[pra_idx],
            )
            println(io, "precip_annual_mean $(mean(pra_series[keep]))")
        end
        # Site climate + C3/C4 dominance, for the param-vs-climate investigation:
        # mean annual air temperature (K) and mean modeled C3 fraction.
        for (short, label) in (
            ("tair", "tair_mean"),
            ("fc3", "fractional_c3_mean"),
            ("a0c3", "a0c3_annual_mean"),
            ("a0c4", "a0c4_annual_mean"),
        )
            idx = findfirst(n -> occursin(short, lowercase(n)), diag_names)
            idx === nothing && continue
            _, series = ClimaLand.Diagnostics.diagnostic_as_vectors(
                diag_writer,
                diag_names[idx],
            )
            println(io, "$label $(mean(series[keep]))")
        end
    end
    @info "LAI vs MODIS" site = site_name beta_in_A0 = beta_in_A0 rmse = rmse bias =
        bias
    println(
        "LAI_RMSE $site_name (beta_in_A0=$beta_in_A0) = $rmse  (bias $bias)",
    )

    comparison_data = (; UTC_datetime = model_dates, lai = modis_lai)
    LandSimVis.make_timeseries(
        root_path,
        diagnostics,
        start_date;
        comparison_data,
        spinup_date,
    )
catch e
    @warn "LAI/MODIS comparison failed (simulation itself succeeded)" exception =
        (e, catch_backtrace())
    LandSimVis.make_timeseries(root_path, diagnostics, start_date; spinup_date)
end

# vpd_gs diagnostic (defensive): dump the STATIC artifact vpd_gs (Pa) used in the
# LAI_max water term. Compared across biomes (semi-arid vs tropical), it reveals
# whether the artifact under-represents the semi-arid growing-season aridity —
# in which case online vpd_gs (item [3]) would strengthen water limitation there.
try
    ic = ClimaLand.Canopy.optimal_lai_initial_conditions(surface_space)
    artifact_vpd_gs = first(Array(parent(ic.vpd_gs)))
    open(joinpath(root_path, "vpd_metrics.txt"), "w") do io
        println(io, "site $site_name")
        println(io, "artifact_vpd_gs $artifact_vpd_gs")
    end
    println("VPD $site_name artifact_vpd_gs=$artifact_vpd_gs")
catch e
    @warn "vpd_gs diagnostic failed" exception = (e, catch_backtrace())
end

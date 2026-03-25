"""
Run the calibration forward model at the prior mean parameters for a short
period, to verify the model setup before running the full calibration.

Configuration via environment variables:
    NEON_SITE_ID     — NEON site ID (default: "NEON-srer")
    NEON_SPINUP_DAYS — Number of spinup days (default: 20)

Usage:
    julia --project=.buildkite experiments/calibrate_neon/run_prior_mean.jl
"""

import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date
import ClimaParams as CP
using Insolation

using Dates
using Statistics
using CairoMakie
using CSV
using DataFrames
CairoMakie.activate!()

const NEON_SITE_METADATA_CSV =
    "/kiwi-data/Data/groupMembers/evametz/ERA5/sitedata/NEON_Field_Site_Metadata_20260324.csv"

function _neon_site_key(SITE_ID)
    return uppercase(replace(string(SITE_ID), "NEON_" => "", "NEON-" => ""))
end

function _get_neon_site_metadata(SITE_ID)
    (data, columns) = DelimitedFiles.readdlm(NEON_SITE_METADATA_CSV, ','; header = true)
    header = vec(String.(columns))

    i_site = findfirst(==("site_id"), header)
    i_lat = findfirst(==("latitude"), header)
    i_long = findfirst(==("longitude"), header)
    i_h = findfirst(==("tower_height_m"), header)

    key = _neon_site_key(SITE_ID)
    row_idx = findfirst(i -> uppercase(string(data[i, i_site])) == key, axes(data, 1))

    lat = parse(Float64, string(data[row_idx, i_lat]))
    long = parse(Float64, string(data[row_idx, i_long]))
    atmos_h = parse(Float64, string(data[row_idx, i_h]))
    return (; lat, long, atmos_h)
end

const FT = Float64
const climaland_dir = pkgdir(ClimaLand)
const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
const DT = Float64(450)
outdir = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/$(SITE_ID)_$(Date(start_date))_$(Date(stop_date))_SpinUp$(SPINUP_DAYS)/"
output_base = joinpath(outdir, "output")
outpath = joinpath(output_base, "calibrate_neon_output", "prior_mean_$(SITE_ID).png")

time_offset = 0
(site_start_date, site_stop_date) =
    FluxnetSimulations.get_data_dates(SITE_ID, time_offset)
start_date =
    DateTime(get(ENV, "NEON_START_DATE", string(Date(site_start_date))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(site_stop_date))))

spinup_date = start_date + Day(SPINUP_DAYS)

time_offset = 0
metadata = _get_neon_site_metadata(SITE_ID)
lat = FT(metadata.lat)
long = FT(metadata.long)
atmos_h = FT(metadata.atmos_h)

# ── Prior mean values ─────────────────────────────────────────────────────────
const PRIOR_TOML = joinpath(output_base, "prior_mean_parameters.toml")
open(PRIOR_TOML, "w") do io
    write(io, """
[soilCO2_pre_exponential_factor]
value = 2000.0
type = "float"
used_in = ["Land"]

[michaelis_constant]
value = 0.046
type = "float"
used_in = ["Land"]

[O2_michaelis_constant]
value = 0.066
type = "float"
used_in = ["Land"]
""")
end
println("Wrote prior mean TOML: $PRIOR_TOML")

# ── Setup ─────────────────────────────────────────────────────────────────────
site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)

#(; time_offset, lat, long) =
#    FluxnetSimulations.get_location(FT, Val(site_ID_val))
#(start_date, stop_date) =
#    FluxnetSimulations.get_data_dates(SITE_ID, time_offset)
#spinup_date = start_date + Day(SPINUP_DAYS)

println("Site: $SITE_ID")
println("Simulating $start_date → $stop_date (spinup until $spinup_date)")

# Domain
dz_bottom = FT(2) #FT(1.5),
dz_top = FT(0.038)
dz_tuple = (dz_bottom, dz_top)
nelements = 24
zmin = FT(-6.2)
zmax = FT(0)
# Domain
#(; dz_tuple, nelements, zmin, zmax) =
#    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
#(; atmos_h) =
#    FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
surface_space = land_domain.space.surface

# Determine target layer for soil CO₂ extraction
z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
z_vals = parent(z_field)[:, 1]
target_depth = FT(-0.06)
target_layer = argmin(abs.(z_vals .- target_depth))
println("Target layer: $target_layer (z = $(z_vals[target_layer]) m)")

# Base and calibrated TOML
toml_dict_base = LP.create_toml_dict(FT)
toml_dict = LP.create_toml_dict(FT; override_files = [PRIOR_TOML])

# ERA5 forcing
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    SITE_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict_base,
    FT,
)

# Custom canopy aerodynamic coefficients
toml_dict_base.data["canopy_d_coeff"]["value"] = FT(0.67)
toml_dict_base.data["canopy_z_0b_coeff"]["value"] = FT(0.013)
toml_dict_base.data["canopy_z_0m_coeff"]["value"] = FT(0.13)

# MODIS LAI
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

# ── Build model ──────────────────────────────────────────────────────────────
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
forcing = (; atmos, radiation)
ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)

photosynthesis = PModel{FT}(land_domain, toml_dict_base)
conductance = PModelConductance{FT}(toml_dict_base)
soil_moisture_stress =
    ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict_base)

canopy = ClimaLand.Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    toml_dict_base;
    prognostic_land_components,
    photosynthesis,
    conductance,
    soil_moisture_stress,
)

α_snow = Snow.ZenithAngleAlbedoModel(toml_dict_base)
snow = Snow.SnowModel(
    FT,
    canopy_domain,
    forcing,
    toml_dict_base,
    DT;
    prognostic_land_components,
    α_snow,
)

land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    land_domain,
    DT;
    prognostic_land_components,
    snow,
    canopy,
)

# ── Initial conditions with SOC profile ──────────────────────────────────────
base_set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    SITE_ID,
    start_date,
    time_offset,
    land,
)

function set_ic!(Y, p, t, model)
    base_set_ic!(Y, p, t, model)
    Y.soilco2.CO2 .= FT(0.000412)
    Y.soilco2.O2_f .= FT(0.21)
    SOC_top = FT(18.0)
    SOC_bot = FT(0.5)
    τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
    z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
    @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
end

# ── Diagnostics ──────────────────────────────────────────────────────────────
output_writer = ClimaDiagnostics.Writers.DictWriter()
output_vars = ["swc", "tsoil", "si", "sco2", "soc", "so2", "sco2_ppm"]
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = output_writer,
    output_vars,
    reduction_period = :halfhourly,
)

simulation = LandSimulation(
    start_date,
    stop_date,
    DT,
    land;
    set_ic! = set_ic!,
    updateat = Second(DT),
    diagnostics = diags,
)

println("Running prior-mean model...")
@time solve!(simulation)
println("Done.")

# ── Extract diagnostics ──────────────────────────────────────────────────────
function get_diag_layer(sim, name, layer)
    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        sim.diagnostics[1].output_writer,
        name;
        layer = layer,
    )
    model_dates = times isa Vector{DateTime} ? times : date.(times)
    df = DataFrame(datetime = model_dates, value = Float64.(data))
    df[!, :date] = Date.(df.datetime)
    df = filter(row -> row.date >= Date(spinup_date), df)
    daily = combine(groupby(df, :date), :value => mean => :daily_mean)
    sort!(daily, :date)
    return daily
end

sco2_daily = get_diag_layer(simulation, "sco2_ppm_30m_average", target_layer)
swc_daily = get_diag_layer(simulation, "swc_30m_average", target_layer)
tsoil_daily = get_diag_layer(simulation, "tsoil_30m_average", target_layer)

# ── Load NEON observations ───────────────────────────────────────────────────
csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(SITE_ID)
obs_df = CSV.read(csv_path, DataFrame)

co2_cols_502 = [
    Symbol("soilCO2concentrationMean_001_502"),
    Symbol("soilCO2concentrationMean_002_502"),
    Symbol("soilCO2concentrationMean_003_502"),
    Symbol("soilCO2concentrationMean_004_502"),
    Symbol("soilCO2concentrationMean_005_502"),
]

function rowmean_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        if !ismissing(v) && !isnan(Float64(v))
            push!(vals, Float64(v))
        end
    end
    return isempty(vals) ? NaN : mean(vals)
end

obs_df[!, :sco2_mean_502] =
    [rowmean_skipinvalid(row, co2_cols_502) for row in eachrow(obs_df)]
obs_df[!, :datetime] =
    DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")
obs_df[!, :date] = Date.(obs_df.datetime)

obs_daily = combine(
    groupby(obs_df, :date),
    :sco2_mean_502 =>
        (x -> begin
            valid = filter(!isnan, x)
            length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean,
)
obs_daily = filter(row -> row.date >= Date(spinup_date), obs_daily)
obs_daily = filter(row -> !isnan(row.daily_mean), obs_daily)
sort!(obs_daily, :date)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = Figure(size=(1200, 1000))

ax1 = Axis(fig[1,1]; ylabel="Soil CO₂ (ppm)", title="$SITE_ID — Prior Mean vs NEON Obs (depth 502)")
lines!(ax1, 1:nrow(sco2_daily), sco2_daily.daily_mean; color=:blue, linewidth=1.5, label="Prior mean model")
lines!(ax1, 1:nrow(obs_daily), Float64.(obs_daily.daily_mean); color=:black, linewidth=1.5, label="NEON obs")
axislegend(ax1; position=:rt, framevisible=false)

ax2 = Axis(fig[2,1]; ylabel="SWC")
lines!(ax2, 1:nrow(swc_daily), swc_daily.daily_mean; color=:blue, linewidth=1.5, label="SWC (model)")
axislegend(ax2; position=:rt, framevisible=false)

ax3 = Axis(fig[3,1]; xlabel="Day index (after spinup)", ylabel="Tsoil (K)")
lines!(ax3, 1:nrow(tsoil_daily), tsoil_daily.daily_mean; color=:red, linewidth=1.5, label="Tsoil (model)")
axislegend(ax3; position=:rt, framevisible=false)

mkpath(dirname(outpath))
CairoMakie.save(outpath, fig)
println("Saved: $outpath")

# Print summary stats
println("\nSoil CO₂ stats: min=$(round(minimum(sco2_daily.daily_mean), digits=1)), max=$(round(maximum(sco2_daily.daily_mean), digits=1)), mean=$(round(mean(sco2_daily.daily_mean), digits=1)) ppm")
println("Obs CO₂ stats:  min=$(round(minimum(obs_daily.daily_mean), digits=1)), max=$(round(maximum(obs_daily.daily_mean), digits=1)), mean=$(round(mean(obs_daily.daily_mean), digits=1)) ppm")

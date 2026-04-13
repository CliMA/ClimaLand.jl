# # Single-column test in the Amazon for soil biogeochemistry diagnostics
#
# Runs a single column at (-60, -3) (central Amazon) with ERA5 forcing
# for one month, then plots timeseries of respiration and soil biogeochemistry
# variables at multiple depths.

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates

using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaLand.Simulations: LandSimulation, solve!
using CairoMakie

const FT = Float64

# --- Domain setup ---
# Single column in the central Amazon
lat = FT(-3.0)
lon = FT(-60.0)
zmin = FT(-15.0)   # Same as global run
zmax = FT(0.0)
nelements = 15
dz_tuple = (FT(3.0), FT(0.05))  # Stretched grid: 3m at bottom, 0.05m at top

domain = Column(;
    zlim = (zmin, zmax),
    nelements,
    dz_tuple,
    longlat = (lon, lat),
)
surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
surface_space = domain.space.surface

# --- Time setup ---
start_date = DateTime("2008-03-01")
stop_date = DateTime("2008-04-01")  # 1 month
Δt = 450.0  # 7.5 minutes

# --- Parameters ---
toml_dict = LP.create_toml_dict(FT)
context = ClimaComms.context()

# --- Forcing ---
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    use_lowres_forcing = true,
    context,
)
forcing = (; atmos, radiation)

# --- LAI ---
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    surface_space,
    start_date,
    stop_date,
)

# --- Build model ---
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

# Use the P model for photosynthesis
photosynthesis = PModel{FT}(surface_domain, toml_dict)
conductance = PModelConductance{FT}(toml_dict)
soil_moisture_stress =
    ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(domain, toml_dict)
maxLAI = FT(6.0)  # Typical Amazon max LAI
biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(
    domain,
    LAI,
    maxLAI,
    toml_dict,
)
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    surface_domain,
    (;
        atmos = forcing.atmos,
        radiation = forcing.radiation,
        ground = ClimaLand.PrognosticGroundConditions{FT}(),
    ),
    LAI,
    toml_dict;
    prognostic_land_components,
    photosynthesis,
    conductance,
    soil_moisture_stress,
    biomass,
)

snow = Snow.SnowModel(
    FT,
    surface_domain,
    forcing,
    toml_dict,
    Δt;
    prognostic_land_components,
)

land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components,
    canopy,
    snow,
)

# --- Diagnostics ---
output_vars = [
    "ra",       # autotrophic respiration (surface, mol CO2 m-2 s-1)
    "hr",       # heterotrophic respiration (surface, mol CO2 m-2 s-1)
    "soc",      # soil organic carbon (depth, kg C m-3)
    "sco2_ppm", # soil CO2 in ppm (depth)
    "so2",      # soil O2 fraction (depth, m3/m3)
    "scms",     # soil CO2 microbial source (depth, kg C m-3 s-1)
]

diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = ClimaDiagnostics.Writers.DictWriter(),
    output_vars,
    reduction_period = :daily,
)

simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    land;
    diagnostics = diags,
)

@info "Running Amazon single-column simulation..."
@time solve!(simulation)

# --- Debug: print key state variables ---
Y = simulation._integrator.u
p = simulation._integrator.p
z_coords = parent(domain.fields.z)[:]
@info "Final state (top → bottom):"
@info "  z:          $(round.(z_coords[end:-1:1]; digits=2))"
@info "  CO2 (kgC/m³ soil): $(round.(parent(Y.soilco2.CO2)[end:-1:1]; sigdigits=4))"
@info "  O2_f:       $(round.(parent(Y.soilco2.O2_f)[end:-1:1]; sigdigits=4))"
@info "  SOC:        $(round.(parent(Y.soilco2.SOC)[end:-1:1]; sigdigits=4))"
@info "  θ_a:        $(round.(parent(p.soilco2.θ_a)[end:-1:1]; sigdigits=4))"
@info "  θ_eff:      $(round.(parent(p.soilco2.θ_eff)[end:-1:1]; sigdigits=4))"
@info "  CO2_air_eq: $(round.(parent(p.soilco2.CO2_air_eq)[end:-1:1]; sigdigits=4))"
@info "  D (diffusivity): $(round.(parent(p.soilco2.D)[end:-1:1]; sigdigits=4))"
@info "  Sm (source): $(round.(parent(p.soilco2.Sm)[end:-1:1]; sigdigits=4))"
@info "  top_bc:     $(round.(parent(p.soilco2.top_bc); sigdigits=4))"
@info "  θ_l:        $(round.(parent(Y.soil.ϑ_l)[end:-1:1]; sigdigits=4))"
# Expected atmos CO2_air_eq: c_co2 * P * M_C / (R * T) ≈ 4.2e-4 * 1e5 * 0.012 / (8.314 * 300) ≈ 2e-4 kgC/m³
@info "  c_co2 (atm): $(parent(p.drivers.c_co2))"
@info "  P_sfc:      $(parent(p.drivers.P))"

# --- Plotting ---
writer = first(diags).output_writer

# Get depth coordinates from the model
z = parent(domain.fields.z)[:]   # z coordinates of cell centers
nlayers = length(z)

# Select a few layers to plot for depth-resolved variables
# Layer numbering: 1 = bottom, nlayers = top
layer_top = nlayers       # near surface
layer_mid1 = nlayers - 2  # ~shallow
layer_mid2 = nlayers - 5  # mid-depth
layer_deep = max(1, nlayers - 10)  # deeper
selected_layers = [layer_top, layer_mid1, layer_mid2, layer_deep]
layer_labels = ["z = $(round(z[l]; digits=2)) m" for l in selected_layers]

# Helper to extract timeseries
# DictWriter keys use the format "{short_name}_1d_average" for daily diagnostics
function get_ts(writer, varname; layer = nothing)
    key = "$(varname)_1d_average"
    return ClimaLand.Diagnostics.diagnostic_as_vectors(
        writer,
        key;
        layer,
    )
end

# Convert ITime vector to days since start
function time_to_days(times)
    t0 = times[1]
    return [Float64((t - t0).counter) / 86400.0 for t in times]
end

# Create figure
fig = Figure(size = (1400, 1800))

# Unit conversion: mol CO2 m-2 s-1 → gC m-2 d-1
# 1 mol CO2 = 12.011 g C, 1 day = 86400 s
const mol_co2_to_gC_per_day = 12.011 * 86400.0

# --- Panel 1: Autotrophic Respiration (surface) ---
ax1 = Axis(
    fig[1, 1],
    xlabel = "Days",
    ylabel = "Ra (gC m⁻² d⁻¹)",
    title = "Autotrophic Respiration",
)
times_ra, vals_ra = get_ts(writer, "ra")
lines!(ax1, time_to_days(times_ra), vals_ra .* mol_co2_to_gC_per_day)

# --- Panel 2: Heterotrophic Respiration (surface) ---
ax2 = Axis(
    fig[1, 2],
    xlabel = "Days",
    ylabel = "HR (gC m⁻² d⁻¹)",
    title = "Heterotrophic Respiration",
)
times_hr, vals_hr = get_ts(writer, "hr")
lines!(ax2, time_to_days(times_hr), vals_hr .* mol_co2_to_gC_per_day)

# --- Panel 3: SOC at different depths ---
ax3 = Axis(
    fig[2, 1],
    xlabel = "Days",
    ylabel = "SOC (kg C m⁻³)",
    title = "Soil Organic Carbon",
)
for (i, l) in enumerate(selected_layers)
    times_soc, vals_soc = get_ts(writer, "soc"; layer = l)
    lines!(ax3, time_to_days(times_soc), vals_soc; label = layer_labels[i])
end
axislegend(ax3; position = :rt)

# --- Panel 4: Soil CO2 in ppm at different depths ---
ax4 = Axis(
    fig[2, 2],
    xlabel = "Days",
    ylabel = "Soil CO₂ (ppm)",
    title = "Soil CO₂",
)
for (i, l) in enumerate(selected_layers)
    times_sco2, vals_sco2 = get_ts(writer, "sco2_ppm"; layer = l)
    lines!(ax4, time_to_days(times_sco2), vals_sco2; label = layer_labels[i])
end
axislegend(ax4; position = :rt)

# --- Panel 5: Soil O2 at different depths ---
ax5 = Axis(
    fig[3, 1],
    xlabel = "Days",
    ylabel = "O₂ fraction (m³/m³)",
    title = "Soil O₂ Volumetric Fraction",
)
for (i, l) in enumerate(selected_layers)
    times_so2, vals_so2 = get_ts(writer, "so2"; layer = l)
    lines!(ax5, time_to_days(times_so2), vals_so2; label = layer_labels[i])
end
axislegend(ax5; position = :rb)

# --- Panel 6: Microbial source at different depths ---
ax6 = Axis(
    fig[3, 2],
    xlabel = "Days",
    ylabel = "Sm (mgC m⁻³ s⁻¹)",
    title = "Microbial CO₂ Source",
)
for (i, l) in enumerate(selected_layers)
    times_sm, vals_sm = get_ts(writer, "scms"; layer = l)
    lines!(ax6, time_to_days(times_sm), vals_sm .* 1e6; label = layer_labels[i])
end
axislegend(ax6; position = :rt)

savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments/integrated/amazon_column_soilco2_out",
)
mkpath(savedir)
save(joinpath(savedir, "amazon_soilco2_timeseries.png"), fig)
@info "Saved plot to $savedir"
display(fig)

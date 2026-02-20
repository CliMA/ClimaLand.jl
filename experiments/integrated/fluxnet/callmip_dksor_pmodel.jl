# ── DK-Sor CalMIP Phase 1a: PModel Configuration ─────────────────────────────
#
# This script runs the DK-Sor (Denmark Sorø) CalMIP Phase 1a site using the
# high-level LandModel constructor with PModel photosynthesis overrides.
#
# Overridden sub-models:
#   - Photosynthesis: PModel (optimality-based)
#   - Conductance: PModelConductance (coupled with PModel)
#   - Soil moisture stress: PiecewiseMoistureStressModel (θ-based, not ψ-based)
#
# Default sub-models (from LandModel constructor):
#   - Radiative transfer: TwoStreamModel
#   - Plant hydraulics: PlantHydraulicsModel
#   - Energy: BigLeafEnergyModel
#   - Snow: SnowModel (with default albedo)
#   - Soil: EnergyHydrology (spatial defaults)
#   - SoilCO2: SoilCO2Model
#
# For more information about the CalMIP project:
# https://github.com/callmip-org/Phase1
# https://callmip-org.github.io

# ── 1. Imports ────────────────────────────────────────────────────────────────
import ClimaLand
import SciMLBase
using ClimaCore
import ClimaComms
import ClimaParams as CP
using Dates
using Insolation

using ClimaLand
using ClimaLand: LandModel
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaDiagnostics
using ClimaUtilities

# ── 2. Configuration ──────────────────────────────────────────────────────────
const FT = Float64
site_ID = "DK-Sor"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
climaland_dir = pkgdir(ClimaLand)
toml_dict = LP.create_toml_dict(FT)

# Simulation parameters
dt = Float64(450) # 7.5 minutes
# CalMIP Phase 1a simulation period (must be within MODIS LAI coverage: 2000-2020)
start_date = DateTime(2008, 1, 1)
stop_date = DateTime(2009, 1, 1) # Full year 2008

# ── 3. Site metadata ──────────────────────────────────────────────────────────
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

# ── 4. Domain ─────────────────────────────────────────────────────────────────
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# ── 5. Atmospheric & radiation forcing ────────────────────────────────────────
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT;
    split_precip = :callmip,
)

# ── 6. LAI from MODIS ─────────────────────────────────────────────────────────
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

# ── 7. Custom canopy components (PModel overrides) ────────────────────────────
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
forcing = (; atmos, radiation)

# PModel photosynthesis (optimality-based, requires PModelConductance)
photosynthesis = PModel{FT}(canopy_domain, toml_dict)
conductance = PModelConductance{FT}(toml_dict)

# Soil moisture stress based on soil moisture (θ) rather than leaf water potential (ψ)
soil_moisture_stress =
    ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)

# Construct canopy with PModel overrides; other components use defaults
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    toml_dict;
    prognostic_land_components,
    photosynthesis,
    conductance,
    soil_moisture_stress,
)

# ── 8. Build LandModel with custom canopy ─────────────────────────────────────
land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    land_domain,
    dt;
    prognostic_land_components,
    canopy,
)

# ── 9. Initial conditions ────────────────────────────────────────────────────
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)

# ── 10. Diagnostics ──────────────────────────────────────────────────────────
outdir = joinpath(
    pkgdir(ClimaLand),
    "experiments/integrated/fluxnet/$(site_ID)/callmip_pmodel/out",
)
mkpath(outdir)
println("Output directory: ", outdir)

output_writer = ClimaDiagnostics.Writers.NetCDFWriter(
    land_domain.space.subsurface,
    outdir;
    start_date,
)

short_names_1D = [
    "sif",
    "ra",
    "gs",
    "gpp",
    "ct",
    "swu",
    "lwu",
    "er",
    "et",
    "msf",
    "shf",
    "lhf",
    "rn",
]
short_names_2D = ["swc", "tsoil", "si", "sco2", "soc", "so2"]
output_vars = [short_names_1D..., short_names_2D..., "swe"]

diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = output_writer,
    output_vars,
    reduction_period = :halfhourly,
)

# ── 11. Run simulation ───────────────────────────────────────────────────────
simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    land;
    set_ic!,
    updateat = Second(dt),
    diagnostics = diags,
)
solve!(simulation)

# ── 12. Post-processing ──────────────────────────────────────────────────────
println("\nSimulation completed successfully!")
println("Start date: ", start_date)
println("Stop date: ", stop_date)

ClimaLand.Diagnostics.close_output_writers(diags)
println("Output files written to: ", outdir)

for file in readdir(outdir; join = false)
    filepath = joinpath(outdir, file)
    if isfile(filepath)
        filesize_mb = round(stat(filepath).size / 1024^2, digits = 2)
        println("  $file ($(filesize_mb) MB)")
    end
end

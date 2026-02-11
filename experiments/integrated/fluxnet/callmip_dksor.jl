"""
This experiment runs the DK-Sor (Denmark Sorø) CalMIP Phase 1a test site.

The DK-Sor site is a deciduous broadleaf forest in Denmark, dominated by beech trees.
This is the test calibration site for CalMIP Phase 1a.

For more information about the CalMIP project and protocols:
https://github.com/callmip-org/Phase1
https://callmip-org.github.io

Citation: CalMIP Phase 1 Protocol v1.1
"""

import ClimaLand
import SciMLBase
using ClimaCore
import ClimaComms
import ClimaParams as CP
using Dates
using Insolation

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaDiagnostics
using ClimaUtilities

using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
# using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
# import ClimaLand.LandSimVis as LandSimVis

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)

# Site identification
site_ID = "DK-Sor"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

# Get site parameters
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_vg_n,
    soil_vg_α,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_α_PAR,
    soil_α_NIR,
    Ω,
    χl,
    α_PAR_leaf,
    λ_γ_PAR,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    ac_canopy,
    g1,
    Drel,
    g0,
    Vcmax25,
    SAI,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    Weibull_param,
    a,
    conductivity_model,
    retention_model,
    plant_ν,
    plant_S_s,
    rooting_depth,
    n_stem,
    n_leaf,
    h_leaf,
    h_stem,
    h_canopy,
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

compartment_midpoints =
    n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
compartment_surfaces = n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

# Construct the ClimaLand domain to run the simulation on
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
prognostic_land_components = (:canopy, :soil, :soilco2)

# Set up the timestepping information for the simulation
dt = Float64(450) # 7.5 minutes
# Define simulation period - short 1-day test
start_date = DateTime(2008, 1, 1)
stop_date = DateTime(2009, 1, 1)  # Full year 2008

# This reads in the data from the CalMIP site and creates
# the atmospheric and radiative driver structs for the model.
# CalMIP-format data for this site are handled via the shared
# Fluxnet/CalMIP data-reading utilities used by `prescribed_forcing_fluxnet`.
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

# Set up soil model
soil_domain = land_domain
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)

runoff = ClimaLand.Soil.Runoff.SurfaceRunoff()
retention_parameters = (;
    ν = soil_ν,
    θ_r,
    K_sat = soil_K_sat,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
)
composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)
soil_forcing = (; atmos, radiation)
soil = Soil.EnergyHydrology{FT}(
    soil_domain,
    soil_forcing,
    toml_dict;
    prognostic_land_components,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
    albedo = soil_albedo,
    runoff,
    retention_parameters,
    composition_parameters,
    S_s = soil_S_s,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
)

# Soil microbes model
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(soil_domain, drivers, toml_dict)

# Set up canopy radiative transfer
radiation_parameters = (;
    Ω,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
)
radiative_transfer = Canopy.TwoStreamModel{FT}(
    canopy_domain,
    toml_dict;
    radiation_parameters,
    ϵ_canopy,
)

# Set up stomatal conductance
conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain, toml_dict; g1)

# Set up photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis =
    FarquharModel{FT}(canopy_domain, toml_dict; photosynthesis_parameters)

# Set up plant hydraulics
surface_space = land_domain.space.surface

# Use prescribed LAI from MODIS data or constant LAI for testing
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)
RAI = maxLAI * f_root_to_shoot

hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain,
    toml_dict;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    retention_model,
)

height = h_stem + h_leaf
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)

# Set up canopy energy model
energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)

# Construct the canopy model
canopy = Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    toml_dict;
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    hydraulics,
    energy,
    biomass,
)

# Integrated land model
land = SoilCanopyModel{FT}(soilco2, soil, canopy)

# Set initial conditions
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)

# Set up diagnostics
# Define output directory for diagnostics
outdir = joinpath(
    pkgdir(ClimaLand),
    "experiments/integrated/fluxnet/$(site_ID)/callmip/out",
)
mkpath(outdir)
println("Creating NetCDF output writer for directory: ", outdir)

# Use NetCDFWriter to save diagnostics to NetCDF files
output_writer = ClimaDiagnostics.Writers.NetCDFWriter(
    soil_domain.space.subsurface,
    outdir;
    start_date,
)

# CalMIP requires specific output variables - ensure these match protocol requirements
short_names_1D = [
    "sif",       # Solar-induced fluorescence
    "ra",        # Autotrophic respiration
    "gs",        # Stomatal conductance
    "gpp",       # Gross primary productivity
    "ct",        # Canopy temperature
    "swu",       # Upwelling shortwave radiation
    "lwu",       # Upwelling longwave radiation
    "er",        # Ecosystem respiration
    "et",        # Evapotranspiration
    "msf",       # Moisture surface flux
    "shf",       # Sensible heat flux
    "lhf",       # Latent heat flux
    "rn",        # Net radiation
]
short_names_2D = [
    "swc",       # Soil water content
    "tsoil",     # Soil temperature
    "si",        # Snow/ice
]
output_vars = [short_names_1D..., short_names_2D...]

diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = output_writer,
    output_vars,
    reduction_period = :halfhourly,
)

# Run the simulation
simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    land;
    set_ic! = set_ic!,
    updateat = Second(dt),
    diagnostics = diags,
)
solve!(simulation)

# Post-processing
println("\nSimulation completed successfully!")
println("Start date: ", start_date)
println("Stop date: ", stop_date)
println("Number of diagnostics: ", length(diags))

# Close the output writers to flush diagnostics to files
ClimaLand.Diagnostics.close_output_writers(diags)
println("Output files written to: ", outdir)

# List output files
println("Output files:")
for file in readdir(outdir; join=false)
    filepath = joinpath(outdir, file)
    if isfile(filepath)
        filesize_mb = round(stat(filepath).size / 1024^2, digits=2)
        println("  $file ($(filesize_mb) MB)")
    end
end

# TODO: Add visualization when CairoMakie is available
# TODO: Add comparison with observation data from Flux NetCDF file

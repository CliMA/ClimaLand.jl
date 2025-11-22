"""
This experiment tests running the Ozark site (US-MOz) using plant parameters
defined by plant functional types instead of fully site-specific parameters.
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
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis
const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
site_ID = "US-MOz"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)


# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
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
    z_0m,
    z_0b,
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
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
# Define the PFT land cover percentages for the Ozark site. Currently we only
# use the dominant PFT, which for Ozark is deciduous broadleaf temperate trees.
pft_pcts = [
    0.0, # NET_Temp
    0.0, # NET_Bor
    0.0, # NDT_Bor
    0.0, # BET_Trop
    0.0, # BET_Temp
    0.0, # BDT_Trop
    1.0, # BDT_Temp
    0.0, # BDT_Bor
    0.0, # BES_Temp
    0.0, # BDS_Temp
    0.0, # BDT_Bor
    0.0, # C3G_A
    0.0, # C3G_NA
    0.0, # C4G
]

# Load the PFT parameters into the namespace
(
    Ω,
    α_PAR_leaf,
    α_NIR_leaf,
    τ_PAR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    χl,
    ac_canopy,
    g1,
    Vcmax25,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    plant_ν,
    rooting_depth,
) = FT.(params_from_pfts(pft_pcts))

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model

(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
)

# Now we set up the model. For the soil model, we pick
# a model type and model args:
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
soil_organic_carbon =
    ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers = Soil.Biogeochemistry.SoilDrivers(
    co2_prognostic_soil,
    soil_organic_carbon,
    atmos,
)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(soil_domain, drivers, toml_dict)

# Canopy model - Set up individual Component arguments with non-default parameters
# Set up radiative transfer
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

# Set up conductance
conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain, toml_dict; g1)

# Set up photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis =
    FarquharModel{FT}(canopy_domain, toml_dict; photosynthesis_parameters)

# Set up optimal LAI model
lai_model =
    Canopy.OptimalLAIModel{FT}(Canopy.OptimalLAIParameters{FT}(toml_dict))

# Set up plant hydraulics
# Read in LAI from MODIS data
surface_space = land_domain.space.surface;

LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long);
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

# Set up energy model
energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
# Construct the canopy model using defaults for autotrophic respiration and SIF models
canopy = Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    toml_dict;
    z_0m,
    z_0b,
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    lai_model,
    hydraulics,
    energy,
    biomass,
)

# Integrated plant hydraulics and soil model
land = SoilCanopyModel{FT}(soilco2, soil, canopy)
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)
# Callbacks
output_writer = ClimaDiagnostics.Writers.DictWriter()

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
short_names_2D = ["swc", "tsoil", "si"]
output_vars = [short_names_1D..., short_names_2D...]

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
    dt,
    land;
    set_ic! = set_ic!,
    updateat = Second(dt), # How often we want to update the drivers
    diagnostics = diags,
)
solve!(simulation)

comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
savedir =
    joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/US-MOz/pft/out")
mkpath(savedir)
LandSimVis.make_diurnal_timeseries(
    land_domain,
    diags,
    start_date;
    savedir,
    short_names = ["gpp", "shf", "lhf", "swu", "lwu"],
    spinup_date = start_date + Day(20),
    comparison_data,
)
LandSimVis.make_timeseries(
    land_domain,
    diags,
    start_date;
    savedir,
    short_names = ["swc", "tsoil"],
    spinup_date = start_date + Day(20),
    comparison_data,
)

"""
This experiment tests running the Ozark site (US-MOz) using plant parameters
defined by plant functional types instead of fully site-specific parameters.
"""

import ClimaLand

import SciMLBase
import ClimaTimeSteppers as CTS
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
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaDiagnostics
using ClimaUtilities

using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Poppler_jll, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis
const FT = Float64
earth_param_set = LP.LandParameters(FT)
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
    z0_m,
    z0_b,
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
t0 = Float64(0)
dt = Float64(450) # 7.5 minutes
N_spinup_days = 15
N_days = N_spinup_days + 340
tf = Float64(t0 + N_days * 3600 * 24)
t_spinup = Float64(t0 + N_spinup_days * 3600 * 24)

timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

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

# For now, replace ac_canopy with larger value in order to increase
# stability of timestepping
ac_canopy = ac_canopy * 3
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
start_date = DateTime(2010) + Hour(time_offset)
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
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
    earth_param_set;
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
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(soil_domain, drivers)

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
radiative_transfer =
    Canopy.TwoStreamModel{FT}(canopy_domain; radiation_parameters, ϵ_canopy)

# Set up conductance
conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain; g1)

# Set up photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis = FarquharModel{FT}(canopy_domain; photosynthesis_parameters)

# Set up plant hydraulics
# Read in LAI from MODIS data
surface_space = land_domain.space.surface;
modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(
    start_date = start_date + Second(t0),
    end_date = start_date + Second(t0) + Second(tf);
    context = ClimaComms.context(surface_space),
);
LAI = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date,
);
# Get the maximum LAI at this site over the first year of the simulation
maxLAI =
    FluxnetSimulations.get_maxLAI_at_site(modis_lai_ncdata_path[1], lat, long);
RAI = maxLAI * f_root_to_shoot
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain,
    LAI;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    SAI,
    RAI,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    retention_model,
    rooting_depth,
)

# Set up energy model
energy = Canopy.BigLeafEnergyModel{FT}(; ac_canopy)

ground = ClimaLand.PrognosticSoilConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
# Construct the canopy model using defaults for autotrophic respiration and SIF models
canopy = Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    earth_param_set;
    z_0m = z0_m,
    z_0b = z0_b,
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    hydraulics,
    energy,
)

# Integrated plant hydraulics and soil model
land = SoilCanopyModel{FT}(soilco2, soil, canopy)

Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land);
jacobian! = make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

FluxnetSimulations.set_fluxnet_ic!(Y, site_ID, start_date, time_offset, land)
set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0);

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
    average_period = :hourly,
)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0, dt = dt);

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler);

## How often we want to update the drivers. Note that this uses the defined `t0` and `tf`
## defined in the simulatons file
data_dt = Float64(FluxnetSimulations.get_data_dt(site_ID))
updateat = Array(t0:data_dt:tf)
model_drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, diag_cb)


prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
);

sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb);

ClimaLand.Diagnostics.close_output_writers(diags)
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
    spinup_date = start_date + Day(N_spinup_days),
    comparison_data,
)
LandSimVis.make_timeseries(
    land_domain,
    diags,
    start_date;
    savedir,
    short_names = ["swc", "tsoil"],
    spinup_date = start_date + Day(N_spinup_days),
    comparison_data,
)

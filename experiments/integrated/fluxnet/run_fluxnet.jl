import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using Insolation

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaDiagnostics
using ClimaUtilities

using DelimitedFiles
FluxnetSimulationsExt =
    Base.get_extension(ClimaLand, :FluxnetSimulationsExt).FluxnetSimulationsExt;
using CairoMakie, ClimaAnalysis, GeoMakie, Poppler_jll, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64
earth_param_set = LP.LandParameters(FT)
climaland_dir = pkgdir(ClimaLand)

# Read in the site to be run from the command line
if length(ARGS) < 1
    error("Must provide site ID on command line")
end

site_ID = ARGS[1]

# Read all site-specific domain parameters from the simulation file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_simulation.jl",
    ),
)

include(
    joinpath(climaland_dir, "experiments/integrated/fluxnet/fluxnet_domain.jl"),
)

# Read all site-specific parameters from the parameter file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_parameters.jl",
    ),
)

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/fluxnet_simulation.jl",
    ),
)
(start_date, end_date) =
    FluxnetSimulationsExt.get_data_dates(site_ID, time_offset)
(; atmos, radiation) = FluxnetSimulationsExt.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT,
)

# Read in LAI from MODIS data
surface_space = land_domain.space.surface
modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
    context = ClimaComms.context(surface_space),
    start_date = start_date + Second(t0),
    end_date = start_date + Second(t0) + Second(tf),
)
LAI = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date,
)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulationsExt.get_maxLAI_at_site(
    modis_lai_ncdata_path[1],
    lat,
    long,
);
RAI = maxLAI * f_root_to_shoot
capacity = plant_ν * maxLAI * h_leaf * FT(1000)

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)
soil_ps = Soil.EnergyHydrologyParameters(
    FT;
    ν = soil_ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    K_sat = soil_K_sat,
    S_s = soil_S_s,
    θ_r,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
    albedo = soil_albedo,
);

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}
soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
soilco2_args = (; domain = soil_domain, parameters = soilco2_ps)

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
    radiative_transfer = Canopy.TwoStreamModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
    energy = Canopy.BigLeafEnergyModel{FT},
)
# Individual Component arguments
# Set up autotrophic respiration
autotrophic_respiration_args =
    (; parameters = AutotrophicRespirationParameters(FT))
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = TwoStreamParameters(
        FT;
        Ω,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        G_Function,
    )
)
# Set up conductance
conductance_args = (; parameters = MedlynConductanceParameters(FT; g1))
# Set up photosynthesis
is_c3 = FT(1) # set the photosynthesis mechanism to C3
photosynthesis_args =
    (; parameters = FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAI, SAI, RAI)

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    rooting_depth = rooting_depth,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
)
plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

# Canopy component args
canopy_component_args = (;
    autotrophic_respiration = autotrophic_respiration_args,
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
    energy = energy_args,
)

# Other info needed
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = canopy_domain)

# Snow model
snow_parameters = SnowParameters{FT}(dt; earth_param_set = earth_param_set);
snow_args = (; parameters = snow_parameters, domain = canopy_domain);
snow_model_type = Snow.SnowModel
# Integrated plant hydraulics and soil model
land_input = (
    atmos = atmos,
    radiation = radiation,
    soil_organic_carbon = Csom,
    runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
)
land = LandModel{FT}(;
    soilco2_type = soilco2_type,
    soilco2_args = soilco2_args,
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
    snow_args = snow_args,
    snow_model_type = snow_model_type,
)

Y, p, cds = initialize(land)

FluxnetSimulationsExt.set_fluxnet_ic!(Y, site_ID, start_date, time_offset, land)
set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0);

exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)
jacobian! = make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);


# Callbacks
output_vars = [
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
    "swe",
    "swc",
    "tsoil",
    "si",
]
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = ClimaDiagnostics.Writers.DictWriter(),
    output_vars,
    average_period = :hourly,
);

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0, dt = dt);

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler);

## How often we want to update the drivers. Note that this uses the defined `t0`, and `tf`
## defined in the simulatons file
data_dt = Float64(FluxnetSimulationsExt.get_data_dt(site_ID));
updateat = Array(t0:data_dt:tf);
model_drivers = ClimaLand.get_drivers(land);
updatefunc = ClimaLand.make_update_drivers(model_drivers);
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc);
cb = SciMLBase.CallbackSet(driver_cb, diag_cb);


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

@time sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb);

ClimaLand.Diagnostics.close_output_writers(diags)
comparison_data =
    FluxnetSimulationsExt.get_comparison_data(site_ID, time_offset)
savedir =
    joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/$(site_ID)/out")
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
    short_names = ["swc", "tsoil", "swe"],
    spinup_date = start_date + Day(N_spinup_days),
    comparison_data,
)

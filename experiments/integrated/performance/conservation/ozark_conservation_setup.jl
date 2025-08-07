# Convert site_ID string to a Val so we can dispatch on it
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
# Use finer dz_top for conservation test
dz_top = FT(0.025)
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val); dz_top)
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
    G_Function,
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

# TIME STEPPING:
t0 = Float64(120 * 3600 * 24)# start mid year
dt = Float64(150)
# Use smaller `tf` for Float32 simulation
tf = (FT == Float64) ? t0 + 3600 * 24 * 10 : t0 + 2 * dt

# Select conv. condition based on float type due to different precision
err = (FT == Float64) ? 1e-8 : 1e-4
norm_condition = CTS.MaximumError(err)
conv_checker = CTS.ConvergenceChecker(; norm_condition)
max_iterations = 20
timestepper = CTS.ARS111()
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = max_iterations,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

# Set up the domain for the simulation
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# Get the atmospheric and radiation forcing data
(start_date, end_date) = FluxnetSimulations.get_data_dates(site_ID, time_offset)
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
maxLAI =
    FluxnetSimulations.get_maxLAI_at_site(modis_lai_ncdata_path[1], lat, long);
RAI = FT(maxLAI) * f_root_to_shoot # convert to float type of simulation
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
)

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
energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

# Set up autotrophic respiration
autotrophic_respiration_args =
    (; parameters = AutotrophicRespirationParameters(FT))
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = TwoStreamParameters(
        FT;
        Ω,
        G_Function,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
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

# Integrated plant hydraulics and soil model
land_input = (atmos = atmos, radiation = radiation, soil_organic_carbon = Csom)
land = SoilCanopyModel{FT}(;
    soilco2_type = soilco2_type,
    soilco2_args = soilco2_args,
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land);
jacobian! = make_jacobian(land);
set_initial_cache! = make_set_initial_cache(land)
Y, p, cds = initialize(land)
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

FluxnetSimulations.set_fluxnet_ic!(Y, site_ID, start_date, time_offset, land)
set_initial_cache!(p, Y, t0)

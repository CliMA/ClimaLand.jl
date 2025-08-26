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

prognostic_land_components = (:canopy, :soil, :soilco2)

# Get the atmospheric and radiation forcing data
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
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

# Soil model
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)
retention_parameters = (;
    ν = soil_ν,
    θ_r,
    K_sat = soil_K_sat,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
)
composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)
soil_forcing = (; atmos, radiation)
soil = Soil.EnergyHydrology{FT}(
    land_domain,
    soil_forcing,
    toml_dict;
    prognostic_land_components,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
    albedo = soil_albedo,
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
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, drivers)

# Canopy model
# Set up radiative transfer
radiation_parameters =
    (; Ω, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf)
radiative_transfer =
    Canopy.TwoStreamModel{FT}(canopy_domain; radiation_parameters, ϵ_canopy)

# Set up conductance
conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain, toml_dict; g1)

# Set up photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis = FarquharModel{FT}(canopy_domain; photosynthesis_parameters)

# Set up plant hydraulics
# Read in LAI from MODIS data
surface_space = land_domain.space.surface
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    surface_space,
    start_date + Second(t0),
    stop_date + Second(tf),
)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long);
RAI = FT(maxLAI) * f_root_to_shoot # convert to float type of simulation

hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain,
    LAI,
    toml_dict;
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
energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
# Construct the canopy model using defaults for autotrophic respiration and SIF models
canopy = Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    toml_dict;
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
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land);
jacobian! = make_jacobian(land);
set_initial_cache! = make_set_initial_cache(land)
Y, p, cds = initialize(land)
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)
set_ic!(Y, p, t0, land)
set_initial_cache!(p, Y, t0)

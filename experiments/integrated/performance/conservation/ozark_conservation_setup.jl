# Column dimensions - separation of layers at the top and bottom of the column:
dz_bottom = FT(1.5)
dz_top = FT(0.025)

# Stem and leaf compartments and their heights:
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9) # m
h_leaf = FT(9.5) # m

# TIME STEPPING:
t0 = Float64(120 * 3600 * 24)# start mid year
dt = Float64(150)
# Use smaller `tf` for Float32 simulation
tf = (FT == Float64) ? t0 + 3600 * 24 * 10 : t0 + 2 * dt

dz_tuple = (dz_bottom, dz_top)
nelements = 20
zmin = FT(-10)
zmax = FT(0)

timestepper = CTS.ARS111()
# Select conv. condition based on float type due to different precision
err = (FT == Float64) ? 1e-8 : 1e-4
norm_condition = CTS.MaximumError(err)
conv_checker = CTS.ConvergenceChecker(; norm_condition)
max_iterations = 20
# Set up timestepper
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = max_iterations,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

# Setup the domain for the simulation
include(
    joinpath(climaland_dir, "experiments/integrated/fluxnet/fluxnet_domain.jl"),
)
# Read in the parameters for the Ozark site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$(site_ID)/$(site_ID)_parameters.jl",
    ),
)

# Set up the model domain
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
prognostic_land_components = (:canopy, :soil, :soilco2)

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
    earth_param_set;
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
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(; domain = land_domain, drivers)

# Canopy model
# Set up radiative transfer
radiation_parameters =
    (; Ω, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf)
radiative_transfer =
    Canopy.TwoStreamModel{FT}(canopy_domain; radiation_parameters, ϵ_canopy)

# Set up conductance
conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain; g1)

# Set up photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis =
    FarquharModel{FT}(canopy_domain; photosynthesis_parameters, sc, pc)

# Set up plant hydraulics
(; LAI, maxLAI) =
    FluxnetSimulationsExt.prescribed_LAI_fluxnet(site_ID, start_date)
RAI = FT(maxLAI) * f_root_to_shoot # convert to float type of simulation
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
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land);
jacobian! = make_jacobian(land);
set_initial_cache! = make_set_initial_cache(land)
Y, p, cds = initialize(land)
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

FluxnetSimulationsExt.set_fluxnet_ic!(Y, site_ID, start_date, time_offset, land)
set_initial_cache!(p, Y, t0)

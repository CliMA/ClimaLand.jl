import SciMLBase
using Plots
using Statistics
using Dates
using Insolation
using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using StaticArrays
using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP


const FT = Float32;
earth_param_set = LP.LandParameters(FT);

nelements = 10
zmin = FT(-2)
zmax = FT(0)
f_root_to_shoot = FT(3.5)
SAI = FT(0.00242)
maxLAI = FT(4.2)
plant_ν = FT(2.46e-4) # kg/m^2
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9)
h_leaf = FT(9.5)
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [zmax, h_stem, h_stem + h_leaf]

include(
    joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/data_tools.jl"),
);

# First provide some information about the site
# Timezone (offset from UTC in hrs)
time_offset = 7

# Site latitude and longitude
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree
land_domain = Point(; z_sfc = FT(0.0))

# Height of the sensor at the site
atmos_h = FT(32)

# Provide the site site ID and the path to the data file:
site_ID = "US-MOz"
data_link = "https://caltech.box.com/shared/static/7r0ci9pacsnwyo0o9c25mhhcjhsu6d72.csv"

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
);

# Populate the SharedCanopyParameters struct, which holds the parameters
# shared between all different components of the canopy model.
z0_m = FT(2)
z0_b = FT(0.2)

shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
);

# For this canopy, we are running in standalone mode, which means we need to
# use a prescribed soil driver, defined as follows:

ψ_soil0 = FT(0.0)

soil_driver = PrescribedGroundConditions(
    FT;
    root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
    ψ = t -> ψ_soil0,
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    T = t -> 298.0,
    ϵ = FT(0.99),
);

# Now, setup the canopy model by component.
# Provide arguments to each component, beginning with radiative transfer:

rt_params = TwoStreamParameters(
    FT;
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    α_NIR_leaf = FT(0.45),
    τ_PAR_leaf = FT(0.05),
    τ_NIR_leaf = FT(0.25),
    Ω = FT(0.69),
    λ_γ_PAR = FT(5e-7),
)

rt_model = TwoStreamModel{FT}(rt_params);

# Arguments for conductance model:
cond_params = PModelConductanceParameters(Drel = FT(1.6))

stomatal_model = PModelConductance{FT}(cond_params);

# Arguments for photosynthesis model:
photo_params = PModelParameters(
    cstar = FT(0.41), 
    β = FT(146), 
    ϕc = FT(0.087),
    ϕ0 = FT(NaN), 
    ϕa0 = FT(0.352),
    ϕa1 = FT(0.022),
    ϕa2 = FT(-0.00034),
    α = FT(0.067),
    sc = FT(2e-6),
    pc = FT(-2e6)
)

photosynthesis_model = PModel{FT}(photo_params);

# Arguments for autotrophic respiration model:
AR_params = AutotrophicRespirationParameters(FT)
AR_model = AutotrophicRespirationModel{FT}(AR_params);

# Arguments for plant hydraulics model are more complicated.

# Begin by providing general plant parameters. For the area
# indices of the canopy, we choose a `PrescribedSiteAreaIndex`,
# which supports LAI as a function of time, with RAI and SAI as constant.
LAI = 4.2
LAIfunction = (t) -> LAI
SAI = FT(0.00242)
f_root_to_shoot = FT(3.5)
RAI = FT((SAI + LAI) * f_root_to_shoot)
ai_parameterization =
    PrescribedSiteAreaIndex{FT}(TimeVaryingInput(LAIfunction), SAI, RAI)
rooting_depth = FT(1.0);


# Create the component conductivity and retention models of the hydraulics
# model. In ClimaLand, a Weibull parameterization is used for the conductivity as
# a function of potential, and a linear retention curve is used.

K_sat_plant = FT(1.8e-8)
ψ63 = FT(-4 / 0.0098)
Weibull_param = FT(4)
a = FT(0.05 * 0.0098)

conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)

retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a);

# Use these values to populate the parameters of the PlantHydraulics model:

ν = FT(0.7)
S_s = FT(1e-2 * 0.0098)

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = ν,
    S_s = S_s,
    rooting_depth = rooting_depth,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
);

# Define the remaining variables required for the plant hydraulics model.

plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_surfaces = compartment_surfaces,
    compartment_midpoints = compartment_midpoints,
);

# Now, instantiate the canopy model, using the atmospheric and radiative
# drivers included from the external file, as well as the soil driver we
# instantiated above. This contains every piece of information needed to
# generate the set of ODEs modeling the canopy biophysics, ready to be passed
# off to a timestepper.

canopy = ClimaLand.Canopy.CanopyModel{FT}(;
    parameters = shared_params,
    domain = land_domain,
    autotrophic_respiration = AR_model,
    radiative_transfer = rt_model,
    photosynthesis = photosynthesis_model,
    conductance = stomatal_model,
    hydraulics = plant_hydraulics,
    boundary_conditions = Canopy.AtmosDrivenCanopyBC(
        atmos,
        radiation,
        soil_driver,
    ),
);

# Initialize the state vectors and obtain the model coordinates, then get the
# explicit time stepping tendency that updates auxiliary and prognostic
# variables that are stepped explicitly.

Y, p, coords = ClimaLand.initialize(canopy)
exp_tendency! = make_exp_tendency(canopy);
imp_tendency! = make_imp_tendency(canopy)
jacobian! = make_jacobian(canopy);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

# Provide initial conditions for the canopy hydraulics model

ψ_stem_0 = FT(-1e5 / 9800)
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        ν,
        S_s,
    )

for i in 1:2
    Y.canopy.hydraulics.ϑ_l.:($i) .= augmented_liquid_fraction.(ν, S_l_ini[i])
end;

# Run the simulation for 10 days with a timestep of 10 minutes 
t0 = FT(0.0)
N_days = 10
tf = t0 + 3600 * 24 * N_days
dt = FT(600.0);
set_initial_cache! = make_set_initial_cache(canopy)
set_initial_cache!(p, Y, t0);

# Save every 16 timesteps 
n = 16
saveat = Array(t0:(n * dt):tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat);

# Create the callback function which updates the forcing variables,
# or drivers.
updateat = Array(t0:1800:tf)
model_drivers = ClimaLand.get_drivers(canopy)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

# p model specific callback. Eventually we will need to make this automatically applied
# if we are using the Pmodel 
pmodel_cb = ClimaLand.make_PModel_callback(FT, UTC_DATETIME[1], dt, canopy, long)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb, pmodel_cb);


# Select a timestepping algorithm and setup the ODE problem.
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

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

# Now, we can solve the problem and store the model data in the saveat array,
# using [`SciMLBase.jl`](https://github.com/SciML/SciMLBase.jl) and
# [`ClimaTimeSteppers.jl`](https://github.com/CliMA/ClimaTimeSteppers.jl).

sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat);

pmodel_vars = [
    "canopy.photosynthesis.GPP",
    "canopy.photosynthesis.OptVars.Vcmax25_opt", 
    "canopy.photosynthesis.OptVars.Jmax25_opt", 
    "canopy.photosynthesis.OptVars.ξ_opt", 
    "canopy.photosynthesis.IntVars.ci",
    "canopy.photosynthesis.Jmax",
    "canopy.photosynthesis.J",
    "canopy.conductance.r_stomata_canopy"
]


# get the pmodel variables and drivers
extracted_vars = extract_variables(sv, pmodel_vars)

# Create time array for plotting
daily = sol.t ./ 3600 ./ 24

# Create plots for all extracted variables
plots_array = []

# GPP plot (convert to μmol/mol)
plt_gpp = Plots.plot(
    daily, 
    extracted_vars.GPP .* 1e6,
    title = "GPP",
    xlabel = "Days",
    ylabel = "GPP",
    linewidth = 2
)
push!(plots_array, plt_gpp)

# Vcmax25_opt plot
plt_vcmax = Plots.plot(
    daily, 
    extracted_vars.Vcmax25_opt,
    title = "Optimal Vcmax25",
    xlabel = "Days", 
    ylabel = "Vcmax25",
    linewidth = 2
)
push!(plots_array, plt_vcmax)

# Jmax25_opt plot
plt_jmax = Plots.plot(
    daily, 
    extracted_vars.Jmax25_opt,
    title = "Optimal Jmax25",
    xlabel = "Days",
    ylabel = "Jmax25", 
    linewidth = 2
)
push!(plots_array, plt_jmax)

# ξ_opt plot
plt_xi = Plots.plot(
    daily, 
    extracted_vars.ξ_opt,
    title = "Optimal ξ",
    xlabel = "Days",
    ylabel = "ξ",
    linewidth = 2
)
push!(plots_array, plt_xi)

# ci plot
plt_ci = Plots.plot(
    daily, 
    extracted_vars.ci,
    title = "Internal CO2 Concentration",
    xlabel = "Days",
    ylabel = "ci [Pa]",
    linewidth = 2
)
push!(plots_array, plt_ci)

# Jmax plot  
plt_jmax_current = Plots.plot(
    daily,
    extracted_vars.Jmax,
    title = "instantaneous Jmax",
    xlabel = "Days",
    ylabel = "Jmax",
    linewidth = 2
)
push!(plots_array, plt_jmax_current)

# J plot
plt_j = Plots.plot(
    daily,
    extracted_vars.J,
    title = "instantaneous electron transport rate",
    xlabel = "Days", 
    ylabel = "J",
    label = "Model",
    linewidth = 2
)
push!(plots_array, plt_j)

# Stomatal resistance plot
plt_resistance = Plots.plot(
    daily,
    extracted_vars.r_stomata_canopy,
    title = "Stomatal resistance",
    xlabel = "Days",
    ylabel = "Resistance [s/m]",
    linewidth = 2
)
push!(plots_array, plt_resistance)

# Combine all plots in a grid layout
final_plot = Plots.plot(plots_array..., layout = (4, 2), size = (1200, 1600))

# Save the output
savefig("ozark_standalone_canopy_test.png")

# Display the plot
display(final_plot)
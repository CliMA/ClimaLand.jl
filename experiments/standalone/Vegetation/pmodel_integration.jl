"""
This script runs a standalone canopy model with the P-model for photosynthesis and stomatal
conductance with prescribed atmosphere and soil. 
"""

import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Statistics
using Dates
using Insolation
using StaticArrays
using StatsBase
using NCDatasets

using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaUtilities.OutputPathGenerator: generate_output_path
using ClimaDiagnostics
using ClimaUtilities
const FT = Float64
earth_param_set = LP.LandParameters(FT)
climaland_dir = pkgdir(ClimaLand)

save_outputs = true

include(joinpath(climaland_dir, "experiments/integrated/fluxnet/data_tools.jl"))
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/plot_utils.jl"))


# Read in the site to be run from the command line
if length(ARGS) < 1
    error("Must provide site ID on command line")
end

site_ID = ARGS[1]
photo_model = ARGS[2] 
@assert photo_model in ("pmodel", "farquhar") "Photo model must be either 'pmodel' or 'farquhar'"

# Read all site-specific domain parameters from the simulation file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_simulation.jl",
    ),
)
h_canopy = h_stem + h_leaf

# Read all site-specific parameters from the parameter file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_parameters.jl",
    ),
)

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/fluxnet_simulation.jl",
    ),
)

include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
)

# Prescribed soil conditions
ψ_soil0 = FT(0.0)
soil_driver = PrescribedGroundConditions(
    FT;
    root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
    ψ = t -> ψ_soil0,
    α_PAR = soil_α_PAR,
    α_NIR = soil_α_NIR,
    T = t -> 298.0,
    ϵ = FT(0.99),
);

# Populate the SharedCanopyParameters struct, which holds the parameters
# shared between all different components of the canopy model.
z0_m = FT(2) # m, roughness length for momentum
z0_b = FT(0.2) # m, roughness length for scalars 
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
);

# Set the radiative transfer model
rt_params = TwoStreamParameters(
    FT;
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = α_PAR_leaf,
    α_NIR_leaf = α_NIR_leaf,
    τ_PAR_leaf = τ_PAR_leaf,
    τ_NIR_leaf = τ_NIR_leaf,
    Ω = Ω,
    λ_γ_PAR = λ_γ_PAR,
)

rt_model = TwoStreamModel{FT}(rt_params);

if photo_model == "pmodel"
    # Set the conductance model 
    cond_params = PModelConductanceParameters(Drel = FT(1.6))
    stomatal_model = PModelConductance{FT}(cond_params);

    # Set the photosynthesis model (P-model currently only supports C3) 
    photo_params = PModelParameters(
        cstar = FT(0.41),
        β = FT(146),
        ϕc = FT(0.087),
        ϕ0 = FT(NaN),
        ϕa0 = FT(0.352),
        ϕa1 = FT(0.022),
        ϕa2 = FT(-0.00034),
        α = FT(0.933),
        sc = FT(2e-6),
        pc = FT(-2e6),
    )
    photosynthesis_model = PModel{FT}(photo_params);
else
    # Set the conductance model 
    cond_params = MedlynConductanceParameters(FT; g1 = g1, g0 = g0, Drel = Drel)
    stomatal_model = MedlynConductanceModel{FT}(cond_params);

    # Set the photosynthesis model 
    is_c3 = FT(1) # set the photosynthesis mechanism to C3
    photo_params = FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25)
    photosynthesis_model = FarquharModel{FT}(photo_params);
end

# Set the autotrophic respiration model
AR_params = AutotrophicRespirationParameters(FT)
AR_model = AutotrophicRespirationModel{FT}(AR_params);

# Set the plant hydraulics model
# Begin by providing general plant parameters. For the area indices of the canopy, we choose a `PrescribedSiteAreaIndex`,
# which supports LAI as a function of time, with RAI and SAI as constant.

# Note: SAI is defined in parameters.jl for each site. maxLAI is calculated when we read in the MODIS LAI data in met_drivers_FLUXNET.jl
# LAIfunction is also defined in met_drivers_FLUXNET.jl
RAI = FT((SAI + maxLAI) * f_root_to_shoot) 
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a);
plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = FT(0.7),
    S_s = FT(1e-2 * 0.0098),
    rooting_depth = rooting_depth,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
);

# Define the remaining variables required for the plant hydraulics model.
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [zmax, h_stem, h_stem + h_leaf]
plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
    parameters = plant_hydraulics_ps,
    n_stem = 1,
    n_leaf = 1,
    compartment_surfaces = compartment_surfaces,
    compartment_midpoints = compartment_midpoints,
);

# instantiate the canopy model with all the components defined above.
canopy_domain = Point(; z_sfc = FT(0.0))
canopy = ClimaLand.Canopy.CanopyModel{FT}(;
    parameters = shared_params,
    domain = canopy_domain,
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
        soil_ν,
        plant_S_s,
    )

for i in 1:2
    Y.canopy.hydraulics.ϑ_l.:($i) .= augmented_liquid_fraction.(soil_ν, S_l_ini[i])
end;

# Run the simulation for 364 days with a timestep of 10 minutes 
# Note that this overrides the different timesteps specified in the {site}_simulation.jl files 
t0 = 0.0
N_days = 364
tf = t0 + 3600 * 24 * N_days
dt = 600.0;
set_initial_cache! = make_set_initial_cache(canopy)
set_initial_cache!(p, Y, t0);

@info "Using timestep of $(dt) seconds"
# outdir = joinpath(pkgdir(ClimaLand), "outputs")
# output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(outdir)
# start_date = DateTime(2010) # 2010-01-01T:00:00:00

# nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
#     canopy_domain.space.surface, 
#     output_dir;
#     start_date,
# )

# diags = ClimaLand.default_diagnostics(
#     canopy,
#     start_date;
#     output_writer = d_writer,
#     output_vars = :long,
#     average_period = :hourly,
# )

# diagnostic_handler =
#     ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0, dt = dt);

# diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler);
# Save every 30 mins 
save_period = 1800.0
saveat = Array(t0:save_period:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat);

# Create the callback function which updates the forcing variables,
# or drivers.
updateat = Array(t0:DATA_DT:tf)
model_drivers = ClimaLand.get_drivers(canopy)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

if photo_model == "pmodel" 
    pmodel_cb = ClimaLand.make_PModel_callback(FT, UTC_DATETIME[1], dt, canopy, long)
    cb = SciMLBase.CallbackSet(saving_cb, driver_cb, pmodel_cb)
else
    cb = SciMLBase.CallbackSet(saving_cb, driver_cb)
end

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
    "canopy.photosynthesis.IntVars.Γstar",
    "canopy.photosynthesis.IntVars.Kmm",
    "canopy.photosynthesis.Jmax",
    "canopy.photosynthesis.J",
    "canopy.photosynthesis.Vcmax",
    "canopy.photosynthesis.Ac",
    "canopy.photosynthesis.Aj",
    "canopy.photosynthesis.An",
    "canopy.photosynthesis.Rd",
    "canopy.conductance.r_stomata_canopy",
    "canopy.radiative_transfer.par.abs",
    "canopy.radiative_transfer.par_d",
    "canopy.hydraulics.area_index.leaf",
]

farquhar_vars = [
    "canopy.photosynthesis.GPP",
    "canopy.photosynthesis.Rd",
    "canopy.photosynthesis.An", 
    "canopy.conductance.r_stomata_canopy",
    "canopy.hydraulics.area_index.leaf"
]

if photo_model == "pmodel"
    extracted_vars = extract_variables(sv, pmodel_vars)
else
    extracted_vars = extract_variables(sv, farquhar_vars)
end

if save_outputs
    # Save outputs to NetCDF files
    # Convert UTC time to local time and create time array
    # Use sol.t which contains the actual simulation times in seconds
    local_datetime =
        UTC_DATETIME[1] - Dates.Hour(time_offset) .+ Dates.Second.(sol.t)
    time_seconds = [Dates.datetime2unix(dt) for dt in local_datetime]

    # Skip the first timestep for canopy integration outputs
    local_datetime_skip_first = local_datetime[2:end]
    time_seconds_skip_first = time_seconds[2:end]

    # Remove existing NetCDF files if they exist
    canopy_file = "outputs/$(site_ID)_canopy_integration_outputs_$(photo_model).nc"
    drivers_file = "outputs/$(site_ID)_fluxnet_drivers.nc"
    if isfile(canopy_file)
        rm(canopy_file)
    end
    if isfile(drivers_file)
        rm(drivers_file)
    end

    # Create canopy integration outputs NetCDF file
    NCDatasets.Dataset(canopy_file, "c") do ds
        NCDatasets.defDim(ds, "time", length(time_seconds_skip_first))
        time_var = NCDatasets.defVar(ds, "time", Float64, ("time",))
        time_var[:] = time_seconds_skip_first

        # Save each extracted variable (excluding first timestep)
        for (var_name, var_data) in pairs(extracted_vars)
            var_clean_name = replace(string(var_name), "." => "_")
            nc_var = NCDatasets.defVar(ds, var_clean_name, Float64, ("time",))
            nc_var[:] = var_data[2:end]
        end
    end

    # Create fluxnet drivers NetCDF file
    NCDatasets.Dataset(drivers_file, "c") do ds
        NCDatasets.defDim(ds, "time", length(time_seconds_skip_first))

        # Create time variable
        time_var = NCDatasets.defVar(ds, "time", Float64, ("time",))
        time_var[:] = time_seconds_skip_first

        # Save each driver variable
        driver_names = fieldnames(typeof(drivers))
        for driver_name in driver_names
            driver_data = getfield(drivers, driver_name)
            # Sample the driver data at the same intervals as the model output
            driver_indices = Int.(round.(sol.t ./ DATA_DT)) .+ 1
            # Ensure indices are within bounds
            driver_indices =
                clamp.(driver_indices, 1, length(driver_data.values))
            var_data = driver_data.values[driver_indices]
            # Skip last timestep to match canopy outputs length
            var_data = var_data[1:(end - 1)]

            nc_var =
                NCDatasets.defVar(ds, string(driver_name), Float64, ("time",))
            nc_var[:] = var_data
        end
    end
end

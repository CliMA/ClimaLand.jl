# # Timestep test for integrated SoilCanopyModel

# This experiment runs the standalone canopy model for 12.5 hours on a spherical domain
# and outputs useful plots related to the timestepping and stability of
# the simulation using different timesteps. This information can be used to
# determine the optimal timestep for a this simulation setup.

# The simulation is run multiple times: first, with a small reference timestep,
# then with increasing timestep sizes. The goal here is to see which timesteps
# are small enough to produce a stable simulation, and which timesteps are too
# large to produce accurate results.

# The runs with larger timesteps are compared with the reference timestep run
# using multiple metrics, including:
# - mean, 99%th, and 95th% error over the course of the simulation
# - T at the end of the simulation (used to test convergence with each dt)
# - T throughout the entire simulation

# Note that the simulation is run for 12.5 hours to allow us to test the
# convergence behavior of the solver.
# Convergence behavior must be assessed when T is actively changing, rather
# than when it is in equilibrium. T changes in the time period after a driver
# update, which occurs every 3 hours (so exactly at 12 hours). Running for
# just longer than 12 hours allows us to observe convergence behavior at a time
# when T is actively changing. See the discussion in
# https://github.com/CliMA/ClimaLand.jl/pull/675 for more information.

# Simulation Setup
# Number of spatial elements: 10 in horizontal, 5 in vertical
# Soil depth: 5 m
# Simulation duration: 12.5 hours
# Timestep: variable, ranging from 0.6s to 3600s
# Timestepper: ARS111
# Fixed number of iterations: 6
# Jacobian update: Every Newton iteration
# Atmospheric and radiation driver updates: Every 3 hours

import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
using Dates
using Insolation
using Statistics
import StatsBase: percentile
using CairoMakie
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

using ClimaLand
using ClimaLand.Domains
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams


PROFILING = false
try
    import Profile, ProfileCanvas
    global PROFILING = true
    @info "ProfileCanvas found, running with profiler"
catch
end

function set_initial_conditions(land, t0)
    Y, p, cds = initialize(land)
    FT = eltype(Y.soil.ϑ_l)
    set_initial_cache! = make_set_initial_cache(land)

    Y.soil.ϑ_l = FT(0.3)
    Y.soil.θ_i = FT(0.0)
    T_0 = FT(297.5)
    ρc_s =
        volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            land.soil.parameters.ρc_ds,
            land.soil.parameters.earth_param_set,
        )
    Y.soil.ρe_int =
        volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T_0,
            land.soil.parameters.earth_param_set,
        )

    Y.soilco2.C = FT(0.000412) # set to atmospheric co2, mol co2 per mol air

    ψ_stem_0 = FT(-1e5 / 9800)
    ψ_leaf_0 = FT(-2e5 / 9800)
    canopy_params = land.canopy.hydraulics.parameters
    S_l_ini =
        inverse_water_retention_curve.(
            canopy_params.retention_model,
            [ψ_stem_0, ψ_leaf_0],
            canopy_params.ν,
            canopy_params.S_s,
        )

    for i in 1:2
        S_l = S_l_ini[i]
        Y.canopy.hydraulics.ϑ_l.:($i) .=
            augmented_liquid_fraction.(canopy_params.ν, S_l)
    end

    Y.canopy.energy.T = FT(297.5)
    set_initial_cache!(p, Y, t0)
    return Y, p
end

context = ClimaComms.context()
device_suffix =
    typeof(ClimaComms.context().device) <: ClimaComms.CPUSingleThreaded ?
    "cpu" : "gpu"
const FT = Float64
earth_param_set = LP.LandParameters(FT)

# Set the model domain
dz_bottom = FT(2.0)
dz_top = FT(0.2)
soil_depth = FT(5)
z_sfc = FT(0)
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9) # m
h_leaf = FT(9.5) # m

h_canopy = h_stem + h_leaf
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [z_sfc, h_stem, h_canopy]

land_domain = ClimaLand.Domains.SphericalShell(;
    radius = FT(6.3781e6),
    depth = soil_depth,
    nelements = (10, 5),
    npolynomial = 1,
    dz_tuple = FT.((dz_bottom, dz_top)),
);
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
sfc_cds = ClimaCore.Fields.coordinate_field(land_domain.space.surface)

# First pick the start date and time of the simulation, since time varying input depends on that
t0 = Float64(0) # start at start date
start_date = DateTime("202001010000", "yyyymmddHHMM")
# Time varying input
LAIfunction = TimeVaryingInput((t) -> 2.0)
# Atmospheric and radiative forcing
precip = TimeVaryingInput((t) -> -1.0e-7)
atmos_q = TimeVaryingInput((t) -> 0.002)
atmos_T = TimeVaryingInput((t) -> 298.0)
atmos_p = TimeVaryingInput((t) -> 101320)
atmos_u = TimeVaryingInput((t) -> 3.0)
LW_IN = TimeVaryingInput((t) -> 5.67e-8 * 298.0^4)
SW_IN = TimeVaryingInput((t) -> 500.0)
snow_precip = TimeVaryingInput((t) -> 0.0)
atmos_h = FT(32)
# Construct the drivers. For the start date we will use the UTC time at the
# start of the simulation
atmos = ClimaLand.PrescribedAtmosphere(
    precip,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    start_date,
    atmos_h,
    earth_param_set;
)
function zenith_angle(
    t,
    start_date;
    cd_field = sfc_cds,
    insol_params::Insolation.Parameters.InsolationParameters{_FT} = earth_param_set.insol_params,
) where {_FT}
    # This should be time in UTC
    current_datetime = start_date + Dates.Second(round(t))

    # Orbital Data uses Float64, so we need to convert to our sim FT
    d, δ, η_UTC =
        _FT.(
            Insolation.helper_instantaneous_zenith_angle(
                current_datetime,
                start_date,
                insol_params,
            )
        )

    Insolation.instantaneous_zenith_angle.(
        d,
        δ,
        η_UTC,
        sfc_cds.long,
        sfc_cds.lat,
    ).:1
end

radiation = ClimaLand.PrescribedRadiativeFluxes(
    FT,
    SW_IN,
    LW_IN,
    start_date,
    θs = zenith_angle,
)

# Model parameters
# Soil parameters
soil_ν = FT(0.5) # m3/m3
soil_K_sat = FT(4e-7) # m/
soil_S_s = FT(1e-3) # 1/m
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
θ_r = FT(0.067)
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
# Energy Balance model
ac_canopy = FT(2.5e3)
# Conductance Model
g1 = FT(141)
#Photosynthesis model
Vcmax25 = FT(9e-5)
# Plant Hydraulics and general plant parameters
K_sat_plant = 5e-9
ψ63 = FT(-4 / 0.0098)
Weibull_param = FT(4)
a = FT(0.05 * 0.0098)
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
plant_ν = FT(1.44e-4)
SAI = FT(1.0)
maxLAI = FT(2.0)
f_root_to_shoot = FT(3.5)
RAI = maxLAI * f_root_to_shoot
capacity = plant_ν * maxLAI * h_leaf * FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m

# Set up model
# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
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
)

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

soilco2_ps = SoilCO2ModelParameters(FT)

# soil microbes args
Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

# Set the soil CO2 BC to being atmospheric CO2
soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> 0.0)
soilco2_sources = (MicrobeProduction{FT}(),)

soilco2_boundary_conditions = (; top = soilco2_top_bc, bottom = soilco2_bot_bc)

soilco2_args = (;
    boundary_conditions = soilco2_boundary_conditions,
    sources = soilco2_sources,
    domain = soil_domain,
    parameters = soilco2_ps,
)

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
radiative_transfer_args = (; parameters = TwoStreamParameters(FT))
# Set up conductance
conductance_args = (; parameters = MedlynConductanceParameters(FT; g1))
# Set up photosynthesis
is_c3 = FT(1) # set the photosynthesis mechanism to C3
photosynthesis_args =
    (; parameters = FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)
plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    rooting_depth = FT(0.5),
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
land_input = (
    atmos = atmos,
    radiation = radiation,
    soil_organic_carbon = Csom,
    runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
)
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

# Define explicit and implicit tendencies, and the jacobian
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land);
jacobian! = make_jacobian(land);

# Timestepping information
N_hours = 8
tf = Float64(t0 + N_hours * 3600.0)
sim_time = round((tf - t0) / 3600, digits = 2) # simulation length in hours

timestepper = CTS.ARS111()
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

# Set up simulation callbacks
drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(drivers)

# Choose timesteps and set up arrays to store information
ref_dt = 50.0
dts = [ref_dt, 100.0, 200.0, 450.0, 900.0, 1800.0, 3600.0]

ref_T = []
mean_err = []
p95_err = []
p99_err = []
sol_err = []
T_states = []
times = []
ΔT = FT(0)
for dt in dts
    @info dt
    # Initialize model and set initial conditions before each simulation
    Y, p = set_initial_conditions(land, t0)
    jac_kwargs = (;
        jac_prototype = ClimaLand.ImplicitEquationJacobian(Y),
        Wfact = jacobian!,
    )

    # Create callback for driver updates
    saveat = vcat(Array(t0:(3 * 3600):tf), tf)
    updateat = vcat(Array(t0:(3 * 3600):tf), tf)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

    # Solve simulation
    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    @time sol = SciMLBase.solve(
        prob,
        ode_algo;
        dt = dt,
        callback = driver_cb,
        saveat = saveat,
    )

    # Save results for comparison
    if dt == ref_dt
        global ref_T =
            [parent(sol.u[k].canopy.energy.T)[1] for k in 1:length(sol.t)]
    else
        T = [parent(sol.u[k].canopy.energy.T)[1] for k in 1:length(sol.t)]

        ΔT = abs.(T .- ref_T)
        push!(mean_err, FT(mean(ΔT)))
        push!(p95_err, FT(percentile(ΔT, 95)))
        push!(p99_err, FT(percentile(ΔT, 99)))
        push!(sol_err, ΔT[end])
        push!(T_states, T)
        push!(times, sol.t)
    end
end

savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "integrated",
    "performance",
    "integrated_timestep_test",
);
!ispath(savedir) && mkpath(savedir)

# Create plot with statistics
# Compare T state computed with small vs large dt
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel = "Timestep (sec)",
    ylabel = "Temperature (K)",
    xscale = log10,
    yscale = log10,
    title = "Error in T over $(sim_time) hour sim, dts in [$(dts[1]), $(dts[end])]",
)
dts = dts[2:end]
lines!(ax, dts, FT.(mean_err), label = "Mean Error")
lines!(ax, dts, FT.(p95_err), label = "95th% Error")
lines!(ax, dts, FT.(p99_err), label = "99th% Error")
axislegend(ax)
save(joinpath(savedir, "errors.png"), fig)

# Create convergence plot
fig2 = Figure()
ax2 = Axis(
    fig2[1, 1],
    xlabel = "log(dt)",
    ylabel = "log(|T[end] - T_ref[end]|)",
    xscale = log10,
    yscale = log10,
    title = "Convergence of T; sim time = $(sim_time) hours, dts in [$(dts[1]), $(dts[end])]",
)
scatter!(ax2, dts, FT.(sol_err))
lines!(ax2, dts, dts / 1e3)
save(joinpath(savedir, "convergence.png"), fig2)

# Plot T throughout full simulation for each run
fig3 = Figure()
ax3 = Axis(
    fig3[1, 1],
    xlabel = "time (hr)",
    ylabel = "T (K)",
    title = "T throughout simulation; length = $(sim_time) hours, dts in [$(dts[1]), $(dts[end])]",
)
times = times ./ 3600.0 # hours
for i in 1:length(times)
    lines!(ax3, times[i], T_states[i], label = "dt $(dts[i]) s")
end
axislegend(ax3, position = :rt)
save(joinpath(savedir, "states.png"), fig3)

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
ClimaComms.@import_required_backends
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.OutputPathGenerator: generate_output_path
import ClimaUtilities.TimeManager: date
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaParams


PROFILING = false
try
    import Profile, ProfileCanvas
    global PROFILING = true
    @info "ProfileCanvas found, running with profiler"
catch
end

function set_ic!(Y, p, t0, model)
    FT = eltype(Y.soil.ϑ_l)
    Y.soil.ϑ_l = FT(0.3)
    Y.soil.θ_i = FT(0.0)
    T_0 = FT(297.5)
    ρc_s =
        volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            model.soil.parameters.ρc_ds,
            model.soil.parameters.earth_param_set,
        )
    Y.soil.ρe_int =
        volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T_0,
            model.soil.parameters.earth_param_set,
        )

    Y.soilco2.C = FT(0.000412) # set to atmospheric co2, mol co2 per mol air

    ψ_stem_0 = FT(-1e5 / 9800)
    ψ_leaf_0 = FT(-2e5 / 9800)
    canopy_params = model.canopy.hydraulics.parameters
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
    return
end

context = ClimaComms.context()
ClimaComms.init(context)
device_suffix =
    typeof(context.device) <: ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
const FT = Float64
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
prognostic_land_components = (:canopy, :soil, :soilco2)

# Set the model domain
dz_bottom = FT(2.0)
dz_top = FT(0.2)
soil_depth = FT(5)
z_sfc = FT(0)

land_domain = ClimaLand.Domains.SphericalShell(;
    radius = FT(6.3781e6),
    depth = soil_depth,
    nelements = (10, 5),
    dz_tuple = FT.((dz_bottom, dz_top)),
);
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
sfc_cds = ClimaCore.Fields.coordinate_field(land_domain.space.surface)

# First pick the start date and time of the simulation, since time varying input depends on that
t0 = Float64(0) # start at start date
start_date = DateTime("202001010000", "yyyymmddHHMM")

# Atmospheric and radiative forcing
precip = TimeVaryingInput((t) -> -1.0e-7)
atmos_q = TimeVaryingInput((t) -> 0.002)
atmos_T = TimeVaryingInput((t) -> 298.0)
atmos_p = TimeVaryingInput((t) -> 101320)
atmos_u = TimeVaryingInput((t) -> 3.0)
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
    toml_dict,
);
LW_IN = TimeVaryingInput((t) -> 5.67e-8 * 298.0^4)
SW_IN = TimeVaryingInput((t) -> 500.0)
insol_params = earth_param_set.insol_params # parameters of Earth's orbit required to compute the insolation
coords = ClimaCore.Fields.coordinate_field(land_domain.space.surface)
zenith_angle =
    (t, s) ->
        default_zenith_angle.(
            t,
            s;
            insol_params,
            longitude = coords.long,
            latitude = coords.lat,
        );
radiation = ClimaLand.PrescribedRadiativeFluxes(
    FT,
    SW_IN,
    LW_IN,
    start_date,
    θs = zenith_angle,
    toml_dict = toml_dict,
);

# Soil model setup
# Soil parameters
soil_ν = FT(0.5) # m3/m3
soil_K_sat = FT(4e-7) # m/
soil_S_s = FT(1e-3)
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
θ_r = FT(0.067)
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);

soil_forcing = (; atmos, radiation)
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = FT(0.2),
    NIR_albedo = FT(0.2),
)
runoff = ClimaLand.Soil.Runoff.SurfaceRunoff()
retention_parameters = (;
    ν = soil_ν,
    θ_r,
    K_sat = soil_K_sat,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
)
composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)
soil = Soil.EnergyHydrology{FT}(
    land_domain,
    soil_forcing,
    toml_dict;
    prognostic_land_components,
    albedo = soil_albedo,
    runoff,
    retention_parameters,
    composition_parameters,
    S_s = soil_S_s,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
)

# Soil microbes model setup
soil_organic_carbon =
    ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers = Soil.Biogeochemistry.SoilDrivers(
    co2_prognostic_soil,
    soil_organic_carbon,
    atmos,
)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, drivers, toml_dict)

# Canopy model setup
# Radiative transfer model
G_Function = ConstantGFunction(ClimaParams.float_type(toml_dict)(0.5))
α_PAR_leaf = 0.3
τ_PAR_leaf = 0.2
α_NIR_leaf = 0.4
τ_NIR_leaf = 0.25
Ω = 1
radiation_parameters =
    (; G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf, Ω)
radiative_transfer = Canopy.TwoStreamModel{FT}(
    Canopy.TwoStreamParameters(toml_dict; radiation_parameters...),
)

# Photosynthesis model
Vcmax25 = FT(9e-5)
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis =
    FarquharModel{FT}(canopy_domain, toml_dict; photosynthesis_parameters)

# Conductance model
g1 = FT(141)
conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain, toml_dict; g1)

# Hydraulics model
LAI = TimeVaryingInput((t) -> 2.0)
SAI = FT(1.0)
maxLAI = FT(2.0)
f_root_to_shoot = FT(3.5)
RAI = maxLAI * f_root_to_shoot
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
conductivity_model = PlantHydraulics.Weibull(toml_dict)
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9) # m
h_leaf = FT(9.5) # m
h_canopy = h_stem + h_leaf
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain,
    LAI,
    toml_dict;
    SAI,
    RAI,
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    rooting_depth = FT(0.5),
)

# Put all the components together to form the canopy model
z_0m = FT(0.13) * h_canopy
z_0b = FT(0.1) * z_0m
canopy_domain = obtain_surface_domain(land_domain)
ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
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
    hydraulics,
)

# Integrated plant and soil model
land = SoilCanopyModel{FT}(soilco2, soil, canopy)

# Timestepping information
N_hours = 8
tf = Float64(t0 + N_hours * 3600.0)
stop_date = start_date + Dates.Second(round(tf))
sim_time = round((tf - t0) / 3600, digits = 2) # simulation length in hours

timestepper = CTS.ARS111()
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

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
for dt in dts
    @info dt
    # Initialize model and set initial conditions before each simulation
    simulation = LandSimulation(
        start_date,
        stop_date,
        dt,
        land;
        solver_kwargs = (;
            saveat = vcat(
                collect(start_date:Dates.Hour(3):stop_date),
                stop_date,
            )
        ),
        diagnostics = [],
        updateat = Dates.Hour(3),
        timestepper = ode_algo,
        set_ic!,
    )
    sol = ClimaLand.Simulations.solve!(simulation)
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
savedir = generate_output_path(
    "experiments/integrated/performance/integrated_timestep_test",
)

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
for i in 1:length(times)
    lines!(ax3, FT.(times[i]) ./ 3600.0, T_states[i], label = "dt $(dts[i]) s")
end
axislegend(ax3, position = :rt)
save(joinpath(savedir, "states.png"), fig3)

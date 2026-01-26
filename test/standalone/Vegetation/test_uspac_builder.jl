# Canopy Model Simulation Script
# Based on MC3Workshop tutorial: https://github.com/CliMA/MC3Workshop

using ClimaLand
using ClimaLand.Canopy
using ClimaCore
using Dates
import ClimaParams as CP
import ClimaTimeSteppers as CTS
import SciMLBase
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

println("="^60)
println("ClimaLand Canopy Model Simulation")
println("="^60)

# Simulation parameters
FT = Float64
t0 = Float64(0)
tf = Float64(86400 * 30)  # 30 days
dt = Float64(1800)         # 30 min timestep
saveat = Float64(3600)     # Save every hour

# Domain setup - single column
zmin = FT(-2.0)  # 2m soil depth
zmax = FT(0.0)
zlim = (zmin, zmax)
nelements = 10

domain = ClimaLand.Domains.Column(;
    zlim = zlim,
    nelements = nelements,
    dz_tuple = FT.((1.0, 0.05)),
)

println("\n✓ Domain created:")
println("  - Soil depth: $(abs(zmin)) m")
println("  - Elements: $nelements")

# Get surface space for spatially-varying parameters
surface_space = domain.space.surface

# Parameter setup
param_dict = CP.create_toml_dict(FT)
earth_params = ClimaLand.Parameters.LandParameters(param_dict)

# Soil parameters
soil_params = ClimaLand.Soil.EnergyHydrologyParameters(
    param_dict;
    ν = FT(0.43),
    ν_ss_om = FT(0.0),
    ν_ss_quartz = FT(0.7),
    ν_ss_gravel = FT(0.0),
    K_sat = FT(1e-5),
    S_s = FT(1e-3),
    θ_r = FT(0.045),
    hydrology_cm = ClimaLand.Soil.vanGenuchten{FT}(;
        α = FT(2.6),
        n = FT(1.56),
    ),
)

println("\n✓ Soil parameters set")

# Atmospheric forcing - use TimeVaryingInput wrappers
# Simple constant values with diurnal solar cycle
liquid_precip = TimeVaryingInput((t) -> 0.0)
snow_precip = TimeVaryingInput((t) -> 0.0)
atmos_T = TimeVaryingInput((t) -> 298.0)  # 25°C
atmos_u = TimeVaryingInput((t) -> 3.0)    # 3 m/s wind
atmos_q = TimeVaryingInput((t) -> 0.01)   # specific humidity
atmos_P = TimeVaryingInput((t) -> 101325.0)  # 1 atm
SW_d = TimeVaryingInput((t) -> 500.0 * max(0.0, sin(2π * t / 86400)))  # diurnal
LW_d = TimeVaryingInput((t) -> 400.0)
co2_concentration = TimeVaryingInput((t) -> 4.2e-4)  # CO2 concentration

start_date = DateTime(2021, 1, 1)

atmos = ClimaLand.PrescribedAtmosphere(
    liquid_precip,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_P,
    start_date,
    FT(15.0),          # canopy height (m)
    param_dict;        # parameter dictionary
    gustiness = FT(1.0),
    c_co2 = co2_concentration,
)

println("\n✓ Atmospheric forcing set (diurnal solar cycle)")

# Plant hydraulics parameters - use positional arguments for Weibull
conductivity_model = ClimaLand.Canopy.PlantHydraulics.Weibull{FT}(
    FT(2e-7),  # K_sat
    FT(-2.5),  # ψ63
    FT(3.0),   # c
)

plant_hydraulics_ps = ClimaLand.Canopy.PlantHydraulics.PlantHydraulicsParameters(;
    conductivity_model = conductivity_model,
    ν = FT(1.44e-4),
    S_s = FT(1e-2 * 0.0098),
    retention_model = ClimaLand.Canopy.PlantHydraulics.LinearRetentionCurve{FT}(
        FT(0.2 * 0.0098)
    ),
)

# Plant hydraulics model
plant_hydraulics = ClimaLand.Canopy.PlantHydraulics.PlantHydraulicsModel{FT}(;
    parameters = plant_hydraulics_ps,
    n_stem = 1,
    n_leaf = 1,
    compartment_surfaces = FT.([0, 10, 15]),
    compartment_midpoints = FT.([5, 12.5]),
)

println("\n✓ Plant hydraulics configured")
println("  - Rooting depth: 1.5 m")

# Stomatal conductance - Medlyn model
stomatal_g_params = ClimaLand.Canopy.MedlynConductanceParameters(
    param_dict; 
    g1 = FT(141)
)
stomatal_model = ClimaLand.Canopy.MedlynConductanceModel{FT}(stomatal_g_params)

println("\n✓ Stomatal conductance: Medlyn model (g1=141)")

# Photosynthesis - Farquhar model with spatially-varying parameters
# Create a C3 field (1.0 = C3, 0.0 = C4)
is_c3_field = ClimaCore.Fields.ones(surface_space)

photosynthesis_params = ClimaLand.Canopy.FarquharParameters(
    param_dict;
    is_c3 = is_c3_field,
    Vcmax25 = FT(50e-6),
)
photosynthesis = ClimaLand.Canopy.FarquharModel{FT}(photosynthesis_params)

println("\n✓ Photosynthesis: Farquhar C3 model")

# Radiative transfer - needs domain and param_dict
radiation = ClimaLand.Canopy.BeerLambertModel{FT}(domain, param_dict)

println("\n✓ Radiation: Beer-Lambert")

# Autotrophic respiration - simple model
resp_params = ClimaLand.Canopy.AutotrophicRespirationParameters(param_dict)
autotrophic_respiration = ClimaLand.Canopy.AutotrophicRespirationModel{FT}(resp_params)

println("\n✓ Autotrophic respiration: Default model")

# Shared canopy parameters
shared_params = ClimaLand.Canopy.SharedCanopyParameters{FT}(
    param_dict;
    z0_m = FT(2.0),
    z0_b = FT(0.2),
)

# Build the canopy model
canopy = ClimaLand.Canopy.CanopyModel{FT}(;
    parameters = shared_params,
    domain = domain,
    plant_hydraulics = plant_hydraulics,
    stomatal_conductance = stomatal_model,
    photosynthesis = photosynthesis,
    radiation = radiation,
    autotrophic_respiration = autotrophic_respiration,
    atmos = atmos,
)

println("\n✓ Canopy model assembled")

println("\n" * "="^60)
println("Initializing simulation...")
println("="^60)

# Initialize model state
Y, p, coords = ClimaLand.initialize(canopy)

# Set initial conditions for plant hydraulics
Y.canopy.hydraulics.ϑ_l .= FT(0.5)  # Initial plant water content

println("\n✓ Initial state set:")
println("  - Plant water content: 0.5")

# Create tendency functions
exp_tendency! = ClimaLand.make_exp_tendency(canopy)
imp_tendency! = ClimaLand.make_imp_tendency(canopy)

# Set up timestepper
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = exp_tendency!, T_imp! = imp_tendency!),
    Y,
    (t0, tf),
    p,
)

# Callback for progress updates
callback = SciMLBase.FunctionCallingCallback(
    (u, t, integrator) -> begin
        if t % (86400.0) < dt  # Print once per day
            day = Int(t / 86400)
            println("  Day $day / $(Int(tf/86400))")
        end
        return nothing
    end;
    func_everystep = true,
)

println("\n" * "="^60)
println("Running simulation...")
println("  Duration: $(tf/86400) days")
println("  Timestep: $(dt/60) minutes")
println("  Outputs: Every $(saveat/3600) hours")
println("="^60)

# Solve
sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = dt,
    saveat = saveat,
    callback = callback,
)

println("\n" * "="^60)
println("Simulation complete!")
println("="^60)
println("  Status: $(sol.retcode)")
println("  Timesteps: $(length(sol.t))")
println("  Final time: $(sol.t[end]/86400) days")

# Extract and display some results
final_state = sol.u[end]
println("\nFinal state:")
println("  - Plant water content: $(mean(final_state.canopy.hydraulics.ϑ_l))")

# Save results (optional)
# using JLD2
# @save "canopy_simulation_results.jld2" sol

println("\n✓ Done!")
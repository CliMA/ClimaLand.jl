# # Coupled heat and water equations tending towards equilibrium

# The [Richards equation tutorial](@ref Soil/richards_equation.md)
# demonstrates how to solve for water flow in soil, without considering
# heat transfer, phase changes, or the effect of temperature and the effect of
# ice on the hydraulic properties of the soil.

# Here we show how to solve the interacting heat and water equations,
# in sand, but without phase changes. This allows us to capture
# behavior that is not present in Richards equation alone.

# The equations
# are:

# ``
# \frac{∂ ρe_{int}}{∂ t} =  ∇ ⋅ κ(θ_l, θ_i; ν, ...) ∇T + ∇ ⋅ ρe_{int_{liq}} K (T,θ_l, θ_i; ν, ...) \nabla h( ϑ_l, z; ν, ...)
# ``

# ``
# \frac{ ∂ ϑ_l}{∂ t} = ∇ ⋅ K (T,θ_l, θ_i; ν, ...) ∇h( ϑ_l, z; ν, ...).
# ``

# Here

# ``t`` is the time (s),

# ``z`` is the location in the vertical (m),

# ``ρe_{int}`` is the volumetric internal energy of the soil (J/m^3),

# ``T`` is the temperature of the soil (K),

# ``κ`` is the thermal conductivity (W/m/K),

# ``ρe_{int_{liq}}`` is the volumetric internal energy of liquid water (J/m^3),

# ``K`` is the hydraulic conductivity (m/s),

# ``h`` is the hydraulic head (m),

# ``ϑ_l`` is the augmented volumetric liquid water fraction,

# ``θ_i`` is the volumetric ice fraction, and

# ``ν, ...`` denotes parameters relating to soil type, such as porosity.


# We will solve this equation in an effectively 1-d domain with ``z ∈ [-1,0]``,
# and with the following boundary and initial conditions:

# ``- κ ∇T(t, z = 0) = 0 ẑ``

# `` -κ ∇T(t, z = -1) = 0 ẑ ``

# `` T(t = 0, z) = T_{min} + (T_{max}-T_{min}) e^{Cz}``

# ``- K ∇h(t, z = 0) = 0 ẑ ``

# `` -K ∇h(t, z = -1) = 0 ẑ``

# `` ϑ(t = 0, z) = ϑ_{min} + (ϑ_{max}-ϑ_{min}) e^{Cz}, ``

# where ``C, T_{min}, T_{max}, ϑ_{min},`` and ``ϑ_{max}`` are
# constants.


# If we evolve this system for times long compared to the dynamical timescales
# of the system, we expect it to reach an equilibrium where
# the LHS of these equations tends to zero.
# Assuming zero fluxes at the boundaries, the resulting equilibrium state
# should satisfy ``∂h/∂z = 0`` and ``∂T/∂z = 0``. Physically, this means that
# the water settles into a vertical profile in which
# the resulting pressure balances gravity and that the temperature
# is constant across the domain.

#  We verify that the system is approaching this equilibrium, and we also sketch out
# an analytic calculation for the final temperature in equilibrium.

# # Import necessary modules
# External (non - CliMA) modules
import SciMLBase
using Statistics
using Plots

# CliMA packages and ClimaLand modules
using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using ClimaLand
using ClimaLand.Domains: Column
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand.Soil

import ClimaLand
import ClimaLand.Parameters as LP

# Choose a floating point precision, and get the parameter set, which holds constants used across CliMA models:
FT = Float32
earth_param_set = LP.LandParameters(FT);


# # Create the model
# Set the values of other parameters required by the model:
ν = FT(0.395);
# Soil solids
# are the components of soil besides water, ice, gases, and air.
# We specify the soil component fractions, relative to all soil solids.
# These do not sum to unity; the remainder is ν_ss_minerals (=0.08, in this case).
ν_ss_quartz = FT(0.92);
ν_ss_om = FT(0.0);
ν_ss_gravel = FT(0.0);
# Other parameters include the hydraulic conductivity at saturation, the specific
# storage, and the van Genuchten parameters for sand.
# We recommend Chapter 8 of Bonan (2019) for finding parameters
# for other soil types.
Ksat = FT(4.42 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.89)
vg_α = FT(7.5) # inverse meters
hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
θ_r = FT(0.0)
params = Soil.EnergyHydrologyParameters(
    FT;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm,
    K_sat = Ksat,
    S_s,
    θ_r,
);

# We also need to pick a domain on which to solve the equations:
zmax = FT(0)
zmin = FT(-1.0)
nelems = 50
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);


# The boundary value problem in this case
# requires a boundary condition at the top and the bottom of the domain
# for each equation being solved. We support conditions on the state (`ϑ_l`
# or `T`), or on the fluxes (`-K∇h` or `-κ∇T`). In the case of fluxes,
# we return the magnitude of the flux, assumed to point along `ẑ`. And, in each case,
# the boundary conditions are supplied in the form of a function of auxiliary variables
# `p` and time `t`.
#  Here we choose flux boundary conditions. The flux boundary condition
# requires a function of the cache and simulation time which returns
# the boundary flux.

# Water boundary conditions:
surface_water_flux = WaterFluxBC((p, t) -> 0.0)
bottom_water_flux = WaterFluxBC((p, t) -> 0.0);

# The boundary conditions for the heat equation:
surface_heat_flux = HeatFluxBC((p, t) -> 0.0)
bottom_heat_flux = HeatFluxBC((p, t) -> 0.0);

# We wrap up all of those in a WaterHeatBC struct:
boundary_fluxes = (;
    top = WaterHeatBC(; water = surface_water_flux, heat = surface_heat_flux),
    bottom = WaterHeatBC(; water = bottom_water_flux, heat = bottom_heat_flux),
);


# We aren't using any sources or sinks in the equations here, but this is where
# freeze/thaw terms, runoff, root extraction, etc. would go.
sources = ();

# Lastly, we can create the [`EnergyHydrology`](@ref
# ClimaLand.Soil.EnergyHydrology) model.
# As always, the
# model encodes and stores all of the information (parameters, continous equations,
# prognostic variables, etc) which are needed to turn the PDE system into a set of ODEs,
# properly spatially discretized for the domain of interest.
soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
);



function set_ic!(Y, p, t0, model)
    params = model.parameters
    z = model.domain.fields.z
    ν = params.ν
    θ_r = params.θ_r
    FT = eltype(Y.soil.ϑ_l)
    zmax = FT(0)
    zmin = FT(-1)

    theta_max = FT(ν * 0.5)
    theta_min = FT(ν * 0.4)
    T_max = FT(289.0)
    T_min = FT(288.0)

    c = FT(20.0)
    @. Y.soil.ϑ_l =
        theta_min +
        (theta_max - theta_min) * exp(-(z - zmax) / (zmin - zmax) * c)
    Y.soil.θ_i .= FT(0.0)

    T = @.(T_min + (T_max - T_min) * exp(-(z - zmax) / (zmin - zmax) * c))

    θ_l = Soil.volumetric_liquid_fraction.(Y.soil.ϑ_l, ν, θ_r)
    ρc_s =
        Soil.volumetric_heat_capacity.(
            θ_l,
            Y.soil.θ_i,
            params.ρc_ds,
            params.earth_param_set,
        )
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            params.earth_param_set,
        )
end

# We choose the initial and final simulation times:
t0 = Float64(0)
tf = Float64(60 * 60 * 72);

# We use [ClimaTimesteppers.jl](https://github.com/CliMA/ClimaTimesteppers.jl) for carrying out the time integration.

# Choose a timestepper and set up the ODE problem:
dt = Float64(1000.0);
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

# By default, it
# only returns Y and t at each time we request output (`saveat`, below). We use
# a callback in order to also get the auxiliary vector `p` back:
saveat = collect(t0:FT(30000):tf);
saved_values = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
);
cb = ClimaLand.NonInterpSavingCallback(saved_values, saveat);

simulation = LandSimulation(
    t0,
    tf,
    dt,
    soil;
    set_ic! = set_ic!,
    solver_kwargs = (; saveat = deepcopy(saveat)),
    timestepper = ode_algo,
    user_callbacks = (cb,),
    diagnostics = (),
);


# Now we can solve the problem.
sol = solve!(simulation);

# Extract output
z = parent(soil.domain.fields.z)
t = parent(FT.(sol.t))
ϑ_l = [parent(sol.u[k].soil.ϑ_l) for k in 1:length(t)]
T = [parent(saved_values.saveval[k].soil.T) for k in 1:length(t)];
# Let's look at the initial and final times:
plot(
    ϑ_l[1],
    z,
    xlabel = "ϑ_l",
    ylabel = "z (m)",
    label = "t = 0d",
    title = "Moisture Equilibration from t = 0d to t = 3d",
)
plot!(ϑ_l[4], z, label = "t = 1.5d")
plot!(ϑ_l[end], z, label = "t = 3d")
savefig("eq_moisture_plot.png");
# ![](eq_moisture_plot.png)

plot(
    T[1],
    z,
    xlabel = "T (K)",
    ylabel = "z (m)",
    label = "t = 0d",
    title = "Temperature Equilibration from t = 0d to t = 3d",
)
plot!(T[4], z, xlabel = "T (K)", ylabel = "z (m)", label = "t = 1.5d")
plot!(T[end], z, xlabel = "T (K)", ylabel = "z (m)", label = "t = 3d")
savefig("eq_temperature_plot.png");
# ![](eq_temperature_plot.png)

# # Analytic Expectations

# We can determine a priori what we expect the final temperature to be in
# equilibrium.

# Regardless of the final water profile in equilibrium, we know that
# the final temperature `T_f` will be a constant across the domain. All
# water that began with a temperature above this point will cool to `T_f`,
# and water that began with a temperature below this point will warm to
# `T_f`. The initial function `T(z)` is equal to `T_f` at a value of
# `z = z̃`. This is the location in space which divides these two groups
# (water that warms over time and water that cools over time) spatially.
# We can solve for `z̃(T_f)` using `T_f = T(z̃)`.

# Next, we can determine the change in energy required to cool
# the water above `z̃` to `T_f`: it is the integral from `z̃` to the surface
# at `z = 0` of ` c θ(z) T(z) `, where `c` is the volumetric heat capacity -
# a constant here - and `θ(z)` is the initial water profile. Compute the energy
# required to warm the water below `z̃` to `T_f` in a similar way, set equal, and solve
# for `T_f`. This results in `T_f = 288.056`, which is very close to the mean `T` we observe
# after 3 days, of `288.054`.

# One could also solve the equation for `ϑ_l` specified by
# ``∂ h/∂ z = 0`` to determine the functional form of the
# equilibrium profile of the liquid water.

# # References
# - Bonan, G.  Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.
# - Balland and Arp, J. Environ. Eng. Sci. 4: 549–558 (2005)

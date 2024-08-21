# This sets up the simulation that mimicks the coarse sand
# lab experiment
# presented in Figures 7 and 8a of
# Lehmann, Assouline, Or  (Phys Rev E 77, 2008).

using CairoMakie
import SciMLBase
import ClimaTimeSteppers as CTS
using Thermodynamics

using ClimaCore
import ClimaParams as CP
using SurfaceFluxes
using StaticArrays
using Dates
using DelimitedFiles: readdlm

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
import SurfaceFluxes.Parameters as SFP

FT = Float64;
earth_param_set = LP.LandParameters(FT)
thermo_params = LP.thermodynamic_parameters(earth_param_set);

# We model evaporation using Monin-Obukhov surface theory.
# In our soil model,
# it is not possible to set the initial condition corresponding to
#  MOST fluxes, but not include radiative fluxes.
# This is because for land surface models does not make sense
# to include atmospheric forcing but not radiative forcing.

# Because of this, we need to supply downward welling short and long wave
# radiation. We chose SW = 0 and LW = σT^4, in order to approximately balance
# out the blackbody emission of the soil which is accounted for by our model.
# Our assumption is that in the lab experiment there was no radiative heating
# or cooling of the soil.

ref_time = DateTime(2005) # required argument, but not used in this case
SW_d = (t) -> 0
LW_d = (t) -> 301.15^4 * 5.67e-8
radiation = PrescribedRadiativeFluxes(
    FT,
    TimeVaryingInput(SW_d),
    TimeVaryingInput(LW_d),
    ref_time,
);
# Set up atmospheric conditions that result in the potential evaporation
# rate obsereved in the experiment.
# Some of these conditions are reported in the paper.
T_air = FT(301.15)
rh = FT(0.38)
esat = Thermodynamics.saturation_vapor_pressure(
    thermo_params,
    T_air,
    Thermodynamics.Liquid(),
)
e = rh * esat
q = FT(0.622 * e / (101325 - 0.378 * e))
precip = (t) -> 0.0
T_atmos = (t) -> T_air
u_atmos = (t) -> 1.0
q_atmos = (t) -> q
h_atmos = FT(0.1)
P_atmos = (t) -> 101325
gustiness = FT(1e-2)
atmos = PrescribedAtmosphere(
    TimeVaryingInput(precip),
    TimeVaryingInput(precip),
    TimeVaryingInput(T_atmos),
    TimeVaryingInput(u_atmos),
    TimeVaryingInput(q_atmos),
    TimeVaryingInput(P_atmos),
    ref_time,
    h_atmos,
    earth_param_set;
    gustiness = gustiness,
);

# Define the boundary conditions
top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(atmos, radiation)
zero_water_flux = WaterFluxBC((p, t) -> 0)
zero_heat_flux = HeatFluxBC((p, t) -> 0)
boundary_fluxes = (;
    top = top_bc,
    bottom = WaterHeatBC(; water = zero_water_flux, heat = zero_heat_flux),
);

# Define the parameters
# n and alpha estimated by matching vG curve.
K_sat = FT(225.1 / 3600 / 24 / 1000)
vg_n = FT(10.0)
vg_α = FT(6.0)
hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
ν = FT(0.43)
θ_r = FT(0.045)
S_s = FT(1e-3)
ν_ss_om = FT(0.0)
ν_ss_quartz = FT(1.0)
ν_ss_gravel = FT(0.0)
emissivity = FT(1.0)
PAR_albedo = FT(0.2)
NIR_albedo = FT(0.4)
z_0m = FT(1e-3)
z_0b = FT(1e-4)
d_ds = FT(0.01)
params = ClimaLand.Soil.EnergyHydrologyParameters(
    FT;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = hcm,
    K_sat,
    S_s,
    θ_r,
    PAR_albedo,
    NIR_albedo,
    emissivity,
    z_0m,
    z_0b,
    earth_param_set,
    d_ds,
);

# Domain - single column
zmax = FT(0)
zmin = FT(-0.35)
nelems = 5
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z;

# Soil model, and create the prognostic vector Y and cache p:
soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = (),
)

Y, p, cds = initialize(soil);

# Set initial conditions
function hydrostatic_equilibrium(z, z_interface, params)
    (; ν, S_s, hydrology_cm) = params
    (; α, n, m) = hydrology_cm
    if z < z_interface
        return -S_s * (z - z_interface) + ν
    else
        return ν * (1 + (α * (z - z_interface))^n)^(-m)
    end
end
function init_soil!(Y, z, params)
    FT = eltype(Y.soil.ϑ_l)
    Y.soil.ϑ_l .= hydrostatic_equilibrium.(z, FT(-0.001), params)
    Y.soil.θ_i .= 0
    T = FT(296.15)
    ρc_s = @. Soil.volumetric_heat_capacity(
        Y.soil.ϑ_l,
        FT(0),
        params.ρc_ds,
        params.earth_param_set,
    )
    Y.soil.ρe_int =
        Soil.volumetric_internal_energy.(FT(0), ρc_s, T, params.earth_param_set)
end
init_soil!(Y, z, soil.parameters);

# Timestepping:
t0 = Float64(0)
tf = Float64(24 * 3600 * 13)
dt = Float64(900.0)

# We also set the initial conditions of the cache here:
set_initial_cache! = make_set_initial_cache(soil)
set_initial_cache!(p, Y, t0);

# Define the tendency functions
exp_tendency! = make_exp_tendency(soil)
imp_tendency! = make_imp_tendency(soil);
jacobian! = ClimaLand.make_jacobian(soil);
jac_kwargs = (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!);

timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

# Define the problem and callbacks:
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
saveat = Array(t0:3600.0:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
updateat = deepcopy(saveat)
model_drivers = ClimaLand.get_drivers(soil)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb);

# Solve
sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat);

# Figures

# Extract the evaporation at each saved step
evap = [
    parent(sv.saveval[k].soil.turbulent_fluxes.vapor_flux_liq)[1] for
    k in 1:length(sol.t)
]
savepath = joinpath(pkgdir(ClimaLand), "docs/tutorials/standalone/Soil/")
evaporation_data =
    ClimaLand.Artifacts.lehmann_assouline_or2008_evaporation_data();
ref_soln_E = readdlm(evaporation_data, ',')
ref_soln_E_350mm = ref_soln_E[2:end, 1:2]
data_dates = ref_soln_E_350mm[:, 1]
data_e = ref_soln_E_350mm[:, 2];

fig = Figure(size = (800, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Day",
    ylabel = "Evaporation rate (mm/d)",
    title = "Bare soil evaporation",
)
CairoMakie.xlims!(minimum(data_dates), maximum(data_dates))
CairoMakie.lines!(
    ax,
    FT.(data_dates),
    FT.(data_e),
    label = "Data",
    color = :blue,
)
CairoMakie.lines!(
    ax,
    sol.t ./ 3600 ./ 24,
    evap .* (1000 * 3600 * 24),
    label = "Model",
    color = :black,
)
CairoMakie.axislegend(ax)

ax = Axis(
    fig[1, 2],
    xlabel = "Mass (g)",
    yticksvisible = false,
    yticklabelsvisible = false,
)
A_col = π * (0.027)^2
mass_0 = sum(sol.u[1].soil.ϑ_l) * 1e6 * A_col
mass_loss =
    [mass_0 - sum(sol.u[k].soil.ϑ_l) * 1e6 * A_col for k in 1:length(sol.t)]
CairoMakie.lines!(
    ax,
    cumsum(FT.(data_e)) ./ (1000 * 24) .* A_col .* 1e6,
    FT.(data_e),
    label = "Data",
    color = :blue,
)
CairoMakie.lines!(
    ax,
    mass_loss,
    evap .* (1000 * 3600 * 24),
    label = "Model",
    color = :black,
)
save("evaporation_lehmann2008_fig8b.png", fig);
# ![](evaporation_lehmann2008_fig8b.png)

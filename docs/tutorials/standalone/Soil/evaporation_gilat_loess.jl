using CairoMakie
import SciMLBase
import ClimaTimeSteppers as CTS
using Thermodynamics

using ClimaCore
import ClimaParams as CP
using SurfaceFluxes
using StaticArrays
using Dates
using ArtifactWrappers
using DelimitedFiles: readdlm

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
import SurfaceFluxes.Parameters as SFP

# Define simulation times
t0 = Float64(0)
tf = Float64(24 * 3600 * 15)
FT = Float64

# Resolution
Δz = 0.01
earth_param_set = LP.LandParameters(FT)
thermo_params = LP.thermodynamic_parameters(earth_param_set)
# Coarse sand experiment described in Figures 7 and 8a
# of Lehmann, Assouline, Or  (Phys Rev E 77, 2008)
K_sat = FT(0.01 / 3600 / 24)
vg_n = FT(1.45)
vg_α = FT(1.5)
hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
ν = FT(0.4)
θ_r = FT(0.04)
S_s = FT(1e-3)
ν_ss_om = FT(0.0)
ν_ss_quartz = FT(0.3)
ν_ss_gravel = FT(0.0)
emissivity = FT(1.0)
PAR_albedo = FT(0.2)
NIR_albedo = FT(0.4)
z_0m = FT(1e-3)
z_0b = FT(1e-4)
d_ds = FT(0.01)# 10mm

ref_time = DateTime(2005)
SW_d = (t) -> 0
LW_d = (t) -> 301.15^4 * 5.67e-8
radiation = PrescribedRadiativeFluxes(
    FT,
    TimeVaryingInput(SW_d),
    TimeVaryingInput(LW_d),
    ref_time,
)
# Atmos
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
    h_atmos;
    gustiness = gustiness,
)
top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(atmos, radiation)
zero_water_flux = WaterFluxBC((p, t) -> 0)
zero_heat_flux = HeatFluxBC((p, t) -> 0)
boundary_fluxes = (;
    top = WaterHeatBC(; water = zero_water_flux, heat = zero_heat_flux),
    bottom = WaterHeatBC(; water = zero_water_flux, heat = zero_heat_flux),
)
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
)

zmax = FT(0)
zmin = FT(-1.6)
nelems = Int((zmax - zmin) / Δz)
dt = Float64(10.0)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z

soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = (),
)

Y, p, cds = initialize(soil) # begins saturated
function estimated_ic(z)
    0.34 / (1 + exp(-(z + 0.165) / 0.005)) + 0.05
end
function init_soil!(Y, z, params)
    FT = eltype(Y.soil.ϑ_l)
    Y.soil.ϑ_l .= estimated_ic.(z)
    Y.soil.θ_i .= 0
    T = FT(294.15)
    ρc_s = @. Soil.volumetric_heat_capacity(Y.soil.ϑ_l, FT(0), params)
    Y.soil.ρe_int =
        Soil.volumetric_internal_energy.(FT(0), ρc_s, T, Ref(params))
end

init_soil!(Y, z, soil.parameters)

# We also set the initial conditions of the auxiliary state here:
set_initial_cache! = make_set_initial_cache(soil)
set_initial_cache!(p, Y, t0)

# Timestepping:
soil_exp_tendency! = make_exp_tendency(soil)
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = soil_exp_tendency!, dss! = ClimaLand.dss!),
    Y,
    (t0, tf),
    p,
)
saveat = Array(t0:3600.0:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
cb = SciMLBase.CallbackSet(saving_cb)
sol_no_evap =
    SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)

# Repeat with evaporation and drainage
boundary_fluxes_evap = (;
    top = top_bc,
    bottom = WaterHeatBC(; water = zero_water_flux, heat = zero_heat_flux),
)
soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes_evap,
    sources = (),
)

Y, p, cds = initialize(soil) # begins saturated
init_soil!(Y, z, soil.parameters)

# We also set the initial conditions of the auxiliary state here:
set_initial_cache! = make_set_initial_cache(soil)
set_initial_cache!(p, Y, t0)

# Timestepping:
soil_exp_tendency! = make_exp_tendency(soil)
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = soil_exp_tendency!, dss! = ClimaLand.dss!),
    Y,
    (t0, tf),
    p,
)
saveat = Array(t0:3600.0:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
updateat = deepcopy(saveat)
updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)
evap = [
    parent(sv.saveval[k].soil.turbulent_fluxes.vapor_flux)[1] for
    k in 1:length(sol.t)
]

## Repeat with no drainage (?), but with evap, in shorter domain
zmax = FT(0)
zmin = FT(-0.16)
nelems = Int((zmax - zmin) / Δz)
dt = Float64(10.0)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
z_no_evap = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z
bottom_water_bc = MoistureStateBC((p, t) -> 0.35)
boundary_fluxes_wet = (;
    top = top_bc,
    bottom = WaterHeatBC(; water = bottom_water_bc, heat = zero_heat_flux),
)

soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes_wet,
    sources = (),
)

Y, p, cds = initialize(soil) # begins saturated
init_soil!(Y, z_no_evap, soil.parameters)

# We also set the initial conditions of the auxiliary state here:
set_initial_cache! = make_set_initial_cache(soil)
set_initial_cache!(p, Y, t0)

# Timestepping:
soil_exp_tendency! = make_exp_tendency(soil)
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = soil_exp_tendency!, dss! = ClimaLand.dss!),
    Y,
    (t0, tf),
    p,
)
saveat = Array(t0:3600.0:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
updateat = deepcopy(saveat)
updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
sol_no_drainage =
    SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)
evap_no_drainage = [
    parent(sv.saveval[k].soil.turbulent_fluxes.vapor_flux)[1] for
    k in 1:length(sol_no_drainage.t)
]



# Figures
savepath = joinpath(pkgdir(ClimaLand), "docs/tutorials/standalone/Soil/")

fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "Day", ylabel = "Evaporation rate (mm/d)")
CairoMakie.lines!(
    ax,
    sol.t ./ 3600 ./ 24,
    evap .* (1000 * 3600 * 24),
    label = "With drainage",
    color = :red,
)
CairoMakie.lines!(
    ax,
    sol_no_drainage.t ./ 3600 ./ 24,
    evap_no_drainage .* (1000 * 3600 * 24),
    label = "No drainage",
    color = :blue,
)

CairoMakie.axislegend(ax)
ax2 = Axis(fig[1, 2], xlabel = "Day", ylabel = "Cumulative evaporation (mm)")
CairoMakie.lines!(
    ax2,
    sol.t ./ 3600 ./ 24,
    cumsum(evap) .* (1000 * 3600),
    color = :red,
)
CairoMakie.lines!(
    ax2,
    sol_no_drainage.t ./ 3600 ./ 24,
    cumsum(evap_no_drainage) .* (1000 * 3600),
    color = :blue,
)
save(joinpath(savepath, "evaporation_lehmann2024_figS6_$Δz.png"), fig)

fig2 = Figure(size = (800, 1200))
ax1 = Axis(
    fig2[1, 1],
    xlabel = "Volumetric Water Content",
    ylabel = "Depth(cm)",
    title = "Drainage only",
)
CairoMakie.ylims!(-0.35, 0)
CairoMakie.xlims!(0.0, 0.4)
linestyles = [:solid, :dash, :dashdot, :dashdotdot, :dot]
hours = [0, 1, 2, 10]
for i in 1:1:4
    CairoMakie.lines!(
        ax1,
        parent(sol_no_evap.u[hours[i] * 24 + 1].soil.ϑ_l)[:],
        parent(z)[:],
        label = "$(hours[i])",
        color = :black,
        linestyle = linestyles[i],
    )
end
ax2 = Axis(fig2[2, 1], title = "Evap+Drainage")
CairoMakie.ylims!(-0.3, 0)
CairoMakie.xlims!(0.0, 0.4)
hours = [0, 1, 2, 5, 13]
for i in 1:1:5
    CairoMakie.lines!(
        ax2,
        parent(sol.u[hours[i] * 24 + 1].soil.ϑ_l)[:],
        parent(z)[:],
        label = "$(hours[i])",
        color = :black,
        linestyle = linestyles[i],
    )
end
ax3 = Axis(fig2[3, 1], title = "Evap only")
CairoMakie.ylims!(-0.15, 0)
CairoMakie.xlims!(0.0, 0.4)
hours = [0, 2, 9, 14]
for i in 1:1:4
    CairoMakie.lines!(
        ax3,
        parent(sol_no_drainage.u[hours[i] * 24 + 1].soil.ϑ_l)[:],
        label = "$(hours[i])",
        parent(z_no_evap)[:],
        color = :black,
        linestyle = linestyles[i],
    )
end

CairoMakie.axislegend(ax3, position = :rb)
CairoMakie.axislegend(ax2, position = :rb)
CairoMakie.axislegend(ax1, position = :rb)
save(joinpath(savepath, "evaporation_gardner_fig1_$Δz.png"), fig2)

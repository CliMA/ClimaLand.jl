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

earth_param_set = LP.LandParameters(FT)
thermo_params = LP.thermodynamic_parameters(earth_param_set)
# Coarse sand experiment described in Figures 7 and 8a
# of Lehmann, Assouline, Or  (Phys Rev E 77, 2008)
K_sat = FT(0.01 / 3600 / 24)
# n and alpha estimated by matching vG curve.
vg_n = FT(1.55)
vg_α = FT(1.5)
hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
ν = FT(0.42)
θ_r = FT(0.075)
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
#top_water_flux = WaterFluxBC((p, t) -> 4.6296296296296295e-8)
boundary_fluxes = (;
    top = top_bc,#WaterHeatBC(; water = zero_water_flux, heat = zero_heat_flux),#top_bc,
    bottom = WaterHeatBC(; water = zero_water_flux, heat = zero_heat_flux),
)
params = ClimaLand.Soil.EnergyHydrologyParameters{FT}(;
    ν = ν,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ν_ss_gravel = ν_ss_gravel,
    hydrology_cm = hcm,
    K_sat = K_sat,
    S_s = S_s,
    θ_r = θ_r,
    PAR_albedo = PAR_albedo,
    NIR_albedo = NIR_albedo,
    emissivity = emissivity,
    z_0m = z_0m,
    z_0b = z_0b,
    earth_param_set = earth_param_set,
    d_ds = d_ds,
)

zmax = FT(0)
zmin = FT(-1.6)
nelems = 64
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
   0.3/(1+exp(-(z+0.15)/0.01))+0.09
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
updateat = deepcopy(saveat)
updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)


# Figures

evap = [
    parent(sv.saveval[k].soil.turbulent_fluxes.vapor_flux)[1] for
    k in 1:length(sol.t)
]

savepath = joinpath(pkgdir(ClimaLand), "docs/tutorials/standalone/Soil/")

fig = Figure(size = (1200, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Day",
    ylabel = "Evaporation rate (mm/d)",
)
CairoMakie.lines!(
    ax,
    sol.t ./ 3600 ./ 24,
    evap .* (1000 * 3600 * 24),
    label = "Model",
    color = :black,
)
CairoMakie.axislegend(ax)
ax2 = Axis(
    fig[1, 2],
    xlabel = "Day",
    ylabel = "Cumulative evaporation (mm)",
)
CairoMakie.lines!(
    ax2,
    sol.t ./ 3600 ./ 24,
    cumsum(evap) .* (1000 * 3600),
    label = "Model",
    color = :black,
)

ax3 = Axis(
    fig[1, 3],
    xlabel = "Volumetric Water Content",
    ylabel = "Depth(cm)",
)
CairoMakie.ylims!(-0.3,0)
CairoMakie.xlims!(0.0,0.4)
for i in [0, 1,2,3,6,13].*24
    CairoMakie.lines!(
        ax3,
        parent(sol.u[i+1].soil.ϑ_l)[:],
        parent(z)[:],
        label = "$(div(i,24))",
        color = :black,
    )
end

CairoMakie.axislegend(ax3, position = :rb)




save(joinpath(savepath, "evaporation_lehmann2024_figS6_sc1.png"), fig)

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
tf = Float64(24 * 3600 * 13)
dt = Float64(2)
FT = Float64

earth_param_set = LP.LandParameters(FT)
thermo_params = LP.thermodynamic_parameters(earth_param_set)
# Coarse sand experiment described in Figures 7 and 8a
# of Lehmann, Assouline, Or  (Phys Rev E 77, 2008)
K_sat = FT(225.1 / 3600 / 24 / 1000)
# n and alpha estimated by matching vG curve.
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
                   top = top_bc,
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

#TODO: Run with higher resolution once we have the implicit stepper
zmax = FT(0)
zmin = FT(-0.35)
nelems = 5
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z

soil = Soil.EnergyHydrology{FT}(;
                                parameters = params,
                                domain = soil_domain,
                                boundary_conditions = boundary_fluxes,
                                sources = (),
                                )

Y, p, cds = initialize(soil) # begins saturated
function init_soil!(Y, z, params)
    ν = params.ν
    FT = eltype(ν)
    Y.soil.ϑ_l .= ν - 1e-2
    Y.soil.θ_i .= 0
    T = FT(301.15)
    ρc_s = Soil.volumetric_heat_capacity(ν, FT(0), params)
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
    CTS.ClimaODEFunction(
        T_exp! = soil_exp_tendency!,
        dss! = ClimaLand.dss!,
    ),
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
sol =
    SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)

evap = [
    parent(sv.saveval[k].soil.turbulent_fluxes.vapor_flux)[1] for
    k in 1:length(sol.t)
]
# Read in reference solution from artifact
evap_dataset = ArtifactWrapper(
    @__DIR__,
    "lehmann2008_fig8_evaporation",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/cgppw3tx6zdz7h02yt28ri44g1j088ju.csv",
        filename = "lehmann2008_fig8_evaporation.csv",
    ),],
)
evap_datapath = get_data_folder(evap_dataset)
ref_soln_E = readdlm(
    joinpath(evap_datapath, "lehmann2008_fig8_evaporation.csv"),
    ',',
)
ref_soln_E_350mm = ref_soln_E[2:end, 1:2]
data_dates = ref_soln_E_350mm[:, 1]
data_e = ref_soln_E_350mm[:, 2]

# Compare our data to Figure 8b of Lehmann, Assouline, Or  (Phys Rev E 77, 2008)
fig = Figure(size = (400, 400))
ax = Axis(fig[1,1],  xlabel = "Day", ylabel = "Evaporation rate (mm/d)",  title = "Bare soil evaporation")
xlims!(minimum(data_dates), maximum(data_dates))
lines!(ax,
       FT.(data_dates),
       FT.(data_e),
       label = "Data",
       color = :blue
       )
lines!(ax,
       sol.t ./ 3600 ./ 24,
       evap .* (1000 * 3600 * 24),
       label = "Model",
       color = :black,
       )
axislegend(ax)
save(joinpath(savepath, "evaporation_lehmann2008_fig8b.png"), fig)

import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
using ClimaLSM.Soil.Biogeochemistry
using ClimaLSM.Soil.Biogeochemistry: MicrobeProduction
using Dates

import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64
earth_param_set = create_lsm_parameters(FT)

# Make soil model args
ν = FT(0.556)
K_sat = FT(0.0443 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(2.0)
vg_α = FT(2.6) # inverse meters
hcm = vanGenuchten(; α = vg_α, n = vg_n)
θ_r = FT(0.1)
ν_ss_om = FT(0.0)
ν_ss_quartz = FT(1.0)
ν_ss_gravel = FT(0.0)
κ_minerals = FT(2.5)
κ_om = FT(0.25)
κ_quartz = FT(8.0)
κ_air = FT(0.025)
κ_ice = FT(2.21)
κ_liq = FT(0.57)
ρp = FT(2.66 / 1e3 * 1e6)
ρc_ds = FT(2e6 * (1.0 - ν))
κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
κ_dry = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)

soil_ps = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = κ_dry,
    κ_sat_frozen = κ_sat_frozen,
    κ_sat_unfrozen = κ_sat_unfrozen,
    ρc_ds = ρc_ds,
    ν = ν,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ν_ss_gravel = ν_ss_gravel,
    hydrology_cm = hcm,
    K_sat = K_sat,
    S_s = S_s,
    θ_r = θ_r,
    earth_param_set = earth_param_set,
)

zmax = FT(0)
zmin = FT(-1)
nelems = 20
Δz = abs(zmax - zmin) / nelems

lsm_domain = Column(; zlim = (zmin, zmax), nelements = nelems)

top_flux_bc_w = Soil.FluxBC((p, t) -> eltype(t)(-0.00001))
bot_flux_bc_w = Soil.FreeDrainage()

top_flux_bc_h = Soil.FluxBC((p, t) -> eltype(t)(0.0))
bot_flux_bc_h = Soil.FluxBC((p, t) -> eltype(t)(0.0))


sources = (PhaseChange{FT}(Δz),)
boundary_fluxes = (;
    top = (water = top_flux_bc_w, heat = bot_flux_bc_h),
    bottom = (water = bot_flux_bc_w, heat = bot_flux_bc_h),
)
soil_args = (;
    boundary_conditions = boundary_fluxes,
    sources = sources,
    domain = lsm_domain,
    parameters = soil_ps,
)

# Make biogeochemistry model args
Csom = (z, t) -> eltype(t)(5.0)

co2_parameters = Soil.Biogeochemistry.SoilCO2ModelParameters{FT}(;
    earth_param_set = earth_param_set,
)
C = FT(100)

co2_top_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> eltype(t)(0.0))
co2_bot_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> eltype(t)(0.0))
co2_sources = (MicrobeProduction{FT}(),)
co2_boundary_conditions =
    (; top = (CO2 = co2_top_bc,), bottom = (CO2 = co2_bot_bc,))

# Make a PrescribedAtmosphere - we only care about atmos_p though
precipitation_function = (t) -> 1.0
snow_precip = (t) -> 1.0
atmos_T = (t) -> 1.0
atmos_u = (t) -> 1.0
atmos_q = (t) -> 1.0
atmos_p = (t) -> 100000.0
UTC_DATETIME = Dates.now()
atmos_h = FT(30)
atmos_co2 = (t) -> 1.0

atmos = ClimaLSM.PrescribedAtmosphere(
    precipitation_function,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    UTC_DATETIME,
    atmos_h;
    c_co2 = atmos_co2,
)

soil_drivers = Soil.Biogeochemistry.SoilDrivers(
    Soil.Biogeochemistry.PrognosticMet(),
    Soil.Biogeochemistry.PrescribedSOC(Csom),
    atmos,
)

soilco2_args = (;
    boundary_conditions = co2_boundary_conditions,
    sources = co2_sources,
    domain = lsm_domain,
    parameters = co2_parameters,
    drivers = soil_drivers,
)

# Create integrated model instance
model = LandSoilBiogeochemistry{FT}(;
    soil_args = soil_args,
    soilco2_args = soilco2_args,
)

Y, p, coords = initialize(model)
set_initial_aux_state! = make_set_initial_aux_state(model);

function init_soil!(Y, z, params)
    ν = params.ν
    FT = eltype(Y.soil.ϑ_l)
    Y.soil.ϑ_l .= FT(0.33)
    Y.soil.θ_i .= FT(0.1)
    T = FT(279.85)
    ρc_s = Soil.volumetric_heat_capacity(FT(0.33), FT(0.1), params)
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(FT(0.0), ρc_s, T, Ref(params))
end

function init_co2!(Y, z)
    function CO2_profile(z::FT) where {FT}
        C = FT(0.0)
        return FT(C)
    end
    Y.soilco2.C .= CO2_profile.(z)
end

z = coords.subsurface.z
init_soil!(Y, z, model.soil.parameters)
init_co2!(Y, z)
t0 = FT(0.0)
set_initial_aux_state!(p, Y, t0);
Soil_bio_exp_tendency! = make_exp_tendency(model)

tf = FT(10000)
dt = FT(10)

timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

saveat = collect(t0:FT(10 * dt):tf)
saved_values = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
cb = ClimaLSM.NonInterpSavingCallback(saved_values, saveat)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = Soil_bio_exp_tendency!),
    Y,
    (t0, tf),
    p,
)
sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb)

# Animation
# You will need to ]add GLMakie to your base Julia Project.toml
#=
using GLMakie
fig = Figure(resolution = (1100, 700))

depth = parent(coords.subsurface.z)[:]
t = Observable(1)

Mcolor = @lift(parent(sol.u[$t].soil.ϑ_l)[:])
Tcolor = @lift(parent(saved_values.saveval[$t].soil.T)[:])
Ccolor = @lift(parent(sol.u[$t].soilco2.C)[:])
Icolor = @lift(parent(sol.u[$t].soil.θ_i)[:])

figtitle = @lift("Time: 10s * " * string($t))

ax_M = Axis(
    fig[1, 1],
    title = "Moisture",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    xtrimspine = true,
    ytrimspine = true,
    ylabel = "depth (m)",
    xticksvisible = false,
    xticklabelsvisible = false,
)

ax_T = Axis(
    fig[1, 3],
    title = "Temperature",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    xtrimspine = true,
    ytrimspine = true,
    yticksvisible = false,
    yticklabelsvisible = false,
    xticksvisible = false,
    xticklabelsvisible = false,
)

ax_C = Axis(
    fig[1, 5],
    title = "CO2",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    xtrimspine = true,
    ytrimspine = true,
    yticksvisible = false,
    yticklabelsvisible = false,
    xticksvisible = false,
    xticklabelsvisible = false,
)

ax_I = Axis(
    fig[1, 7],
    title = "Ice",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    xtrimspine = true,
    ytrimspine = true,
    yticksvisible = false,
    yticklabelsvisible = false,
    xticksvisible = false,
    xticklabelsvisible = false,
)

lines!(
    ax_M,
    ones(nelems),
    depth,
    color = Mcolor,
    linewidth = 140,
    colormap = :deep,
    colorrange = [0.31, 0.6],
)
lines!(
    ax_T,
    ones(nelems),
    depth,
    color = Tcolor,
    linewidth = 140,
    colormap = Reverse(:autumn1),
    colorrange = [274, 280],
)
lines!(
    ax_C,
    ones(nelems),
    depth,
    color = Ccolor,
    linewidth = 140,
    colormap = Reverse(:grays),
    colorrange = [0.0, 0.006],
)
lines!(
    ax_I,
    ones(nelems),
    depth,
    color = Icolor,
    linewidth = 140,
    colormap = :Blues,
    colorrange = [0.0, 0.1],
)

Colorbar(
    fig[1, 2],
    limits = (0.31, 0.35),
    label = "Soil Moisture",
    colormap = :deep,
)
Colorbar(
    fig[1, 4],
    limits = (274, 280),
    label = "Soil Temperature",
    colormap = Reverse(:autumn1),
)
Colorbar(
    fig[1, 6],
    limits = (0.0, 0.006),
    label = "Soil CO2",
    colormap = Reverse(:grays),
)
Colorbar(fig[1, 8], limits = (0.0, 0.1), label = "Ice", colormap = :Blues)

xlims!(ax_M, (0, 2))
xlims!(ax_T, (0, 2))
xlims!(ax_C, (0, 2))
xlims!(ax_I, (0, 2))

hidespines!(ax_M, :b)
hidespines!(ax_T, :b, :l)
hidespines!(ax_C, :b, :l)
hidespines!(ax_I, :b, :l)

Label(fig[0, :], text = figtitle, fontsize = 30)

framerate = 100
timestamps = range(1, 1000, step = 1)

filename = "./soil_biogeochem_animation.mp4"
record(fig, filename, timestamps; framerate = framerate) do dayt
    t[] = dayt
end
=#

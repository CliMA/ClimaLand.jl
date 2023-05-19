using Plots
using OrdinaryDiffEq: ODEProblem, solve, RK4, Euler, step!, init
using ClimaCore
import CLIMAParameters as CP
using Thermodynamics
using Insolation

using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil

import ClimaLSM
import ClimaLSM.Parameters as LSMP
import SurfaceFluxes.Parameters as SFP
using RootSolvers
using SurfaceFluxes
using StaticArrays

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
FT = Float64
earth_param_set = create_lsm_parameters(FT)
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
# Coarse sand experiment described in B of Lehmann (2008)
K_sat = FT(225.1 / 3600 / 24 / 1000)
# n and alpha estimated by matching vG curve.
vg_n = FT(10.0)
vg_α = FT(6.0)
ν = FT(0.43)
θ_r = FT(0.045)
S_s = FT(1e-3)
vg_m = FT(1) - FT(1) / vg_n
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
κ_dry_soil = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)
emissivity = FT(1.0)
albedo = FT(0.2)
z_0m = 1e-4
z_0b = 1e-5

SW_d = (t) -> eltype(t)(0)
LW_d = (t) -> eltype(t)(301.15^4 * 5.67e-8)
radiation = PrescribedRadiativeFluxes(
    FT,
    SW_d,
    LW_d;
    orbital_data = Insolation.OrbitalData(
        joinpath(pkgdir(ClimaLSM), "artifacts"),
    ),
)
# Atmos
T_air = 301.15
rh = 0.38
esat = Thermodynamics.saturation_vapor_pressure(
    thermo_params,
    T_air,
    Thermodynamics.Liquid(),
)
e = rh * esat
q = FT(0.622 * e / (101325 - 0.378 * e))
precip = (t) -> eltype(t)(0.0)
T_atmos = (t) -> eltype(t)(T_air)
u_atmos = (t) -> eltype(t)(0.1)
q_atmos = (t) -> eltype(t)(q)
h_atmos = FT(0.1)
P_atmos = (t) -> eltype(t)(101325)
atmos = PrescribedAtmosphere(
    precip,
    precip,
    T_atmos,
    u_atmos,
    q_atmos,
    P_atmos,
    h_atmos,
)
top_bc = ClimaLSM.Soil.AtmosDrivenFluxBC(atmos, radiation)
zero_flux = FluxBC((p, t) -> eltype(t)(0.0))
boundary_fluxes =
    (; top = top_bc, bottom = (water = zero_flux, heat = zero_flux))
params = ClimaLSM.Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = κ_dry_soil,
    κ_sat_frozen = κ_sat_frozen,
    κ_sat_unfrozen = κ_sat_unfrozen,
    ρc_ds = ρc_ds,
    ν = ν,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ν_ss_gravel = ν_ss_gravel,
    vg_α = vg_α,
    vg_n = vg_n,
    K_sat = K_sat,
    S_s = S_s,
    θ_r = θ_r,
    albedo = albedo,
    emissivity = emissivity,
    z_0m = z_0m,
    z_0b = z_0b,
    earth_param_set = earth_param_set,
)

#TODO: Run with higher resolution once we have the implicit stepper
zmax = FT(0)
zmin = FT(-1.0)
nelems = 10
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
z = ClimaCore.Fields.coordinate_field(soil_domain.space).z

soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = (),
)

Y, p, cds = initialize(soil) # begins saturated
function init_soil!(Y, z, params)
    ν = params.ν
    FT = typeof(ν)
    Y.soil.ϑ_l .= ν - 1e-2
    Y.soil.θ_i .= 0
    T = FT(301.15)
    ρc_s = Soil.volumetric_heat_capacity(ν, FT(0), params)
    Y.soil.ρe_int =
        Soil.volumetric_internal_energy.(FT(0), ρc_s, T, Ref(params))
end

t = FT(0)
t0 = FT(0)
tf = FT(24 * 3600 * 15)
dt = FT(100)

init_soil!(Y, cds.z, soil.parameters)
soil_exp_tendency! = make_exp_tendency(soil)

prob = ODEProblem(soil_exp_tendency!, Y, (t0, tf), p)
sol = solve(prob, RK4(); dt = dt, saveat = 3600);

(; ν, vg_m, vg_n, θ_r, d_ds) = soil.parameters
_D_vapor = FT(LSMP.D_vapor(soil.parameters.earth_param_set))
update_aux! = make_update_aux(soil)
S_c::FT = (1 + ((vg_n - 1) / vg_n)^(1 - 2 * vg_n))^(-vg_m)
evap = []
evap_0 = []
r_ae = []
r_soil = []
dsl = []
T_soil = []
q_soil = []
ρ_soil = []
ρ_atmos = []
L_MO = []
surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)

for i in 1:length(sol.t)
    time = sol.t[i]
    u = sol.u[i]
    update_aux!(p, u, time)
    conditions = surface_fluxes(top_bc.atmos, soil, u, p, time)
    τ_a = ClimaLSM.Domains.top_center_to_surface(
        @. max(eps(FT), (ν - p.soil.θ_l - u.soil.θ_i)^(FT(5 / 2)) / ν)
    )
    S_l_sfc = ClimaLSM.Domains.top_center_to_surface(
        effective_saturation.(ν, u.soil.ϑ_l, θ_r),
    )
    layer_thickness =
        parent(Soil.dry_soil_layer_thickness.(S_l_sfc, S_c, d_ds))[1]
    push!(dsl, layer_thickness)
    push!(r_soil, parent(@. layer_thickness / (_D_vapor * τ_a))[1])
    push!(r_ae, parent(@. 1 / (conditions.Ch * abs(top_bc.atmos.u(time))))[1])
    push!(evap, parent(conditions.vapor_flux)[1])
    T_sfc = parent(p.soil.T)[end]
    push!(T_soil, T_sfc)
    ts_in = ClimaLSM.construct_atmos_ts(top_bc.atmos, time, thermo_params)
    push!(ρ_atmos, Thermodynamics.air_density(thermo_params, ts_in))
    ρ_sfc = compute_ρ_sfc(thermo_params, ts_in, T_sfc)
    push!(ρ_soil, ρ_sfc)
    q_sat = Thermodynamics.q_vap_saturation_generic(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Liquid(),
    )
    q_sfc =
        parent(ClimaLSM.surface_specific_humidity(soil, Y, p, T_sfc, ρ_sfc))[1]
    push!(q_soil, q_sfc)
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)
    state_sfc = SurfaceFluxes.SurfaceValues(0.0, SVector{2, FT}(0, 0), ts_sfc)
    state_in = SurfaceFluxes.InteriorValues(
        h_atmos,
        SVector{2, FT}(u_atmos(time), 0),
        ts_in,
    )

    # State containers
    sc = SurfaceFluxes.ValuesOnly{FT}(;
        state_in,
        state_sfc,
        z0m = z_0m,
        z0b = z_0b,
        beta = 1.0,
    )
    potential_conditions = SurfaceFluxes.surface_conditions(
        surface_flux_params,
        sc;
        tol_neutral = SFP.cp_d(surface_flux_params) / 100000,
    )
    vapor_flux =
        SurfaceFluxes.evaporation(
            surface_flux_params,
            sc,
            potential_conditions.Ch,
        ) / 1000.0
    push!(evap_0, vapor_flux)
    push!(L_MO, potential_conditions.L_MO)
end
plt1 = Plots.plot()
Plots.plot!(
    plt1,
    sol.t ./ 3600 ./ 24,
    (evap .* r_ae ./ (r_ae .+ r_soil)) ./ evap_0,
    xlabel = "Days",
    ylabel = "E/E₀",
    label = "",
)
total_moisture_in_mm =
    [sum(sol.u[k].soil.ϑ_l) for k in 1:length(sol.t)] * 1000.0

plt2 = Plots.plot()
Plots.plot!(
    plt2,
    sol.t ./ 3600 ./ 24,
    total_moisture_in_mm,
    xlabel = "Days",
    ylabel = "∫θdz (mm)",
    label = "",
)

plt3 = Plots.plot()
Plots.plot!(
    plt3,
    sol.t ./ 3600 ./ 24,
    r_soil,
    xlabel = "Days",
    ylabel = "R_soil (m/s)",
    label = "",
)

plt4 = Plots.plot()
Plots.plot!(
    plt4,
    sol.t ./ 3600 ./ 24,
    (evap .* r_ae ./ (r_ae .+ r_soil)) .* (1000 * 3600 * 24),
    xlabel = "Days",
    ylabel = "E (mm/d)",
    label = "",
)
plt5 = Plots.plot()
Plots.plot!(
    plt5,
    sol.t ./ 3600 ./ 24,
    T_soil,
    xlabel = "Days",
    ylabel = "T_sfc (K)",
    label = "",
)

plt6 = Plots.plot()
Plots.plot!(
    plt6,
    sol.t ./ 3600 ./ 24,
    q_soil,
    xlabel = "Days",
    ylabel = "q_sfc",
    label = "",
)
top = [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:length(sol.t)]
S = Soil.effective_saturation.(ν, top, θ_r)
ψ = Soil.matric_potential.(vg_α, vg_n, vg_m, S)
plt7 = Plots.plot()
Plots.plot!(
    plt7,
    sol.t ./ 3600 ./ 24,
    top,
    xlabel = "Days",
    ylabel = "Soil Moisture",
    label = "5cm",
)

Plots.plot!(
    plt7,
    sol.t ./ 3600 ./ 24,
    [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:length(sol.t)],
    xlabel = "Days",
    ylabel = "Soil Moisture",
    label = "15cm",
)
Plots.plot!(
    plt7,
    sol.t ./ 3600 ./ 24,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:length(sol.t)],
    xlabel = "Days",
    ylabel = "Soil Moisture",
    label = "25cm",
)

plt8 = Plots.plot()
Plots.plot!(
    plt8,
    sol.t ./ 3600 ./ 24,
    ψ,
    xlabel = "Days",
    ylabel = "ψ_sfc",
    label = "",
)
Plots.plot(plt1, plt2, plt3, plt4; layout = (2, 2))

savepath = joinpath(pkgdir(ClimaLSM), "experiments/Standalone/Soil")
Plots.savefig(joinpath(savepath, "evaporation_from_coarse_sand1.png"))
Plots.plot(plt5, plt6, plt7, plt8; layout = (2, 2))
Plots.savefig(joinpath(savepath, "evaporation_from_coarse_sand2.png"))

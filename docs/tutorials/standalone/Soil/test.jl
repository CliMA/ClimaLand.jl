import SciMLBase
using Statistics
using Plots

using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil

import ClimaLand
import ClimaLand.Parameters as LP

FT = Float32
earth_param_set = LP.LandParameters(FT);

ν = FT(0.53)
ν_ss_quartz = FT(0.36)
ν_ss_om = FT(0.05)
ν_ss_gravel = FT(0.58)
Ksat = FT(4.4e-5) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.44)
vg_α = FT(3.12) # inverse meters
hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
θ_r = FT(0.14)
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

zmax = FT(0)
zmin = FT(-10.0)
nelems = 15
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems, dz_tuple = FT.((2.0, 0.05)));

surface_water_flux = WaterFluxBC((p, t) -> -2e-7)
bottom_water_flux = FreeDrainage();
surface_heat_flux = HeatFluxBC((p, t) -> -24.0)
bottom_heat_flux = HeatFluxBC((p, t) -> 0.0);
boundary_fluxes = (;
    top = WaterHeatBC(; water = surface_water_flux, heat = surface_heat_flux),
    bottom = EnergyFreeDrainage()
);

sources = ()#(PhaseChange{FT}(),);

soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
);
exp_tendency! = make_exp_tendency(soil);
imp_tendency! = make_imp_tendency(soil);
jacobian! = ClimaLand.make_jacobian(soil);

Y, p, coords = initialize(soil);
Y.soil.ϑ_l .= 0.35
Y.soil.θ_i .= 0.0# 0.07
θ_l = Soil.volumetric_liquid_fraction.(Y.soil.ϑ_l, ν, θ_r)
ρc_s =
    Soil.volumetric_heat_capacity.(
        θ_l,
        Y.soil.θ_i,
        params.ρc_ds,
        params.earth_param_set,
    )
T = FT(288)
Y.soil.ρe_int .=
    Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T,
        params.earth_param_set,
        )

t0 = Float64(0)
tf = Float64(60 * 60 * 24*500);

set_initial_cache! = make_set_initial_cache(soil);
set_initial_cache!(p, Y, t0);

dt = Float64(1000.0);
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

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
saveat = collect(t0:FT(30000):tf)
saved_values = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
);
cb = ClimaLand.NonInterpSavingCallback(saved_values, saveat);
sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = saveat, callback = cb);
t = sol.t;
Np = 6
plot([parent(sol.u[k].soil.ϑ_l)[1] for k in 1:length(t)], xlabel = "ϑ", label= "1")
for i in 2:Np
    plot!([parent(sol.u[k].soil.ϑ_l)[i] for k in 1:length(t)], label = "$i")
end
savefig("eq_moisture_plot_bc_bot.png");


#plot([parent(sol.u[k].soil.θ_i)[1] for k in 1:length(t)], xlabel = "θ_i", label= "1")
#for i in 2:15
#    plot!([parent(sol.u[k].soil.θ_i)[i] for k in 1:length(t)], label = "$i")
#end

#savefig("eq_ic_plot.png");
# ![](eq_moisture_plot.png)

plot([parent(saved_values.saveval[k].soil.T)[1] for k in 1:length(t)], xlabel = "T", label = "1")
for i in 2:Np
    plot!([parent(saved_values.saveval[k].soil.T)[i] for k in 1:length(t)], label = "$i")
end
savefig("eq_temperature_plot_bc_bot.png");

plot([parent(saved_values.saveval[k].soil.K)[1] for k in 1:length(t)], xlabel = "K")
for i in 2:Np
    plot!([parent(saved_values.saveval[k].soil.K)[i] for k in 1:length(t)], label = "$i")
end
savefig("eq_K_plot_bc_bot.png");

plot([parent(sol.u[k].soil.ρe_int)[1] for k in 1:length(t)], xlabel = "ρe")
for i in 2:Np
    plot!([parent(sol.u[k].soil.ρe_int)[i] for k in 1:length(t)], label = "$i")
end
savefig("eq_energy_plot_bc_bot.png");
# ![](eq_temperature_plot.png)
#=
# Compute fluxes:
z = parent(soil.domain.fields.z)[:];
function compute_water_flux(Y, p, z)
    top_bc = FT(-2e-7)
    h = parent(p.soil.ψ)[:] .+ z
    ∂h∂z = (h[2:end] .- h[1:end-1]) ./ (z[2:end] .- z[1:end-1])
    K = parent(p.soil.K)[:]
    K_face = (K[2:end].+ K[1:end-1])./2
    F = @. -K_face *  ∂h∂z
    bottom_bc = -K[1]
    F_complete = [bottom_bc, F..., top_bc]
    return F_complete
end

function compute_water_energy_flux(Y, p, z)
    h = parent(p.soil.ψ)[:] .+ z
    ∂h∂z = (h[2:end] .- h[1:end-1]) ./ (z[2:end] .- z[1:end-1])
    K = parent(p.soil.K)[:]
    K_face = (K[2:end].+ K[1:end-1])./2
    T = parent(p.soil.T)[:]
    ρe_int_liq = ClimaLand.Soil.volumetric_internal_energy_liq.(
        T,
        soil.parameters.earth_param_set,
    )
    ρe_int_liq_face = (ρe_int_liq[2:end].+ ρe_int_liq[1:end-1])./2
    F = @. -K_face *  ρe_int_liq_face * ∂h∂z
    bottom_bc = 0
    F_complete = [bottom_bc, F..., top_bc]
    return F_complete
end

function compute_gradT_energy_flux(Y, p, z)
    top_bc = FT(-24)
    T = parent(p.soil.T)[:]
    ∂T∂z = (T[2:end] .- T[1:end-1]) ./ (z[2:end] .- z[1:end-1])
    κ = parent(p.soil.κ)[:]
    κ_face = (κ[2:end].+ κ[1:end-1])./2
    F = @. -κ_face * ∂T∂z
    bottom_bc = 0
    F_complete = [bottom_bc, F..., top_bc]
    return F_complete
end

F_liq = [compute_water_flux(sol.u[k], saved_values.saveval[k],z) for k in 1:length(sol.t)]
F_liq_energy = [compute_water_energy_flux(sol.u[k], saved_values.saveval[k],z) for k in 1:length(sol.t)]
F_gradT = [compute_gradT_energy_flux(sol.u[k], saved_values.saveval[k],z) for k in 1:length(sol.t)]

plot([F_gradT[k][1] for k in 1:length(t)], xlabel = "-κ∂T∂z", label = "1-1/2")
for i in 2:(Np+1)
    plot!([F_gradT[k][i] for k in 1:length(t)], label = "$i-1/2")
end
savefig("eq_gradT_plot_bc_bot.png");

plot([F_liq_energy[k][1] for k in 1:length(t)], xlabel = "-Kρe_int_liq ∂h∂z", label = "1-1/2")
for i in 2:(Np+1)
    plot!([F_liq_energy[k][i] for k in 1:length(t)], label = "$i-1/2")
end
savefig("eq_liq_energy_plot_bc_bot.png");

plot([F_liq[k][1] for k in 1:length(t)], xlabel = "-K∂h∂z", label = "1-1/2")
for i in 2:(Np+1)
    plot!([F_liq[k][i] for k in 1:length(t)], label = "$i-1/2")
end
savefig("eq_liq_flux_plot_bc_bot.png");
=#

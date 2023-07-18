using ClimaCore
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using Plots
using Statistics
using DelimitedFiles
using ArtifactWrappers

using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: Column

FT = Float64
rmse(v1, v2) = sqrt(mean((v1 .- v2) .^ 2))

# Read in reference solutions from artifacts
flux_dataset = ArtifactWrapper(
    @__DIR__,
    "richards_flux_bc",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/bsfokpg0wvxoq04e8na0t3o0u6x5yw9n.csv",
        filename = "ref_soln_flux.csv",
    ),],
)
flux_datapath = get_data_folder(flux_dataset)
ref_soln_flux = readdlm(joinpath(flux_datapath, "ref_soln_flux.csv"))

dirichlet_dataset = ArtifactWrapper(
    @__DIR__,
    "richards_dirichlet_bc",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/w6q30flbgj68lr0ncvoco10okrupmab1.csv",
        filename = "ref_soln_dirichlet.csv",
    ),],
)
dirichlet_datapath = get_data_folder(dirichlet_dataset)
ref_soln_dirichlet =
    readdlm(joinpath(dirichlet_datapath, "ref_soln_dirichlet.csv"))

stepper = CTS.ARS111()
convergence_cond = CTS.MaximumError(FT(1e-8))
conv_checker = CTS.ConvergenceChecker(norm_condition = convergence_cond)
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 50,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

t_start = FT(0)
t_end = FT(1e6)
dts = FT.([10, 100, 1000, 10000])

# van Genuchten parameters for clay (from Bonan 2019 supplemental program 8.2)
ν = FT(0.495)
K_sat = FT(0.0443 / 3600 / 100) # m/s
vg_n = FT(1.43)
vg_α = FT(0.026 * 100) # inverse meters
θ_r = FT(0.124)
S_s = FT(1e-3) #inverse meters
hcm = vanGenuchten(; α = vg_α, n = vg_n)

zmax = FT(0)
zmin = FT(-0.5)
nelems = 50

params = Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
sources = ()

# Set flux boundary conditions (used for calculating mass balance)
flux_in = FT(-1e-7)
top_bc = Soil.FluxBC((p, t) -> eltype(t)(flux_in))
flux_out = FT(0)
bot_bc = Soil.FluxBC((p, t) -> eltype(t)(flux_out))

boundary_fluxes = (; top = (water = top_bc,), bottom = (water = bot_bc,))

soil = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
)

exp_tendency! = make_exp_tendency(soil)
imp_tendency! = make_imp_tendency(soil)
update_jacobian! = make_update_jacobian(soil)

rmses = Array{FT}(undef, length(dts))
mass_errors = Array{FT}(undef, length(dts))
for i in eachindex(dts)
    dt = dts[i]

    Y, p, coords = initialize(soil)
    @. Y.soil.ϑ_l = FT(0.24)

    jac_kwargs =
        (; jac_prototype = RichardsTridiagonalW(Y), Wfact = update_jacobian!)

    prob = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = ODE.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLSM.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )

    sol = ODE.solve(prob, ode_algo; dt = dt, saveat = dt)

    # Calculate water mass balance over entire simulation
    mass_end = sum(sol.u[end].soil.ϑ_l)
    mass_start = sum(sol.u[1].soil.ϑ_l)
    t_sim = sol.t[end] - sol.t[1]
    # Flux changes water content every timestep (assumes constant flux_in, flux_out)
    mass_change_exp = -(flux_in - flux_out) * t_sim
    mass_change_actual = mass_end - mass_start
    relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp
    @assert relerr < FT(1e-9)
    mass_errors[i] = relerr

    # Compute RMSE vs reference solution (found using dt = 1s)
    rmse_flux = rmse(ref_soln_flux, parent(sol.u[end].soil.ϑ_l))
    @assert rmse_flux < FT(1e-2)
    rmses[i] = rmse_flux
end

# Save flux BC mass conservation error and RMSE as artifact
savedir = joinpath(pkgdir(ClimaLSM), "experiments/standalone/Soil")
plt = Plots.plot()
plt_twin = twinx(plt)
Plots.plot!(
    plt,
    dts,
    rmses,
    label = "RMSE",
    color = "red",
    linewidth = 3,
    xlabel = "Δt (s)",
    ylabel = "RMSE ϑ",
    legend = :bottomleft,
    background_color_legend = nothing,
    xaxis = :log10,
    yaxis = :log10,
    xticks = dts,
    title = "RMSE and Water Conservation with Flux BCs",
)
Plots.plot!(
    plt_twin,
    dts,
    mass_errors,
    label = "Water mass error",
    color = "purple",
    linewidth = 3,
    ylabel = "|∑ϑ-∑ϑ(0)|/∫ΔFdt",
    xlabel = "",
    legend = :bottomright,
    background_color_legend = nothing,
)
Plots.savefig(joinpath(savedir, "water_conservation_flux.png"))

# Uncomment to recreate reference solution artifact (using small dt)
# soln_file = joinpath(savedir, "ref_soln_flux.csv")
# open((soln_file), "w") do io
#     writedlm(io, parent(sol.u[end].soil.ϑ_l), ',')
# end


# Perform simulation with Dirichlet boundary conditions
top_state_bc = MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
flux_out = FT(0)
bot_flux_bc = Soil.FluxBC((p, t) -> eltype(t)(flux_out))
boundary_conds =
    (; top = (water = top_state_bc,), bottom = (water = bot_flux_bc,))

soil_dirichlet = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_conds,
    sources = sources,
)

exp_tendency! = make_exp_tendency(soil_dirichlet)
imp_tendency! = make_imp_tendency(soil_dirichlet)
update_jacobian! = make_update_jacobian(soil_dirichlet)
update_aux! = make_update_aux(soil_dirichlet)

rmses_dirichlet = Array{FT}(undef, length(dts))
mass_errors_dirichlet = Array{FT}(undef, length(dts))
for i in eachindex(dts)
    dt = dts[i]

    Y, p, coords = initialize(soil_dirichlet)
    @. Y.soil.ϑ_l = FT(0.24)

    jac_kwargs =
        (; jac_prototype = RichardsTridiagonalW(Y), Wfact = update_jacobian!)

    prob = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = ODE.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLSM.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )
    sol = ODE.solve(prob, ode_algo; dt = dt, saveat = dt)

    # Calculate water mass balance over entire simulation
    # Convert Dirichlet BC to flux
    z = ClimaCore.Fields.coordinate_field(soil_domain.space).z
    Δz_top, Δz_bottom = get_Δz(z)

    times = collect(t_start:dt:t_end)
    flux_in_sim = Array{FT}(undef, length(times) - 1)
    # Because we use Backward Euler, compute fluxes at times[2:end]
    for j in 2:length(times)
        update_aux!(p, sol.u[j], times[j])
        top_flux_bc = boundary_flux(
            top_state_bc,
            TopBoundary(),
            soil_dirichlet,
            Δz_top,
            sol.u[j],
            p,
            times[j],
        )
        flux_in_sim[j - 1] = parent(top_flux_bc)[1]
    end

    mass_end = sum(sol.u[end].soil.ϑ_l)
    mass_start = sum(sol.u[1].soil.ϑ_l)
    t_sim = sol.t[end] - sol.t[1]
    mass_change_exp = -(sum(flux_in_sim) * dt - flux_out * t_sim)
    mass_change_actual = mass_end - mass_start
    relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp
    @assert relerr < FT(1e-5)
    mass_errors_dirichlet[i] = relerr

    # Compute RMSE vs reference solution (found using small dt = 1s)
    rmse_dirichlet = rmse(ref_soln_dirichlet, parent(sol.u[end].soil.ϑ_l))
    @assert rmse_dirichlet < FT(1e-2)
    rmses_dirichlet[i] = rmse_dirichlet
end

# Save Dirichlet BC mass conservation error and RMSE as artifact
plt = Plots.plot()
plt_twin = twinx(plt)
Plots.plot!(
    plt,
    dts,
    rmses_dirichlet,
    label = "RMSE",
    color = "red",
    linewidth = 3,
    xlabel = "Δt (s)",
    ylabel = "RMSE ϑ",
    legend = :bottomleft,
    background_color_legend = nothing,
    xaxis = :log10,
    yaxis = :log10,
    xticks = dts,
    title = "RMSE and Water Conservation with Dirichlet BCs",
)
Plots.plot!(
    plt_twin,
    dts,
    mass_errors_dirichlet,
    label = "Water mass error",
    color = "purple",
    linewidth = 3,
    ylabel = "|∑ϑ-∑ϑ(0)|/∫ΔFdt",
    xlabel = "",
    legend = :bottomright,
    background_color_legend = nothing,
)
Plots.savefig(joinpath(savedir, "water_conservation_dirichlet.png"))

# Uncomment to recreate Dirichlet BC reference solution artifact (using small dt)
# soln_file_dirichlet = joinpath(savedir, "ref_soln_dirichlet.csv")
# open((soln_file_dirichlet), "w") do io
#     writedlm(io, parent(sol.u[end].soil.ϑ_l), ',')
# end

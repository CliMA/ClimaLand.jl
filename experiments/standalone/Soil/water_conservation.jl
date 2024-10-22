using ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
import SciMLBase
import ClimaTimeSteppers as CTS
using Plots
using Statistics
using DelimitedFiles

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Domains: Column
import ClimaUtilities.OutputPathGenerator: generate_output_path

rmse(v1, v2) = sqrt(mean((v1 .- v2) .^ 2))

# Read in reference solutions from artifacts
flux_datapath, dirichlet_datapath =
    ClimaLand.Artifacts.water_conservation_test_data_path()
ref_soln_flux = readdlm(joinpath(flux_datapath, "ref_soln_flux.csv"))
ref_soln_dirichlet =
    readdlm(joinpath(dirichlet_datapath, "ref_soln_dirichlet.csv"))

# Define simulation times
t_start = Float64(0)
dt_min = Float64(10)
t_end = Float64(1e6)

for FT in (Float32, Float64)
    stepper = CTS.ARS111()
    # Select conv. condition based on float type due to different precision
    err = (FT == Float64) ? 1e-8 : 1e-4
    convergence_cond = CTS.MaximumError(err)
    conv_checker = CTS.ConvergenceChecker(norm_condition = convergence_cond)
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 50,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
            convergence_checker = conv_checker,
        ),
    )

    dts =
        (FT == Float64) ? Float64.([10, 100, 1000, 10000]) : Float64.([dt_min])

    # van Genuchten parameters for clay (from Bonan 2019 supplemental program 8.2)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    vg_n = FT(1.43)
    vg_α = FT(0.026 * 100) # inverse meters
    θ_r = FT(0.124)
    S_s = FT(1e-3) #inverse meters
    hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)

    zmax = FT(0)
    zmin = FT(-0.5)
    nelems = 50

    params = Soil.RichardsParameters(;
        ν = ν,
        hydrology_cm = hcm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    sources = ()

    # Set flux boundary conditions (used for calculating mass balance)
    flux_in = FT(-1e-7)
    top_bc = Soil.WaterFluxBC((p, t) -> flux_in)
    flux_out = FT(0)
    bot_bc = Soil.WaterFluxBC((p, t) -> flux_out)

    boundary_fluxes = (; top = top_bc, bottom = bot_bc)

    soil = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    exp_tendency! = make_exp_tendency(soil)
    set_initial_cache! = make_set_initial_cache(soil)
    imp_tendency! = make_imp_tendency(soil)
    jacobian! = make_jacobian(soil)

    rmses = Array{FT}(undef, length(dts))
    mass_errors = Array{FT}(undef, length(dts))
    for i in eachindex(dts)
        dt = dts[i]

        Y, p, coords = initialize(soil)
        @. Y.soil.ϑ_l = FT(0.24)
        set_initial_cache!(p, Y, FT(0.0))

        jac_kwargs =
            (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

        prob = SciMLBase.ODEProblem(
            CTS.ClimaODEFunction(
                T_exp! = exp_tendency!,
                T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
                dss! = ClimaLand.dss!,
            ),
            Y,
            (t_start, t_end),
            p,
        )

        sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = dt)

        # Check that simulation still has correct float type
        @assert eltype(sol.u[end].soil) == FT

        if FT == Float64
            # Calculate water mass balance over entire simulation
            mass_end = sum(sol.u[end].soil.ϑ_l)
            mass_start = sum(sol.u[1].soil.ϑ_l)
            t_sim = sol.t[end] - sol.t[1]
            # Flux changes water content every timestep (assumes constant flux_in, flux_out)
            mass_change_exp = -(flux_in - flux_out) * t_sim
            mass_change_actual = mass_end - mass_start
            relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp
            @show relerr
            @assert relerr < FT(1e-9)
            mass_errors[i] = relerr

            # Compute RMSE vs reference solution (found using dt = 1s)
            rmse_flux = rmse(ref_soln_flux, parent(sol.u[end].soil.ϑ_l))
            @assert rmse_flux < FT(1e-2)
            rmses[i] = rmse_flux
        end
    end

    if FT == Float64
        # Save flux BC mass conservation error and RMSE as artifact
        savedir = generate_output_path(
            "experiments/standalone/Soil/water_conservation",
        )
        plt = Plots.plot(margin = 10Plots.mm)
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
    end


    # Perform simulation with Dirichlet boundary conditions
    top_state_bc = Soil.MoistureStateBC((p, t) -> ν - 1e-3)
    flux_out = FT(0)
    bot_flux_bc = Soil.WaterFluxBC((p, t) -> flux_out)
    boundary_conds = (; top = top_state_bc, bottom = bot_flux_bc)

    soil_dirichlet = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_conds,
        sources = sources,
    )

    exp_tendency! = make_exp_tendency(soil_dirichlet)
    set_initial_cache! = make_set_initial_cache(soil_dirichlet)
    imp_tendency! = make_imp_tendency(soil_dirichlet)
    jacobian! = make_jacobian(soil_dirichlet)
    update_aux! = make_update_aux(soil_dirichlet)

    rmses_dirichlet = Array{FT}(undef, length(dts))
    mass_errors_dirichlet = Array{FT}(undef, length(dts))
    for i in eachindex(dts)
        dt = dts[i]

        Y, p, coords = initialize(soil_dirichlet)
        @. Y.soil.ϑ_l = FT(0.24)
        set_initial_cache!(p, Y, FT(0.0))

        jac_kwargs =
            (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

        prob = SciMLBase.ODEProblem(
            CTS.ClimaODEFunction(
                T_exp! = exp_tendency!,
                T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
                dss! = ClimaLand.dss!,
            ),
            Y,
            (t_start, t_end),
            p,
        )
        saveat = Array(t_start:dt:t_end)
        sv = (;
            t = Array{Int64}(undef, length(saveat)),
            saveval = Array{NamedTuple}(undef, length(saveat)),
        )
        cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
        sol = SciMLBase.solve(
            prob,
            ode_algo;
            dt = dt,
            callback = cb,
            adaptive = false,
            saveat = saveat,
        )
        # Check that simulation still has correct float type
        @assert eltype(sol.u[end].soil) == FT

        # Calculate water mass balance over entire simulation
        # Because we use Backward Euler, compute fluxes at times[2:end]
        flux_in_sim =
            [parent(sv.saveval[k].soil.top_bc)[1] for k in 2:length(sv.saveval)]


        mass_end = sum(sol.u[end].soil.ϑ_l)
        mass_start = sum(sol.u[1].soil.ϑ_l)
        t_sim = sol.t[end] - sol.t[1]
        mass_change_exp = -(sum(flux_in_sim) * dt - flux_out * t_sim)
        mass_change_actual = mass_end - mass_start
        relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp
        @assert relerr < 1e11 * eps(FT)
        mass_errors_dirichlet[i] = relerr

        # Compute RMSE vs reference solution (found using small dt = 1s)
        rmse_dirichlet = rmse(ref_soln_dirichlet, parent(sol.u[end].soil.ϑ_l))
        @assert rmse_dirichlet < 1e14 * eps(FT)
        rmses_dirichlet[i] = rmse_dirichlet
    end

    # Save Dirichlet BC mass conservation error and RMSE as artifact
    if FT == Float64
        plt = Plots.plot(margin = 10Plots.mm)
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
        #     writedlm(io, parent(sol.u[end].soil.ϑ_l), ',')
        # end
    end
end

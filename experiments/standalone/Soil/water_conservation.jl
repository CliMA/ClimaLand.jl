using ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
import SciMLBase
import ClimaTimeSteppers as CTS
using CairoMakie
using Statistics
using DelimitedFiles

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Domains: Column
import ClimaUtilities.OutputPathGenerator: generate_output_path
import ClimaLand.Simulations: LandSimulation, solve!

rmse(v1, v2) = sqrt(mean((v1 .- v2) .^ 2))

# Read in reference solutions from test artifacts
include("../../../test/Artifacts.jl")
flux_datapath, dirichlet_datapath = water_conservation_test_data_path()
ref_soln_flux = readdlm(flux_datapath)
ref_soln_dirichlet = readdlm(dirichlet_datapath)

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


    rmses = Array{FT}(undef, length(dts))
    mass_errors = Array{FT}(undef, length(dts))
    function set_ic!(Y, p, t0, model)
        @. Y.soil.ϑ_l = FT(0.24)
    end
    for i in eachindex(dts)
        dt = dts[i]
        simulation = LandSimulation(
            t_start,
            t_end,
            dt,
            soil;
            set_ic! = set_ic!,
            timestepper = ode_algo,
            solver_kwargs = (; saveat = collect(t_start:dt:t_end)),
        )
        p = simulation._integrator.p
        p_init = deepcopy(p)
        mass_start = p_init.soil.total_water
        sol = solve!(simulation)

        # Check that simulation still has correct float type
        @assert eltype(sol.u[end].soil) == FT

        if FT == Float64
            # Calculate water mass balance over entire simulation
            mass_end = p.soil.total_water
            ∫Fdt_end = sol.u[end].soil.∫F_vol_liq_water_dt
            ∫Fdt_start = sol.u[1].soil.∫F_vol_liq_water_dt
            mass_change_exp = Array(parent(∫Fdt_end .- ∫Fdt_start))[1]
            mass_change_actual = Array(parent(mass_end .- mass_start))[1]
            relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp
            @assert relerr < sqrt(eps(FT))
            mass_errors[i] = relerr

            # Compute RMSE vs reference solution (found using dt = 0.01s)
            rmse_flux = rmse(ref_soln_flux, parent(sol.u[end].soil.ϑ_l))
            @assert rmse_flux < 1e-2
            rmses[i] = rmse_flux
        end
    end

    if FT == Float64
        # Save flux BC mass conservation error and RMSE as artifact
        savedir = generate_output_path(
            "experiments/standalone/Soil/water_conservation",
        )
        fig = CairoMakie.Figure(
            title = "RMSE and Water Conservation with Flux BCs",
        )
        ax1 = Axis(
            fig[1, 1],
            xlabel = "Δt (s)",
            ylabel = "RMSE ϑ",
            xscale = log10,
            yscale = log10,
            xticks = dts,
        )
        ax2 = Axis(
            fig[1, 1],
            yaxisposition = :right,
            ylabel = "|∑ϑ-∑ϑ(0)|/∫ΔFdt",
            xscale = log10,
            yscale = log10,
            xticks = dts,
        )
        hidespines!(ax2)
        hidexdecorations!(ax2)

        l1 = lines!(
            ax1,
            dts,
            rmses,
            label = "RMSE",
            color = "red",
            linewidth = 3,
        )

        l2 = lines!(
            ax2,
            dts,
            mass_errors,
            label = "Water mass error",
            color = "purple",
            linewidth = 3,
        )

        axislegend(
            ax1,
            [l1, l2],
            ["RMSE", "Water mass error"],
            position = :rb,
            orientation = :vertical,
        )

        CairoMakie.save(joinpath(savedir, "water_conservation_flux.png"), fig)

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

    rmses_dirichlet = Array{FT}(undef, length(dts))
    mass_errors_dirichlet = Array{FT}(undef, length(dts))
    for i in eachindex(dts)
        dt = dts[i]

        function set_ic!(Y, p, t0, model)
            @. Y.soil.ϑ_l = FT(0.24)
        end
        simulation = LandSimulation(
            t_start,
            t_end,
            dt,
            soil_dirichlet;
            set_ic! = set_ic!,
            timestepper = ode_algo,
            solver_kwargs = (; saveat = collect(t_start:dt:t_end)),
        )
        p = simulation._integrator.p
        p_init = deepcopy(p)
        mass_start = p_init.soil.total_water
        sol = solve!(simulation)
        # Check that simulation still has correct float type
        @assert eltype(sol.u[end].soil) == FT

        mass_end = p.soil.total_water
        ∫Fdt_end = sol.u[end].soil.∫F_vol_liq_water_dt
        ∫Fdt_start = sol.u[1].soil.∫F_vol_liq_water_dt
        mass_change_exp = Array(parent(∫Fdt_end .- ∫Fdt_start))[1]
        mass_change_actual = Array(parent(mass_end .- mass_start))[1]
        relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp
        @assert relerr < 1e9 * eps(FT)
        mass_errors_dirichlet[i] = relerr

        # Compute RMSE vs reference solution (found using small dt = 1s)
        rmse_dirichlet = rmse(ref_soln_dirichlet, parent(sol.u[end].soil.ϑ_l))
        @assert rmse_dirichlet < 1e14 * eps(FT)
        rmses_dirichlet[i] = rmse_dirichlet
    end

    # Save Dirichlet BC mass conservation error and RMSE as artifact
    if FT == Float64
        fig = CairoMakie.Figure(
            title = "RMSE and Water Conservation with Dirichlet BCs",
        )
        ax1 = Axis(
            fig[1, 1],
            xlabel = "Δt (s)",
            ylabel = "RMSE ϑ",
            xscale = log10,
            yscale = log10,
            xticks = dts,
        )
        ax2 = Axis(
            fig[1, 1],
            yaxisposition = :right,
            ylabel = "|∑ϑ-∑ϑ(0)|/∫ΔFdt",
            xscale = log10,
            yscale = log10,
            xticks = dts,
        )
        hidespines!(ax2)
        hidexdecorations!(ax2)

        l1 = lines!(
            ax1,
            dts,
            rmses_dirichlet,
            label = "RMSE",
            color = "red",
            linewidth = 3,
        )

        l2 = lines!(
            ax2,
            dts,
            mass_errors_dirichlet,
            label = "Water mass error",
            color = "purple",
            linewidth = 3,
        )

        axislegend(
            ax1,
            [l1, l2],
            ["RMSE", "Water mass error"],
            position = :rb,
            orientation = :vertical,
        )

        CairoMakie.save(
            joinpath(savedir, "water_conservation_dirichlet.png"),
            fig,
        )
        # Uncomment to recreate Dirichlet BC reference solution artifact (using small dt)
        # soln_file_dirichlet = joinpath(savedir, "ref_soln_dirichlet.csv")
        # open((soln_file_dirichlet), "w") do io
        #     writedlm(io, parent(sol.u[end].soil.ϑ_l), ',')
        #     writedlm(io, parent(sol.u[end].soil.ϑ_l), ',')
        # end
    end
end

using Test
import ClimaComms
ClimaComms.@import_required_backends
using CairoMakie
using DelimitedFiles
using Statistics
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaUtilities.OutputPathGenerator: generate_output_path

import ClimaLand
import ClimaLand.Parameters as LP

bonan_data_folder = ClimaLand.Artifacts.richards_eqn_bonan_data_path()
clay_datapath = joinpath(bonan_data_folder, "bonan_data_clay.txt")
sand_datapath = joinpath(bonan_data_folder, "bonan_data_sand.txt")

context = ClimaComms.context()
ClimaComms.init(context)
device_suffix =
    typeof(context.device) <: ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = generate_output_path(
    joinpath("experiments", "standalone", "Soil", device_suffix),
)
!ispath(outdir) && mkpath(outdir)

@testset "Richards comparison to Bonan; clay" begin
    # Define simulation times
    t0 = Float64(0)
    dt = Float64(1e3)
    tf = Float64(1e6)

    for (FT, tf) in ((Float32, 2 * dt), (Float64, tf))
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(1.43)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.124)
        zmax = FT(0)
        zmin = FT(-1.5)
        nelems = 150
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z

        top_state_bc = MoistureStateBC((p, t) -> ν - 1e-3)
        bot_flux_bc = FreeDrainage()
        sources = ()
        boundary_states = (; top = top_state_bc, bottom = bot_flux_bc)
        params = Soil.RichardsParameters(;
            ν = ν,
            hydrology_cm = hcm,
            K_sat = K_sat,
            S_s = S_s,
            θ_r = θ_r,
        )

        soil = Soil.RichardsModel{FT}(;
            parameters = params,
            domain = soil_domain,
            boundary_conditions = boundary_states,
            sources = sources,
        )
        set_initial_cache! = make_set_initial_cache(soil)

        Y, p, coords = initialize(soil)

        # specify ICs
        Y.soil.ϑ_l .= FT(0.24)
        exp_tendency! = make_exp_tendency(soil)
        imp_tendency! = ClimaLand.make_imp_tendency(soil)
        jacobian! = ClimaLand.make_jacobian(soil)
        set_initial_cache!(p, Y, t0)

        stepper = CTS.ARS111()
        norm_condition = CTS.MaximumError(FT(1e-8))
        conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
        ode_algo = CTS.IMEXAlgorithm(
            stepper,
            CTS.NewtonsMethod(
                max_iters = 50,
                update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
                convergence_checker = conv_checker,
            ),
        )

        # set up jacobian info
        jac_kwargs =
            (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

        prob = SciMLBase.ODEProblem(
            CTS.ClimaODEFunction(
                T_exp! = exp_tendency!,
                T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
                dss! = ClimaLand.dss!,
            ),
            Y,
            (t0, tf),
            p,
        )
        sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = 10000)

        # Check that simulation still has correct float type
        @assert eltype(sol.u[end].soil) == FT

        # Check results and plot for Float64 simulation on CPU
        if FT == Float64 && context.device isa ClimaComms.CPUSingleThreaded
            N = length(sol.t)
            ϑ_l = parent(sol.u[N].soil.ϑ_l)
            ds_bonan = readdlm(clay_datapath)
            bonan_moisture = reverse(ds_bonan[:, 1])
            bonan_z = reverse(ds_bonan[:, 2]) ./ 100.0
            @test sqrt.(mean((bonan_moisture .- ϑ_l) .^ 2.0)) < FT(1e-3)
            fig = CairoMakie.Figure()
            ax = CairoMakie.Axis(fig[1, 1])
            lines!(ax, ϑ_l[:], parent(z)[:], label = "Clima")
            lines!(ax, bonan_moisture, bonan_z, label = "Bonan's Matlab code")
            axislegend()
            CairoMakie.save(
                joinpath(outdir, "comparison_clay_bonan_matlab.png"),
                fig,
            )
        end
    end
end


@testset "Richards comparison to Bonan; sand" begin
    # Define simulation times
    # Note, we can use a bigger step and still conserve mass.
    t0 = Float64(0)
    dt = Float64(1)
    tf = Float64(60 * 60 * 0.8)

    for (FT, tf) in ((Float32, 2 * dt), (Float64, tf))
        ν = FT(0.287)
        K_sat = FT(34 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(3.96)
        vg_α = FT(2.7) # inverse meters
        vg_m = FT(1)
        S_c = ν + 1
        hcm = vanGenuchten(vg_α, vg_n, vg_m, S_c)
        θ_r = FT(0.075)
        zmax = FT(0)
        zmin = FT(-1.5)
        nelems = 150
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z

        top_state_bc = MoistureStateBC((p, t) -> 0.267)
        bot_flux_bc = FreeDrainage()
        sources = ()
        boundary_states = (; top = top_state_bc, bottom = bot_flux_bc)

        params =
            Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)

        soil = Soil.RichardsModel{FT}(;
            parameters = params,
            domain = soil_domain,
            boundary_conditions = boundary_states,
            sources = sources,
        )

        Y, p, coords = initialize(soil)
        set_initial_cache! = make_set_initial_cache(soil)

        # specify ICs
        Y.soil.ϑ_l .= FT(0.1)
        exp_tendency! = make_exp_tendency(soil)
        imp_tendency! = ClimaLand.make_imp_tendency(soil)
        jacobian! = ClimaLand.make_jacobian(soil)
        set_initial_cache!(p, Y, t0)

        stepper = CTS.ARS111()
        norm_condition = CTS.MaximumError(FT(1e-8))
        conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
        ode_algo = CTS.IMEXAlgorithm(
            stepper,
            CTS.NewtonsMethod(
                max_iters = 50,
                update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
                convergence_checker = conv_checker,
            ),
        )
        # set up jacobian info
        jac_kwargs =
            (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

        prob = SciMLBase.ODEProblem(
            CTS.ClimaODEFunction(
                T_exp! = exp_tendency!,
                T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
                dss! = ClimaLand.dss!,
            ),
            Y,
            (t0, tf),
            p,
        )
        sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = 60 * dt)

        # Check that simulation still has correct float type
        @assert eltype(sol.u[end].soil) == FT

        # Check results and plot for Float64 simulation on CPU
        if FT == Float64 && context.device isa ClimaComms.CPUSingleThreaded
            N = length(sol.t)
            ϑ_l = parent(sol.u[N].soil.ϑ_l)
            ds_bonan = readdlm(sand_datapath)
            bonan_moisture = reverse(ds_bonan[:, 1])
            bonan_z = reverse(ds_bonan[:, 2]) ./ 100.0
            @test sqrt.(mean((bonan_moisture .- ϑ_l) .^ 2.0)) < FT(1e-3)
            fig = CairoMakie.Figure()
            ax = CairoMakie.Axis(fig[1, 1])
            lines!(ax, ϑ_l[:], parent(z)[:], label = "Clima")
            lines!(ax, bonan_moisture, bonan_z, label = "Bonan's Matlab code")
            axislegend()
            CairoMakie.save(
                joinpath(outdir, "comparison_sand_bonan_matlab.png"),
                fig,
            )
        end
    end
end

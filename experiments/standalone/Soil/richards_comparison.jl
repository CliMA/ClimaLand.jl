using Test
using Plots
using DelimitedFiles
using Statistics
using ArtifactWrappers
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil

import ClimaLSM
import ClimaLSM.Parameters as LSMP

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

@testset "Richards comparison to Bonan; clay" begin
    FT = Float64
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(1.43)
    vg_α = FT(2.6) # inverse meters
    hcm = vanGenuchten(; α = vg_α, n = vg_n)
    θ_r = FT(0.124)
    zmax = FT(0)
    zmin = FT(-1.5)
    nelems = 150
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    z = ClimaCore.Fields.coordinate_field(soil_domain.space).z

    top_state_bc = MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
    bot_flux_bc = FreeDrainage()
    sources = ()
    boundary_states =
        (; top = (water = top_state_bc,), bottom = (water = bot_flux_bc,))
    params = Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)

    soil = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_states,
        sources = sources,
    )
    set_initial_aux_state! = make_set_initial_aux_state(soil)

    Y, p, coords = initialize(soil)

    # specify ICs
    Y.soil.ϑ_l .= FT(0.24)
    exp_tendency! = make_exp_tendency(soil)
    imp_tendency! = ClimaLSM.make_imp_tendency(soil)
    update_jacobian! = ClimaLSM.make_update_jacobian(soil)

    t0 = FT(0)
    set_initial_aux_state!(p, Y, t0)
    tf = FT(1e6)
    dt = FT(1e3)

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
        (; jac_prototype = RichardsTridiagonalW(Y), Wfact = update_jacobian!)

    prob = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = ODE.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLSM.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    sol = ODE.solve(prob, ode_algo; dt = dt, saveat = 10000)

    N = length(sol.t)
    ϑ_l = parent(sol.u[N].soil.ϑ_l)
    bonan_clay_dataset = ArtifactWrapper(
        @__DIR__,
        "richards_clay",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/nk89znth59gcsdb65lnywnzjnuno3h6k.txt",
            filename = "clay_bonan_sp801_22323.txt",
        ),],
    )
    datapath = get_data_folder(bonan_clay_dataset)
    data = joinpath(datapath, "clay_bonan_sp801_22323.txt")
    ds_bonan = readdlm(data)
    bonan_moisture = reverse(ds_bonan[:, 1])
    bonan_z = reverse(ds_bonan[:, 2]) ./ 100.0
    @test sqrt.(mean((bonan_moisture .- ϑ_l) .^ 2.0)) < FT(1e-3)
    plot(ϑ_l, parent(z), label = "Clima")
    plot!(bonan_moisture, bonan_z, label = "Bonan's Matlab code")
    savefig("./experiments/standalone/Soil/comparison_clay_bonan_matlab.png")
end


@testset "Richards comparison to Bonan; sand" begin
    FT = Float64
    ν = FT(0.287)
    K_sat = FT(34 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(3.96)
    vg_α = FT(2.7) # inverse meters
    hcm = vanGenuchten(; α = vg_α, n = vg_n)
    θ_r = FT(0.075)
    zmax = FT(0)
    zmin = FT(-1.5)
    nelems = 150
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    z = ClimaCore.Fields.coordinate_field(soil_domain.space).z

    top_state_bc = MoistureStateBC((p, t) -> eltype(t)(0.267))
    bot_flux_bc = FreeDrainage()
    sources = ()
    boundary_states =
        (; top = (water = top_state_bc,), bottom = (water = bot_flux_bc,))

    params = Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)

    soil = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_states,
        sources = sources,
    )

    Y, p, coords = initialize(soil)
    set_initial_aux_state! = make_set_initial_aux_state(soil)

    # specify ICs
    Y.soil.ϑ_l .= FT(0.1)
    exp_tendency! = make_exp_tendency(soil)
    imp_tendency! = ClimaLSM.make_imp_tendency(soil)
    update_jacobian! = ClimaLSM.make_update_jacobian(soil)

    t0 = FT(0)
    set_initial_aux_state!(p, Y, t0)
    tf = FT(60 * 60 * 0.8)
    dt = FT(1)
    # Note, we can use a bigger step and still conserve mass.

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
        (; jac_prototype = RichardsTridiagonalW(Y), Wfact = update_jacobian!)

    prob = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = ODE.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLSM.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    sol = ODE.solve(prob, ode_algo; dt = dt, saveat = 60 * dt)

    N = length(sol.t)
    ϑ_l = parent(sol.u[N].soil.ϑ_l)
    bonan_sand_dataset = ArtifactWrapper(
        @__DIR__,
        "richards_sand",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/2vk7bvyjah8xd5b7wxcqy72yfd2myjss.csv",
            filename = "sand_bonan_sp801.csv",
        ),],
    )
    datapath = get_data_folder(bonan_sand_dataset)
    data = joinpath(datapath, "sand_bonan_sp801.csv")
    ds_bonan = readdlm(data, ',')
    bonan_moisture = reverse(ds_bonan[:, 1])
    bonan_z = reverse(ds_bonan[:, 2]) ./ 100.0
    @test sqrt.(mean((bonan_moisture .- ϑ_l) .^ 2.0)) < FT(1e-3)
    plot(ϑ_l, parent(z), label = "Clima")
    plot!(bonan_moisture, bonan_z, label = "Bonan's Matlab code")
    savefig("./experiments/standalone/Soil/comparison_sand_bonan_matlab.png")
end

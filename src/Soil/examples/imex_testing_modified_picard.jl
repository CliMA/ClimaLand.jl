using ClimaCore
using DiffEqBase
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: Column
include("./TridiagonalJacobian.jl")
using .TridiagonalJacobian: TridiagonalW, make_Wfact, make_implicit_tendency, explicit_tendency!

is_imex_CTS_algo(::CTS.IMEXAlgorithm) = true
is_imex_CTS_algo(::DiffEqBase.AbstractODEAlgorithm) = false

is_implicit(::ODE.OrdinaryDiffEqImplicitAlgorithm) = true
is_implicit(::ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
is_implicit(ode_algo) = is_imex_CTS_algo(ode_algo)

is_rosenbrock(::ODE.Rosenbrock23) = true
is_rosenbrock(::ODE.Rosenbrock32) = true
is_rosenbrock(::DiffEqBase.AbstractODEAlgorithm) = false
use_transform(ode_algo) =
    !(is_imex_CTS_algo(ode_algo) || is_rosenbrock(ode_algo))
stepper = CTS.ARS111()
norm_condition = CTS.MaximumError(
    Float64(1e-6),
)
conv_checker = CTS.ConvergenceChecker(; norm_condition)
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 500,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
        verbose = CTS.Verbose(),
    ),
)

#function main(ode_algo, t_end::Float64, dt::Float64; explicit=false)

    FT = Float64


    t_end = FT(130000)
    dt = FT(1000)
    # t_end = 1 * dt
    explicit = false
    # parameters for clay from Bonan 2019 supplemental program 8.2
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    vg_n = FT(1.43)
    vg_α = FT(0.026 * 100) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.124)
    S_s = FT(1e-3) #inverse meters

    zmax = FT(0)
    zmin = FT(-1.5)
    nelems = 150

    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

    top_bc = Soil.MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
    bot_bc = Soil.FreeDrainage()
    sources = ()
    boundary_fluxes = (; top = (water = top_bc,), bottom = (water = bot_bc,))
    params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

    soil = Soil.RichardsModel{FT}(;
                                  parameters = params,
                                  domain = soil_domain,
                                  boundary_conditions = boundary_fluxes,
                                  sources = sources,
                                  )

    Y, p, coords = initialize(soil)
    @. Y.soil.ϑ_l = FT(0.24)
    # update_aux! = make_update_aux(soil)
    t_start = FT(0)
    # update_aux!(p, Y, t_start)

    @. p.soil.K = FT(NaN)
    @. p.soil.ψ = FT(NaN)

    if !explicit
        #@show "in if statement"
        #@show Y.soil.ϑ_l
        #@show p.soil.K
        #@show p.soil.ψ
        # #@show W.temp1
        # #@show W.temp2

        transform = use_transform(ode_algo)

        W = TridiagonalW(Y, transform)
        #@show "after TridiagonalW"
        #@show Y.soil.ϑ_l
        #@show p.soil.K
        #@show p.soil.ψ
        #@show W.temp1
        #@show W.temp2

        implicit_tendency! = make_implicit_tendency(soil)
        Wfact! = make_Wfact(soil)
        Wfact!(W, Y, p, dt, t_start)
        #@show "after init Wfact!"
        #@show Y.soil.ϑ_l
        #@show p.soil.K
        #@show p.soil.ψ
        jac_kwargs = if use_transform(ode_algo)
            (; jac_prototype = W, Wfact_t = Wfact!)
        else
            (; jac_prototype = W, Wfact = Wfact!)
        end

        implicit_problem = ODEProblem(
            CTS.ClimaODEFunction(
                T_exp! = explicit_tendency!,
                T_imp! = ODEFunction(implicit_tendency!; jac_kwargs...),
                dss! = ClimaLSM.Soil.dss!,
            ),
            Y,
            (t_start, t_end),
            p,
        )

        #@show "after implicit problem"
        #@show Y.soil.ϑ_l
        #@show p.soil.K
        #@show p.soil.ψ
        #@show implicit_problem.p
        #@show implicit_problem.u0.soil.ϑ_l

        integrator = init(
            implicit_problem,
            ode_algo;
            dt = dt,
            adaptive = false,
            progress = true,
            saveat = t_start:1:t_end,
        )
        #@show "after integrator init"
        #@show Y.soil.ϑ_l
        #@show p.soil.K
        #@show p.soil.ψ
        #@show integrator.sol.u[1].soil.ϑ_l
    else
        implicit_tendency! = make_implicit_tendency(soil)
        explicit_problem = ODE.ODEProblem(
            CTS.ClimaODEFunction(
                T_exp! = implicit_tendency!,
                T_imp! = nothing,
               dss! = ClimaLSM.Soil.dss!,
            ),
            Y,
            (t_start, t_end),
            p,
        )
        integrator = init(
            explicit_problem,
            ode_algo;
            dt = dt,
            adaptive = false,
            progress = true,
            saveat = t_start:1000:t_end,
        )
    end
    #@show "before stepping loop"
    #@show Y.soil.ϑ_l
    #@show p.soil.K
    #@show p.soil.ψ
    #@show integrator.sol.u[1].soil.ϑ_l
    for step in 1:118
        #@show step
        ODE.step!(integrator)
        #@show "end of step $step"
        #@show integrator.sol.u[1].soil.ϑ_l
    end
    # ODE.solve!(integrator)

    # sol = ODE.solve(implicit_problem, ode_algo; dt = dt, adaptive = false)

#end

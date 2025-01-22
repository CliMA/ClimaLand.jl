using ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
import SciMLBase
import ClimaTimeSteppers as CTS
using CairoMakie
using Statistics
using DelimitedFiles
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Domains: Column
import ClimaUtilities.OutputPathGenerator: generate_output_path

rmse(v1, v2) = sqrt(mean((v1 .- v2) .^ 2))
FT = Float64
ν = FT(0.395)
ν_ss_quartz = FT(0.92)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.0)
Ksat = FT(4.42 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.89)
vg_α = FT(7.5) # inverse meters
hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
θ_r = FT(0.0)
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
zmin = FT(-1.0)
nelems = 50
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);
surface_water_flux = WaterFluxBC((p, t) -> 0.0)
bottom_water_flux = WaterFluxBC((p, t) -> 0.0);

surface_heat_flux = HeatFluxBC((p, t) -> 0.0)
bottom_heat_flux = HeatFluxBC((p, t) -> 0.0);

boundary_fluxes = (;
    top = WaterHeatBC(; water = surface_water_flux, heat = surface_heat_flux),
    bottom = WaterHeatBC(; water = bottom_water_flux, heat = bottom_heat_flux),
);

sources = (PhaseChange{FT}(),);

soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
);

exp_tendency! = make_exp_tendency(soil);
imp_tendency! = make_imp_tendency(soil);
jacobian! = ClimaLand.make_jacobian(soil);


function init_soil!(Y, z, Trange, params)
    ν = params.ν
    θ_r = params.θ_r
    FT = eltype(Y.soil.ϑ_l)
    zmax = FT(0)
    zmin = FT(-1)

    theta_max = FT(ν * 0.5)
    theta_min = FT(ν * 0.4)
    (T_max, T_min) = Trange

    c = FT(20.0)
    @. Y.soil.ϑ_l =
        theta_min +
        (theta_max - theta_min) * exp(-(z - zmax) / (zmin - zmax) * c)
    Y.soil.θ_i .= FT(0.0)

    T = @.(T_min + (T_max - T_min) * exp(-(z - zmax) / (zmin - zmax) * c))

    θ_l = Soil.volumetric_liquid_fraction.(Y.soil.ϑ_l, ν, θ_r)
    ρc_s =
        Soil.volumetric_heat_capacity.(
            θ_l,
            Y.soil.θ_i,
            params.ρc_ds,
            params.earth_param_set,
        )
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            params.earth_param_set,
        )
end
set_initial_cache! = make_set_initial_cache(soil);

stepper = CTS.ARS111()
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


no_phase_change = (;
    t0 = Float64(0),
    tf = Float64(60 * 60 * 24),
    dt_ref = Float64(10),
    Trange = (FT(289.0), FT(288)),
    dts = Float64.([100, 500, 1000]),
    exp_name = "no_pc",
)
phase_change = (;
    t0 = Float64(0),
    tf = Float64(60 * 60 * 24),
    dt_ref = Float64(0.1),
    Trange = (FT(269.0), FT(270.0)),
    dts = Float64.([1, 10, 100]),
    exp_name = "pc",
)

savedir = generate_output_path(
    "experiments/standalone/Soil/water_energy_conservation",
)

for experiment in [no_phase_change, phase_change]
    (; t0, tf, dt_ref, Trange, dts, exp_name) = experiment
    Y, p, coords = initialize(soil)
    init_soil!(Y, coords.subsurface.z, Trange, soil.parameters)
    set_initial_cache!(p, Y, t0)
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
    ref_sol = SciMLBase.solve(prob, ode_algo; dt = dt_ref, saveat = [t0, tf])

    # Now try several dts
    rmses_water = Array{FT}(undef, length(dts))
    rmses_energy = Array{FT}(undef, length(dts))
    mass_errors = Array{FT}(undef, length(dts))
    energy_errors = Array{FT}(undef, length(dts))
    for i in eachindex(dts)
        dt = dts[i]

        Y, p, coords = initialize(soil)
        init_soil!(Y, coords.subsurface.z, Trange, soil.parameters)
        set_initial_cache!(p, Y, t0)
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

        sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = [t0, tf])
        # Calculate water mass balance over entire simulation
        mass_end = sum(sol.u[end].soil.ϑ_l .+ sol.u[end].soil.θ_i)
        mass_start = sum(sol.u[1].soil.ϑ_l .+ sol.u[1].soil.θ_i)
        t_sim = sol.t[end] - sol.t[1]
        # We used zero flux BC, so we expect no mass change.
        mass_change_exp = FT(0)
        mass_change_actual = abs(mass_end - mass_start) + eps(FT)
        relerr = abs(mass_change_actual) / mass_start
        mass_errors[i] = relerr

        # Calculate energy balance over entire simulation
        energy_end = sum(sol.u[end].soil.ρe_int)
        energy_start = sum(sol.u[1].soil.ρe_int)
        # We used zero flux BC, so we expect no energy change.
        energy_change_exp = FT(0)
        energy_change_actual = abs(energy_end - energy_start) + eps(FT)
        relerr = abs(energy_change_actual) / abs(energy_start)
        energy_errors[i] = relerr

        # Compute RMSE vs reference solution (found using dt_ref)
        rmse_energy =
            rmse(
                parent(ref_sol.u[end].soil.ρe_int),
                parent(sol.u[end].soil.ρe_int),
            ) / mean(abs.(parent(sol.u[1].soil.ρe_int)))
        rmses_energy[i] = rmse_energy

        rmse_water =
            rmse(parent(ref_sol.u[end].soil.ϑ_l), parent(sol.u[end].soil.ϑ_l)) /
            mean(parent(ref_sol.u[end].soil.ϑ_l))
        rmses_water[i] = rmse_water
    end

    # Save flux BC mass conservation error and RMSE as artifact
    fig = CairoMakie.Figure(
        title = "RMSE and Water Conservation with Zero Flux BC",
    )
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Δt (s)",
        ylabel = "RMSE(θ_l)/θ̄_l",
        xscale = log10,
        yscale = log10,
        xticks = dts,
    )
    ax2 = Axis(
        fig[1, 1],
        yaxisposition = :right,
        ylabel = "Error",
        xscale = log10,
        yscale = log10,
        xticks = dts,
    )
    hidespines!(ax2)
    hidexdecorations!(ax2)

    l1 = lines!(
        ax1,
        dts,
        rmses_water,
        label = "Fractional RMSE θ_l",
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
        ["Fractional RMSE θ_l", "Water mass error"],
        position = :rb,
        orientation = :vertical,
    )

    CairoMakie.save(
        joinpath(savedir, "water_conservation_flux_full_soil_$(exp_name).png"),
        fig,
    )

    fig = CairoMakie.Figure(
        title = "RMSE and Energy Conservation with Zero Flux BC",
    )
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Δt (s)",
        ylabel = "RMSE(ρe_int)/ρ̄e_int",
        xscale = log10,
        yscale = log10,
        xticks = dts,
    )
    ax2 = Axis(
        fig[1, 1],
        yaxisposition = :right,
        ylabel = "Error",
        xscale = log10,
        yscale = log10,
        xticks = dts,
    )
    hidespines!(ax2)
    hidexdecorations!(ax2)

    l1 = lines!(
        ax1,
        dts,
        rmses_energy,
        label = "Fractional RMSE ρe_int",
        color = "red",
        linewidth = 3,
    )

    l2 = lines!(
        ax2,
        dts,
        energy_errors,
        label = "Energy error",
        color = "purple",
        linewidth = 3,
    )

    axislegend(
        ax1,
        [l1, l2],
        ["Fractional RMSE ρe_int", "Energy error"],
        position = :rb,
        orientation = :vertical,
    )

    CairoMakie.save(
        joinpath(savedir, "energy_conservation_flux_full_soil_$(exp_name).png"),
        fig,
    )

end

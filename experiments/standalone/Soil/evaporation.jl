using Plots
import SciMLBase
import ClimaTimeSteppers as CTS
using Thermodynamics

using ClimaCore
import CLIMAParameters as CP
using RootSolvers
using SurfaceFluxes
using StaticArrays
using Dates
using ArtifactWrappers
using DelimitedFiles: readdlm

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
import SurfaceFluxes.Parameters as SFP

# Define simulation times
t0 = Float64(0)
tf = Float64(24 * 3600 * 13)
dt = Float64(2)

for (FT, tf) in ((Float32, 2 * dt), (Float64, tf))
    earth_param_set = LP.LandParameters(FT)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    # Coarse sand experiment described in Figures 7 and 8a
    # of Lehmann, Assouline, Or  (Phys Rev E 77, 2008)
    K_sat = FT(225.1 / 3600 / 24 / 1000)
    # n and alpha estimated by matching vG curve.
    vg_n = FT(10.0)
    vg_α = FT(6.0)
    hcm = vanGenuchten(; α = vg_α, n = vg_n)
    # Alternative parameters for Brooks Corey water retention model
    #ψb = FT(-0.14)
    #c = FT(5.5)
    #hcm = BrooksCorey(; ψb = ψb, c = c);
    ν = FT(0.43)
    θ_r = FT(0.045)
    S_s = FT(1e-3)
    ν_ss_om = FT(0.0)
    ν_ss_quartz = FT(1.0)
    ν_ss_gravel = FT(0.0)
    emissivity = FT(1.0)
    PAR_albedo = FT(0.2)
    NIR_albedo = FT(0.4)
    z_0m = FT(1e-3)
    z_0b = FT(1e-4)
    d_ds = FT(0.005)

    ref_time = DateTime(2005)
    SW_d = (t) -> 0
    LW_d = (t) -> 301.15^4 * 5.67e-8
    radiation = PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(SW_d),
        TimeVaryingInput(LW_d),
        ref_time,
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
    precip = (t) -> 0.0
    T_atmos = (t) -> T_air
    u_atmos = (t) -> 1
    q_atmos = (t) -> q
    h_atmos = FT(0.1)
    P_atmos = (t) -> 101325
    gustiness = FT(1e-2)
    atmos = PrescribedAtmosphere(
        TimeVaryingInput(precip),
        TimeVaryingInput(precip),
        TimeVaryingInput(T_atmos),
        TimeVaryingInput(u_atmos),
        TimeVaryingInput(q_atmos),
        TimeVaryingInput(P_atmos),
        ref_time,
        h_atmos;
        gustiness = gustiness,
    )
    top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(atmos, radiation)
    zero_flux = FluxBC((p, t) -> 0)
    boundary_fluxes =
        (; top = top_bc, bottom = (water = zero_flux, heat = zero_flux))
    params = ClimaLand.Soil.EnergyHydrologyParameters{FT}(;
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        hydrology_cm = hcm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
        PAR_albedo = PAR_albedo,
        NIR_albedo = NIR_albedo,
        emissivity = emissivity,
        z_0m = z_0m,
        z_0b = z_0b,
        earth_param_set = earth_param_set,
        d_ds = d_ds,
    )

    #TODO: Run with higher resolution once we have the implicit stepper
    zmax = FT(0)
    zmin = FT(-0.35)
    nelems = 5
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z

    soil = Soil.EnergyHydrology{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = (),
    )

    Y, p, cds = initialize(soil) # begins saturated
    function init_soil!(Y, z, params)
        ν = params.ν
        FT = eltype(ν)
        Y.soil.ϑ_l .= ν - 1e-2
        Y.soil.θ_i .= 0
        T = FT(301.15)
        ρc_s = Soil.volumetric_heat_capacity(ν, FT(0), params)
        Y.soil.ρe_int =
            Soil.volumetric_internal_energy.(FT(0), ρc_s, T, Ref(params))
    end

    init_soil!(Y, z, soil.parameters)

    # We also set the initial conditions of the auxiliary state here:
    set_initial_cache! = make_set_initial_cache(soil)
    set_initial_cache!(p, Y, t0)

    # Timestepping:
    soil_exp_tendency! = make_exp_tendency(soil)
    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = soil_exp_tendency!,
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    saveat = Array(t0:3600.0:tf)
    sv = (;
        t = Array{Float64}(undef, length(saveat)),
        saveval = Array{NamedTuple}(undef, length(saveat)),
    )
    saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
    updateat = deepcopy(saveat)
    updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
    sol =
        SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)

    # Check that simulation still has correct float type
    @assert eltype(sol.u[end].soil) == FT

    # Post processing for Float64 simulation
    if FT == Float64
        (; ν, θ_r, d_ds) = soil.parameters
        _D_vapor = FT(LP.D_vapor(soil.parameters.earth_param_set))
        evap = [
            parent(sv.saveval[k].soil.turbulent_fluxes.vapor_flux)[1] for
            k in 1:length(sol.t)
        ]
        r_ae = [
            parent(sv.saveval[k].soil.turbulent_fluxes.r_ae)[1] for
            k in 1:length(sol.t)
        ]
        T_soil = [parent(sv.saveval[k].soil.T)[end] for k in 1:length(sol.t)]
        S_c = hcm.S_c
        N = length(sol.t)
        potential_evap = zeros(N)
        r_soil = zeros(N)
        dsl = zeros(N)
        q_soil = zeros(N)
        surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)

        for i in 1:length(sol.t)
            time = sol.t[i]
            u = sol.u[i]
            p = sv.saveval[i]
            τ_a = ClimaLand.Domains.top_center_to_surface(
                @. max(eps(FT), (ν - p.soil.θ_l - u.soil.θ_i)^(FT(5 / 2)) / ν)
            )
            S_l_sfc = ClimaLand.Domains.top_center_to_surface(
                effective_saturation.(ν, u.soil.ϑ_l, θ_r),
            )
            layer_thickness = Soil.dry_soil_layer_thickness.(S_l_sfc, S_c, d_ds)
            T_sfc = T_soil[i]
            ts_in = ClimaLand.construct_atmos_ts(top_bc.atmos, p, thermo_params)
            ρ_sfc = compute_ρ_sfc.(thermo_params, ts_in, T_sfc)
            q_sat =
                Thermodynamics.q_vap_saturation_generic.(
                    thermo_params,
                    T_sfc,
                    ρ_sfc,
                    Thermodynamics.Liquid(),
                )
            q_sfc =
                ClimaLand.surface_specific_humidity(soil, Y, p, T_sfc, ρ_sfc)
            ts_sfc =
                Thermodynamics.PhaseEquil_ρTq.(
                    thermo_params,
                    ρ_sfc,
                    T_sfc,
                    q_sfc,
                )
            state_sfc =
                SurfaceFluxes.StateValues.(
                    0.0,
                    Ref(SVector{2, FT}(0, 0)),
                    ts_sfc,
                )
            state_in =
                SurfaceFluxes.StateValues.(
                    h_atmos,
                    Ref(SVector{2, FT}(u_atmos(time), 0)),
                    ts_in,
                )

            # State containers
            sc = SurfaceFluxes.ValuesOnly.(state_in, state_sfc, z_0m, z_0b)
            potential_conditions =
                SurfaceFluxes.surface_conditions.(
                    surface_flux_params,
                    sc;
                    tol_neutral = SFP.cp_d(surface_flux_params) / 100000,
                )
            vapor_flux =
                SurfaceFluxes.evaporation.(
                    surface_flux_params,
                    sc,
                    potential_conditions.Ch,
                ) ./ 1000.0

            # Store the values we want to plot later.
            # Because these are Fields, convert to scalars using `parent`
            dsl[i] = parent(layer_thickness)[1]
            r_soil[i] = parent(@. layer_thickness / (_D_vapor * τ_a))[1]
            q_soil[i] = parent(q_sfc)[1]
            potential_evap[i] = parent(vapor_flux)[1]
        end

        plt1 = Plots.plot()
        Plots.plot!(
            plt1,
            sol.t ./ 3600 ./ 24,
            evap ./ potential_evap,
            xlabel = "Days",
            ylabel = "E/E₀",
            label = "",
            margins = 6Plots.mm,
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
            margins = 6Plots.mm,
        )

        plt3 = Plots.plot()
        Plots.plot!(
            plt3,
            sol.t ./ 3600 ./ 24,
            r_soil,
            xlabel = "Days",
            ylabel = "R_soil (m/s)",
            label = "",
            margins = 6Plots.mm,
        )

        plt4 = Plots.plot()
        Plots.plot!(
            plt4,
            sol.t ./ 3600 ./ 24,
            evap .* (1000 * 3600 * 24),
            xlabel = "Days",
            ylabel = "E (mm/d)",
            label = "",
            margins = 6Plots.mm,
        )
        plt5 = Plots.plot()
        Plots.plot!(
            plt5,
            sol.t ./ 3600 ./ 24,
            T_soil,
            xlabel = "Days",
            ylabel = "T_sfc (K)",
            label = "",
            margins = 6Plots.mm,
        )

        plt6 = Plots.plot()
        Plots.plot!(
            plt6,
            sol.t ./ 3600 ./ 24,
            q_soil,
            xlabel = "Days",
            ylabel = "q_sfc",
            label = "",
            margins = 6Plots.mm,
        )
        top = [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:length(sol.t)]
        S = Soil.effective_saturation.(ν, top, θ_r)
        ψ = Soil.matric_potential.(hcm, S)
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
            margins = 6Plots.mm,
        )

        plt8 = Plots.plot()
        Plots.plot!(
            plt8,
            sol.t ./ 3600 ./ 24,
            ψ,
            xlabel = "Days",
            ylabel = "ψ_sfc",
            label = "",
            margins = 6Plots.mm,
        )

        savepath = joinpath(pkgdir(ClimaLand), "experiments/standalone/Soil")
        Plots.plot(plt1, plt2, plt3, plt4; layout = (2, 2))
        Plots.savefig(joinpath(savepath, "evaporation_from_coarse_sand1.png"))

        Plots.plot(plt5, plt6, plt7, plt8; layout = (2, 2))
        Plots.savefig(joinpath(savepath, "evaporation_from_coarse_sand2.png"))

        # Read in reference solution from artifact
        evap_dataset = ArtifactWrapper(
            @__DIR__,
            "lehmann2008_fig8_evaporation",
            ArtifactFile[ArtifactFile(
                url = "https://caltech.box.com/shared/static/cgppw3tx6zdz7h02yt28ri44g1j088ju.csv",
                filename = "lehmann2008_fig8_evaporation.csv",
            ),],
        )
        evap_datapath = get_data_folder(evap_dataset)
        ref_soln_E = readdlm(
            joinpath(evap_datapath, "lehmann2008_fig8_evaporation.csv"),
            ',',
        )
        ref_soln_E_350mm = ref_soln_E[2:end, 1:2]
        data_dates = ref_soln_E_350mm[:, 1]
        data_evaprate = ref_soln_E_350mm[:, 2]

        # Compare our data to Figure 8b of Lehmann, Assouline, Or  (Phys Rev E 77, 2008)
        plt_fig8b = Plots.plot(size = (600, 400))
        Plots.plot!(
            plt_fig8b,
            data_dates,
            data_evaprate,
            xlabel = "Day",
            ylabel = "Evaporation rate (mm/d)",
            label = "Data",
            linewidth = 3,
            margins = 6Plots.mm,
            title = "Soil evaporation rate over time",
            xlim = [minimum(data_dates), maximum(data_dates)],
        )
        Plots.plot!(
            plt_fig8b,
            sol.t ./ 3600 ./ 24,
            evap .* (1000 * 3600 * 24),
            label = "Model",
            color = :black,
            linewidth = 3,
        )
        Plots.plot(plt_fig8b)
        Plots.savefig(joinpath(savepath, "evaporation_lehmann2008_fig8b.png"))
    end
end

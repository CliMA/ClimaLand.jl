import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Soil.Biogeochemistry: MicrobeProduction
import ClimaLand.Simulations: LandSimulation, solve!
using Dates

import ClimaLand.Parameters as LP
import ClimaParams

# Define simulation times
t0 = Float64(0)
tf = Float64(10000)
dt = Float64(10)

for (FT, tf) in ((Float32, 2 * dt), (Float64, tf))
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)

    # Make soil model
    ν = FT(0.556)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    θ_r = FT(0.1)
    ν_ss_om = FT(0.0)
    ν_ss_quartz = FT(1.0)
    ν_ss_gravel = FT(0.0)
    soil_parameters = Soil.EnergyHydrologyParameters(
        toml_dict;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n),
        K_sat,
        S_s,
        θ_r,
    )

    zmax = FT(0)
    zmin = FT(-1)
    nelems = 20

    domain = Column(; zlim = (zmin, zmax), nelements = nelems)

    top_flux_bc_w = Soil.WaterFluxBC((p, t) -> -0.00001)
    bot_flux_bc_w = Soil.FreeDrainage()

    top_flux_bc_h = Soil.HeatFluxBC((p, t) -> 0.0)
    bot_flux_bc_h = Soil.HeatFluxBC((p, t) -> 0.0)


    sources = (PhaseChange{FT}(),)
    boundary_conditions = (;
        top = WaterHeatBC(; water = top_flux_bc_w, heat = bot_flux_bc_h),
        bottom = WaterHeatBC(; water = bot_flux_bc_w, heat = bot_flux_bc_h),
    )
    soil = EnergyHydrology{FT}(;
        boundary_conditions,
        sources,
        domain,
        parameters = soil_parameters,
    )

    # Make a PrescribedAtmosphere - we only care about atmos_p though
    precipitation_function = (t) -> 1.0
    snow_precip = (t) -> 1.0
    atmos_T = (t) -> 1.0
    atmos_u = (t) -> 1.0
    atmos_q = (t) -> 1.0
    atmos_p = (t) -> 100000.0
    UTC_DATETIME = Dates.now()
    atmos_h = FT(30)
    atmos_co2 = (t) -> 1.0

    atmos = ClimaLand.PrescribedAtmosphere(
        TimeVaryingInput(precipitation_function),
        TimeVaryingInput(snow_precip),
        TimeVaryingInput(atmos_T),
        TimeVaryingInput(atmos_u),
        TimeVaryingInput(atmos_q),
        TimeVaryingInput(atmos_p),
        UTC_DATETIME,
        atmos_h,
        earth_param_set;
        c_co2 = TimeVaryingInput(atmos_co2),
    )

    # Make biogeochemistry model
    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

    co2_parameters = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
    C = FT(100)

    co2_top_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> 0.0)
    co2_bot_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> 0.0)
    co2_boundary_conditions = (; top = co2_top_bc, bottom = co2_bot_bc)
    drivers = Soil.Biogeochemistry.SoilDrivers(
        PrognosticMet(soil_parameters),
        Csom,
        atmos,
    )
    soilco2 = SoilCO2Model{FT}(
        domain,
        drivers;
        boundary_conditions = co2_boundary_conditions,
        parameters = co2_parameters,
    )

    model = LandSoilBiogeochemistry{FT}(soil, soilco2)

    function init_soil!(Y, z, params)
        ν = params.ν
        FT = eltype(Y.soil.ϑ_l)
        Y.soil.ϑ_l .= FT(0.33)
        Y.soil.θ_i .= FT(0.1)
        T = FT(279.85)
        ρc_s = Soil.volumetric_heat_capacity(
            FT(0.33),
            FT(0.1),
            params.ρc_ds,
            params.earth_param_set,
        )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                FT(0.0),
                ρc_s,
                T,
                params.earth_param_set,
            )
    end

    function init_co2!(Y, z)
        function CO2_profile(z::FT) where {FT}
            C = FT(0.0)
            return FT(C)
        end
        Y.soilco2.C .= CO2_profile.(z)
    end

    function set_ic!(Y, p, t0, model)
        z = model.soil.domain.fields.z
        init_co2!(Y, z)
        init_soil!(Y, z, model.soil.parameters)
    end

    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)

    saveat = collect(t0:FT(10 * dt):tf)
    sv = (;
        t = Array{Float64}(undef, length(saveat)),
        saveval = Array{NamedTuple}(undef, length(saveat)),
    )
    saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
    updateat = deepcopy(saveat)

    simulation = LandSimulation(
        t0,
        tf,
        dt,
        model;
        diagnostics = [],
        timestepper = ode_algo,
        set_ic!,
        user_callbacks = (saving_cb,),
        updateat = updateat,
        solver_kwargs = (; saveat = deepcopy(saveat)),
    )
    sol = ClimaLand.Simulations.solve!(simulation)

    # Check that simulation still has correct float type
    @assert eltype(sol.u[end].soil) == FT
    @assert eltype(sol.u[end].soilco2) == FT

    # Animation
    # You will need to ]add GlMakie to your base Julia Project.toml
    #=
    if FT == Float64
        using GLMakie
        fig = Figure(resolution = (1100, 700))

        depth = parent(coords.subsurface.z)[:]
        t = Observable(1)

        Mcolor = @lift(parent(sol.u[$t].soil.ϑ_l)[:])
        Tcolor = @lift(parent(saved_values.saveval[$t].soil.T)[:])
        Ccolor = @lift(parent(sol.u[$t].soilco2.C)[:])
        Icolor = @lift(parent(sol.u[$t].soil.θ_i)[:])

        figtitle = @lift("Time: 10s * " * string($t))

        ax_M = Axis(
            fig[1, 1],
            title = "Moisture",
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            xtrimspine = true,
            ytrimspine = true,
            ylabel = "depth (m)",
            xticksvisible = false,
            xticklabelsvisible = false,
        )

        ax_T = Axis(
            fig[1, 3],
            title = "Temperature",
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            xtrimspine = true,
            ytrimspine = true,
            yticksvisible = false,
            yticklabelsvisible = false,
            xticksvisible = false,
            xticklabelsvisible = false,
        )

        ax_C = Axis(
            fig[1, 5],
            title = "CO2",
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            xtrimspine = true,
            ytrimspine = true,
            yticksvisible = false,
            yticklabelsvisible = false,
            xticksvisible = false,
            xticklabelsvisible = false,
        )

        ax_I = Axis(
            fig[1, 7],
            title = "Ice",
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            xtrimspine = true,
            ytrimspine = true,
            yticksvisible = false,
            yticklabelsvisible = false,
            xticksvisible = false,
            xticklabelsvisible = false,
        )

        lines!(
            ax_M,
            ones(nelems),
            depth,
            color = Mcolor,
            linewidth = 140,
            colormap = :deep,
            colorrange = [0.31, 0.6],
        )
        lines!(
            ax_T,
            ones(nelems),
            depth,
            color = Tcolor,
            linewidth = 140,
            colormap = Reverse(:autumn1),
            colorrange = [274, 280],
        )
        lines!(
            ax_C,
            ones(nelems),
            depth,
            color = Ccolor,
            linewidth = 140,
            colormap = Reverse(:grays),
            colorrange = [0.0, 0.006],
        )
        lines!(
            ax_I,
            ones(nelems),
            depth,
            color = Icolor,
            linewidth = 140,
            colormap = :Blues,
            colorrange = [0.0, 0.1],
        )

        Colorbar(
            fig[1, 2],
            limits = (0.31, 0.35),
            label = "Soil Moisture",
            colormap = :deep,
        )
        Colorbar(
            fig[1, 4],
            limits = (274, 280),
            label = "Soil Temperature",
            colormap = Reverse(:autumn1),
        )
        Colorbar(
            fig[1, 6],
            limits = (0.0, 0.006),
            label = "Soil CO2",
            colormap = Reverse(:grays),
        )
        Colorbar(fig[1, 8], limits = (0.0, 0.1), label = "Ice", colormap = :Blues)

        xlims!(ax_M, (0, 2))
        xlims!(ax_T, (0, 2))
        xlims!(ax_C, (0, 2))
        xlims!(ax_I, (0, 2))

        hidespines!(ax_M, :b)
        hidespines!(ax_T, :b, :l)
        hidespines!(ax_C, :b, :l)
        hidespines!(ax_I, :b, :l)

        Label(fig[0, :], text = figtitle, fontsize = 30)

        framerate = 100
        timestamps = range(1, 1000, step = 1)

        filename = "./soil_biogeochem_animation.mp4"
        record(fig, filename, timestamps; framerate = framerate) do dayt
            t[] = dayt
        end
    end
    =#
end

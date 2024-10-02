import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
using CairoMakie
using Statistics
using Dates

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Snow
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams
import ClimaUtilities.OutputPathGenerator: generate_output_path

global climaland_dir = pkgdir(ClimaLand)
global site_ID = "US-MOz"

for float_type in (Float32, Float64)
    # Make these global so we can use them in other ozark files
    global FT = float_type
    global earth_param_set = LP.LandParameters(FT)

    # Create model, set initial conditions, and setup most simulation parameters
    include(
        joinpath(
            climaland_dir,
            "experiments/integrated/performance/conservation/ozark_conservation_setup.jl",
        ),
    )

    # Use smaller `tf` for Float32 simulation
    tf = (FT == Float64) ? t0 + 3600 * 24 * 10 : t0 + 2 * dt
    saveat = Array(t0:dt:tf)
    sv = (;
        t = Array{Float64}(undef, length(saveat)),
        saveval = Array{NamedTuple}(undef, length(saveat)),
    )
    saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)

    updateat = deepcopy(saveat)
    drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
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
    cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
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
    @assert eltype(sol.u[end].soilco2) == FT
    @assert eltype(sol.u[end].canopy) == FT

    # Plotting for Float64 simulation
    if FT == Float64
        # Check that the driver variables stored in `p` are
        # updated correctly - just check one from radiation,
        # and one from atmos.
        cache_θs = [parent(sv.saveval[k].drivers.θs)[1] for k in 1:length(sv.t)]
        cache_Tair =
            [parent(sv.saveval[k].drivers.T)[1] for k in 1:length(sv.t)]
        @assert mean(
            abs.(radiation.θs.(sv.t, radiation.start_date) .- cache_θs),
        ) < eps(FT)
        T_mutable = Vector{FT}(undef, 1)
        atmos_T = map(sv.t) do time
            ClimaLand.evaluate!(T_mutable, atmos.T, time)
            return T_mutable[]
        end |> collect

        @assert mean(abs.(atmos_T .- cache_Tair)) < eps(FT)

        daily = sol.t ./ 3600 ./ 24
        savedir = generate_output_path(
            "experiments/integrated/performance/conservation",
        )

        # Assess conservation
        _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
        _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
        fig = CairoMakie.Figure(size = (1600, 1200), fontsize = 26)
        ax1 = CairoMakie.Axis(
            fig[2, 1],
            ylabel = "ΔEnergy (J/A)",
            xlabel = "Days",
        )
        ΔE_expected =
            cumsum(
                -1 .* [
                    parent(
                        sv.saveval[k].atmos_energy_flux .-
                        sv.saveval[k].soil.bottom_bc.heat,
                    )[end] for k in 1:1:(length(sv.t) - 1)
                ],
            ) * (sv.t[2] - sv.t[1])
        E_measured = [
            sum(sol.u[k].soil.ρe_int) +
            parent(sol.u[k].snow.U)[end] +
            parent(
                sol.u[k].canopy.energy.T .* ac_canopy .*
                max.(
                    sv.saveval[k].canopy.hydraulics.area_index.leaf .+
                    sv.saveval[k].canopy.hydraulics.area_index.stem,
                    eps(FT),
                ),
            )[end] for k in 1:1:length(sv.t)
        ]
        ΔW_expected =
            cumsum(
                -1 .* [
                    parent(
                        sv.saveval[k].atmos_water_flux .-
                        sv.saveval[k].soil.bottom_bc.water .+
                        sv.saveval[k].soil.R_s,
                    )[end] for k in 1:1:(length(sv.t) - 1)
                ],
            ) * (sv.t[2] - sv.t[1])

        W_measured = [
            sum(sol.u[k].soil.ϑ_l) +
            sum(sol.u[k].soil.θ_i) * _ρ_i / _ρ_l +
            sum(parent(sol.u[k].canopy.hydraulics.ϑ_l) .* [h_stem, h_leaf]) +
            parent(sol.u[k].snow.S)[end] for k in 1:1:length(sv.t)
        ]
        CairoMakie.lines!(
            ax1,
            daily[2:end],
            E_measured[2:end] .- E_measured[1],
            label = "Simulated",
        )
        CairoMakie.lines!(ax1, daily[2:end], ΔE_expected, label = "Expected")
        CairoMakie.axislegend(ax1, position = :rt)

        # Temp
        ax4 = CairoMakie.Axis(fig[1, 1], ylabel = "ΔWater (m)")
        CairoMakie.hidexdecorations!(ax4, ticks = false)
        CairoMakie.lines!(
            ax4,
            daily[2:end],
            W_measured[2:end] .- W_measured[1],
            label = "Simulated",
        )

        CairoMakie.lines!(ax4, daily[2:end], ΔW_expected, label = "Expected")
        CairoMakie.axislegend(ax4, position = :rt)


        ax3 = CairoMakie.Axis(
            fig[2, 2],
            ylabel = "ΔE/E",
            xlabel = "Days",
            yscale = log10,
        )
        CairoMakie.lines!(
            ax3,
            daily[2:end],
            abs.(E_measured[2:end] .- E_measured[1] .- ΔE_expected) ./
            mean(E_measured),
        )

        ax2 = CairoMakie.Axis(fig[1, 2], ylabel = "ΔW/W", yscale = log10)
        CairoMakie.hidexdecorations!(ax2, ticks = false)
        CairoMakie.lines!(
            ax2,
            daily[2:end],
            abs.(W_measured[2:end] .- W_measured[1] .- ΔW_expected) ./
            mean(W_measured),
        )

        CairoMakie.save(joinpath(savedir, "results_conservation.png"), fig)
    end
end

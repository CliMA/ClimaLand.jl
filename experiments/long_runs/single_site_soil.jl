import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date
using Statistics
import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaCore
using ClimaLand
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

using Dates

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis
const FT = Float64;
# If you want to do a very long run locally, you can enter `export
# LONGER_RUN=""` in the terminal and run this script. If you want to do a very
# long run on Buildkite manually, then make a new build and pass `LONGER_RUN=""`
# as an environment variable. In both cases, the value of `LONGER_RUN` does not
# matter.
const LONGER_RUN = haskey(ENV, "LONGER_RUN") ? true : false
# If you want to do run the simulation with uncalibrated parameters, type
# `export UNCALIBRATED=""` in the terminal and run this script, or
# pass `UNCALIBRATED=""` as an environment variable on buildkite.
const UNCALIBRATED = haskey(ENV, "UNCALIBRATED") ? true : false
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "soil_main_no_water_fluxes_long_eastern_kersten"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
) where {FT}
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data - high resolution
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
        use_lowres_forcing = true,
    )
    forcing = (; atmos, radiation)
    # Construct the land model with all default components except for snow
    land = ClimaLand.Soil.EnergyHydrology{FT}(domain,
                                              forcing,
                                              toml_dict,
                                              )
    return land
end

# If not LONGER_RUN, run for 2 years; note that the forcing from 2008 is repeated.
# If LONGER run, run for 19 years, with the correct forcing each year.
# Note that since the Northern hemisphere's winter season is defined as DJF,
# we simulate from and until the beginning of
# March so that a full season is included in seasonal metrics.
start_date = DateTime("2000-09-01")
stop_date = DateTime("2040-09-01")
Δt = 1800.0
longlat = FT.((140.0, 58.0))#FT.((-128.3, 69.1))
zlim = FT.((-15, 0))
nelements = 15
dz_tuple = FT.((3, 0.05))
domain = ClimaLand.Domains.Column(; zlim, longlat, nelements, dz_tuple)

if UNCALIBRATED
    override_params_path = "toml/uncalibrated_parameters.toml"
    toml_dict = LP.create_toml_dict(FT, override_files = [override_params_path])
else
    toml_dict = LP.create_toml_dict(FT)
end

model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    outdir;
    reduction_period = :monthly,
    output_vars = [
        "tsoil",
        "swc",
        "si",
        "sr",
        "ssr",
        "tair",
        "precip",
    ],
)
set_ic! =
    ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
        model,
    )
simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    model;
    outdir,
    diagnostics,
    set_ic!,
)

@info "Run: Global Soil-Canopy-Snow Model"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation)

#LandSimVis.make_timeseries(simulation; savedir = root_path)

times = collect(keys(diagnostics[1].output_writer.dict["tsoil_1M_average"]))
fig = CairoMakie.Figure(size = (800, 1200))
mean_T_air = mean([
            parent(diagnostics[1].output_writer.dict["tair_1M_average"][t])[1]
            for t in times
        ])
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Tsoil -mean(T_air)")
for i in 1:1:nelements
    lines!(
        ax1,
        [
            parent(diagnostics[1].output_writer.dict["tsoil_1M_average"][t])[i]
            for t in times
        ] .- mean_T_air,
        label = "$i",
    )
end

ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Ice")
for i in 1:1:nelements
    lines!(
        ax2,
        [
            parent(diagnostics[1].output_writer.dict["si_1M_average"][t])[i] for
            t in times
        ],
        label = "$i",
    )
end
fig[2, 2] = Legend(fig, ax2)

ax3 =
    CairoMakie.Axis(fig[3, 1], xlabel = "Time", ylabel = "Effective Saturation")
for i in 1:1:nelements
    lines!(
        ax3,
        [
            (
                parent(
                    diagnostics[1].output_writer.dict["si_1M_average"][t],
                )[i] .+ parent(
                    diagnostics[1].output_writer.dict["swc_1M_average"][t],
                )[i] - parent(model.parameters.θ_r)[i]
            ) / (
                parent(model.parameters.ν)[i] -
                parent(model.parameters.θ_r)[i]
            ) for t in times
        ],
        label = "$i",
    )
end
CairoMakie.save(joinpath(root_path, "soil.png"), fig)
fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "SSR")
lines!(
    ax1,
    [
        parent(diagnostics[1].output_writer.dict["ssr_1M_average"][t])[1] for
        t in times
    ],
)
ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "SR")
lines!(
    ax2,
    [
        parent(diagnostics[1].output_writer.dict["sr_1M_average"][t])[1] for
        t in times
    ],
)
CairoMakie.save(joinpath(root_path, "runoff.png"), fig)

fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Tsoil -mean(T_air)")
for i in 1:3:nelements
    lines!(
        ax1,
        [
            parent(diagnostics[1].output_writer.dict["tsoil_1M_average"][t])[i]
            for t in times
        ] .- mean_T_air,
        label = "$i",
    )
end
ylims!(ax1, -4, 4)

ax2 =
    CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Effective Saturation")
for i in 1:3:nelements
    lines!(
        ax2,
        [
            (
                parent(
                    diagnostics[1].output_writer.dict["si_1M_average"][t],
                )[i] .+ parent(
                    diagnostics[1].output_writer.dict["swc_1M_average"][t],
                )[i] - parent(model.parameters.θ_r)[i]
            ) / (
                parent(model.parameters.ν)[i] -
                parent(model.parameters.θ_r)[i]
            ) for t in times
        ],
        label = "$i",
    )
end
fig[2, 2] = Legend(fig, ax2)
CairoMakie.save(joinpath(root_path, "zoom_in.png"), fig)
avg = zeros(nelements)
for i in 1:nelements
    @show i
    avg[i] = mean([parent(diagnostics[1].output_writer.dict["tsoil_1M_average"][t])[i]
            for t in times[length(times)-11-12:end]
        ])
end

fig = CairoMakie.Figure(size = (800, 800))
ax = Axis(fig[1,1])
lines!(ax, avg, parent(model.domain.fields.z)[:])
lines!(ax, zeros(nelements) .+ mean_T_air, parent(model.domain.fields.z)[:])
lines!(ax, zeros(nelements) .+ 273.15, parent(model.domain.fields.z)[:])
CairoMakie.save(joinpath(root_path, "mean_profile.png"), fig)

@show sum(  [
               parent(diagnostics[1].output_writer.dict["ssr_1M_average"][t])[1] for
               t in times
           ])
@show sum(  [
               parent(diagnostics[1].output_writer.dict["sr_1M_average"][t])[1] for
               t in times
           ])
@show sum(  [
               parent(diagnostics[1].output_writer.dict["precip_1M_average"][t])[1] for
               t in times
           ])
#=
Main
i = 1
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.4967704094152663
i = 2
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.49684954577961593
i = 3
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.49697736637343426
i = 4
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.49712305337673995
i = 5
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.49746657837329744
i = 6
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.5066077847641213
i = 7
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.8054803260481647
i = 8
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -1.0381149168846031
i = 9
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -1.2389596320764082
i = 10
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -1.4039122076472292
i = 11
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -1.4101910904067319
i = 12
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -1.1147485285221352
i = 13
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.5189086614952271
i = 14
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = 0.07154352007040106
i = 15
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = 0.414779930781491
sum([(parent(((diagnostics[1]).output_writer.dict["ssr_1d_average"])[t]))[1] for t = times]) = 3.046659720842384e-7
sum([(parent(((diagnostics[1]).output_writer.dict["sr_1d_average"])[t]))[1] for t = times]) = 3.0314972990072472e-5
sum([(parent(((diagnostics[1]).output_writer.dict["precip_1d_average"])[t]))[1] for t = times]) = -0.22537400948916383
-0.22537400948916383

# Relax
i = 1
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.07411125273485988
i = 2
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.07487787587205648
i = 3
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.076187951685527
i = 4
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.07758379553405721
i = 5
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.07912071170514602
i = 6
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.08132713043172715
i = 7
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.08575634279868251
i = 8
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.09566780366748598
i = 9
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.12415024612517515
i = 10
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.19931272923213905
i = 11
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.2293658945122783
i = 12
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = -0.10674319535939904
i = 13
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = 0.1746314270944652
i = 14
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = 0.4957676176258909
i = 15
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = 0.7521587345051377
mean([(parent(((diagnostics[1]).output_writer.dict["tsoil_1d_average"])[t]))[i] for t = times[(end - 364) - 365:end]] .- mean_T_air) = 0.7521587345051377
sum([(parent(((diagnostics[1]).output_writer.dict["ssr_1d_average"])[t]))[1] for t = times]) = 8.222772717963096e-5
sum([(parent(((diagnostics[1]).output_writer.dict["sr_1d_average"])[t]))[1] for t = times]) = 6.428602704200076e-6
sum([(parent(((diagnostics[1]).output_writer.dict["precip_1d_average"])[t]))[1] for t = times]) = -0.22537400948916383
=#
# Δs compared to mean air temp are much smaller for all layers but bigger at the surface. SR lower, SSR a lot higher (7e-9) at end of simulation. Probably from the bottom layer

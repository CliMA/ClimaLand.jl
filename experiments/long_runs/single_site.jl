import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
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
root_path = "snowy_land_pmodel_longrun_$(device_suffix)_main_short"
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

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    prognostic_land_components = (:canopy, :snow, :soil)

    # Construct the P model manually since it is not a default
    photosynthesis = PModel{FT}(domain, toml_dict)
    conductance = PModelConductance{FT}(toml_dict)
    # Use the soil moisture stress function based on soil moisture only
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(domain, toml_dict)
    biomass =
        ClimaLand.Canopy.PrescribedBiomassModel{FT}(domain, LAI, toml_dict)
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        biomass,
    )

    # Snow model setup
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
    horz_degree_res = FT(1)
    scf = Snow.WuWuSnowCoverFractionModel(toml_dict, horz_degree_res)
    surf_temp = Snow.EquilibriumGradientTemperatureModel{FT}()
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components,
        α_snow,
        scf,
        #        surf_temp,
    )

    # Construct the land model with all default components except for snow
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        snow,
        canopy,
    )
    return land
end

# If not LONGER_RUN, run for 2 years; note that the forcing from 2008 is repeated.
# If LONGER run, run for 19 years, with the correct forcing each year.
# Note that since the Northern hemisphere's winter season is defined as DJF,
# we simulate from and until the beginning of
# March so that a full season is included in seasonal metrics.
start_date = DateTime("2008-01-01")
stop_date = DateTime("2010-01-01")
Δt = 450.0
longlat = FT.((-74.4, 47.7))
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
    reduction_period = :daily,
    output_vars = [
        "tsoil",
        "swc",
        "si",
        "snowtb",
        "snowtsfc",
        "snowc",
        "swe",
        "snd",
        "sr",
        "ssr",
        "infc",
    ],
)

simulation =
    LandSimulation(start_date, stop_date, Δt, model; outdir, diagnostics)

@info "Run: Global Soil-Canopy-Snow Model"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_timeseries(simulation; savedir = root_path)

times = collect(keys(diagnostics[1].output_writer.dict["tsoil_1d_average"]))
fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Tsoil")
for i in 1:2:15
    lines!(
        ax1,
        [
            parent(diagnostics[1].output_writer.dict["tsoil_1d_average"][t])[i]
            for t in times
        ],
        label = "$i",
    )
end

ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Ice")
for i in 1:2:15
    lines!(
        ax2,
        [
            parent(diagnostics[1].output_writer.dict["si_1d_average"][t])[i] for
            t in times
        ],
        label = "$i",
    )
end
fig[2, 2] = Legend(fig, ax2)

ax3 = CairoMakie.Axis(fig[3, 1], xlabel = "Time", ylabel = "Liquid Water")
for i in 1:2:15
    lines!(
        ax3,
        [
            parent(diagnostics[1].output_writer.dict["swc_1d_average"][t])[i]
            for t in times
        ],
        label = "$i",
    )
end
CairoMakie.save(joinpath(root_path, "soil.png"), fig)

fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Snow depth")
lines!(
    ax1,
    [
        parent(diagnostics[1].output_writer.dict["snd_1d_average"][t])[1] for
        t in times
    ],
)
ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Snow cover")
lines!(
    ax2,
    [
        parent(diagnostics[1].output_writer.dict["snowc_1d_average"][t])[1] for
        t in times
    ],
)
ax3 = CairoMakie.Axis(fig[3, 1], xlabel = "Time", ylabel = "Temps")
lines!(
    ax3,
    [
        parent(diagnostics[1].output_writer.dict["snowtb_1d_average"][t])[1] for
        t in times
    ],
    label = "Bulk Snow",
)
lines!(
    ax3,
    [
        parent(diagnostics[1].output_writer.dict["snowtsfc_1d_average"][t])[1]
        for t in times
    ],
    label = "Surface Snow",
)
lines!(
    ax3,
    [
        parent(diagnostics[1].output_writer.dict["tsoil_1d_average"][t])[15] for
        t in times
    ],
    label = "Surface Soil",
)
fig[3, 2] = Legend(fig, ax3)
CairoMakie.save(joinpath(root_path, "snow.png"), fig)

fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Inf c")
lines!(
    ax1,
    [
        parent(diagnostics[1].output_writer.dict["infc_1d_average"][t])[1] for
        t in times
    ],
)
ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "SR")
lines!(
    ax2,
    [
        parent(diagnostics[1].output_writer.dict["sr_1d_average"][t])[1] for
        t in times
    ],
)
CairoMakie.save(joinpath(root_path, "runoff.png"), fig)
@show sum([
    parent(diagnostics[1].output_writer.dict["sr_1d_average"][t])[1] for
    t in times
])
@show sum([
    parent(diagnostics[1].output_writer.dict["ssr_1d_average"][t])[1] for
    t in times
])

# # Global spin-up run of land model

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

const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "snowy_land_pmodel_spinup_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(
    ::Type{FT},
    Δt,
    domain,
    toml_dict,
) where {FT}
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data - high resolution
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        DateTime("1979"),
        DateTime("1998"),
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
    ); # Will repeat these 20 years
    forcing = (; atmos, radiation);
    LAI = prescribed_climatological_lai_modis(surface_space);

    # Construct the land model with default components
    prognostic_land_components = (:canopy, :snow, :soil,);
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
    );
    return land
end

Δt = 900.0
domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
toml_dict = LP.create_toml_dict(FT)

start_date = DateTime("2000-03-01")
stop_date = DateTime("2099-03-01")

model = setup_model(FT, Δt, domain, toml_dict);
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    outdir;
    reduction_period=:monthly,
    output_vars = [
        "tsoil",
        "sie",
        "swc",
        "si",
        "lwp",
        "tair",
        "swe",
    ],
)
set_ic! =
    ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
        model,
    )
simulation = LandSimulation(start_date, stop_date, Δt, model; outdir, set_ic!, diagnostics)
@info "Run: Global Soil-Canopy-Snow Model"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation)




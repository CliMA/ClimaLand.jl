# Soil-Canopy-Snow-SoilCO2 column at IIT Kanpur (26.51 N, 80.23 E), in the
# Indo-Gangetic plain of Uttar Pradesh, India.
#
# This run needs high resolution ERA5 forcing covering the site. That artifact
# (`forty_yrs_era5_land_forcing_data`) is not downloadable; supply it yourself by
# placing files named `era5_<year>_1.0x1.0.nc` in a directory and pointing
# `~/.julia/artifacts/Overrides.toml` at it. The low resolution alternative is
# unusable here: its 8 degree grid snaps this site to a Tibetan Plateau cell whose
# annual mean temperature is -9.7 C.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
import ClimaUtilities
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

using Dates

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
root_path = "iitk_site"
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
    surface_space = domain.space.surface
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
    )
    forcing = (; atmos, radiation)

    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    # Every component is left at its default. This matters for the canopy in
    # particular: the default ties its piecewise moisture stress thresholds to the
    # soil model's own retention parameters, and building the two separately leaves
    # them inconsistent, which LandModel rejects.
    return LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    )
end

start_date = DateTime("2010-01-01")
# Ending inside 2010 keeps the forcing to a single year file: find_era5_year_paths
# requires a file for every year in year(start_date):year(stop_date)
stop_date = DateTime("2010-12-31")
Δt = 450.0
longlat = FT.((80.23, 26.51))
# The 15 m extent matches the default spun-up initial conditions artifact
zlim = FT.((-15, 0))
nelements = 15
dz_tuple = FT.((3, 0.05))
domain = ClimaLand.Domains.Column(; zlim, longlat, nelements, dz_tuple);
toml_dict = LP.create_toml_dict(FT)

model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict);
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    outdir;
    reduction_period = :daily,
    output_vars = [
        "swc",
        "tsoil",
        "lhf",
        "shf",
        "et",
        "trans",
        "gpp",
        "lai",
        "precip",
        "tair",
        "swe",
        "infil",
        "nee",
    ],
);
simulation =
    LandSimulation(start_date, stop_date, Δt, model; outdir, diagnostics);

@info "Run: Soil-Canopy-Snow-SoilCO2 Model at IIT Kanpur"
@info "Location: $(longlat[2]) N, $(longlat[1]) E"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
solve!(simulation);

LandSimVis.make_timeseries(simulation; savedir = root_path)

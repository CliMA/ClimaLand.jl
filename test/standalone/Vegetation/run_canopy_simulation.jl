ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
required_pkgs = ["CairoMakie", "ClimaAnalysis", "GeoMakie", "Printf", "StatsBase", "ClimaComms", "ClimaParams", "Dates", "ClimaUtilities", "Interpolations", "ClimaDiagnostics", "ClimaLand", "TOML"]
import Pkg
Pkg.Registry.update()
Pkg.add(required_pkgs)

# Create a directory structure in the artifacts directory to point to your local data
using TOML

# The artifact hash that ClimaLand is looking for
artifact_hash = "f269a0b057b9f438b4caafdef17da73746310787"
artifacts_dir = joinpath(homedir(), ".julia", "artifacts")
artifact_path = joinpath(artifacts_dir, artifact_hash)

# Create the artifacts directory if it doesn't exist
mkpath(artifacts_dir)

# Remove old symlink or directory if it exists
if ispath(artifact_path)
    rm(artifact_path; recursive=true, force=true)
    println("Removed old artifact path")
end

# Create the artifact directory
mkpath(artifact_path)
println("Created artifact directory: $artifact_path")

# Now create symlinks for the files with the expected naming convention
# ClimaLand expects era5_YYYY_1.0x1.0.nc but you have era5_YYYY_0.25x0.25.nc
local_data_path = "/net/sampo/data1/era5/high_res_2008_era5_land_forcing_data"
src_file = joinpath(local_data_path, "era5_2008_0.25x0.25.nc")
dst_file = joinpath(artifact_path, "era5_2008_1.0x1.0.nc")
if isfile(src_file)
    symlink(src_file, dst_file)
    println("Created file symlink: era5_2008_1.0x1.0.nc -> $(src_file)")
end

# Do the same for LAI file if needed
src_lai = joinpath(local_data_path, "era5_2008_0.25x0.25_lai.nc")
dst_lai = joinpath(artifact_path, "era5_2008_1.0x1.0_lai.nc")
if isfile(src_lai)
    symlink(src_lai, dst_lai)
    println("Created LAI symlink: era5_2008_1.0x1.0_lai.nc -> $(src_lai)")
end

import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using Dates
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import Interpolations
using ClimaDiagnostics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

pkgversion(ClimaLand)

const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "land_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir);
toml_dict = CP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict);

Δt = 450.0
start_date = DateTime(2008)
stop_date = DateTime(2008, 12, 31);  # Changed from DateTime(2009)

nelements = (101, 15)
domain = ClimaLand.Domains.global_domain(FT; context, nelements);

forcing = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = false,
    max_wind_speed = 25.0,
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
    context = context,
);

# Use LAI from ERA5 forcing data
lai_file = joinpath(artifact_path, "era5_2008_1.0x1.0_lai.nc")

LAI = ClimaLand.prescribed_lai_modis(
    domain.space.surface,  # Space comes first
    start_date,             # Then start date
    stop_date;              # Then stop date (as keyword or positional)
    modis_lai_ncdata_path = [lai_file],  # Pass file path as keyword argument
    time_interpolation_method = LinearInterpolation(),
    context = context,
);

model = ClimaLand.LandModel{FT}(forcing, LAI, earth_param_set, domain, Δt);

subsurface_space = domain.space.subsurface
nc_writer =
   ClimaDiagnostics.Writers.NetCDFWriter(subsurface_space, outdir; start_date);
diagnostics = ClimaLand.default_diagnostics(
   model,
   start_date;
   output_writer = nc_writer,
   average_period = :daily,
);

simulation = ClimaLand.Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt,
    model;
    outdir,
    diagnostics,
    user_callbacks = (ClimaLand.ReportCallback(10000),)
);

ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation;date = stop_date, savedir = root_path)
LandSimVis.make_leaderboard_plots(simulation, savedir = root_path)
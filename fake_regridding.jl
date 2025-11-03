const FT = Float64
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput,
    AbstractTimeVaryingInput,
    LinearInterpolation,
    PeriodicCalendar
import Interpolations: Constant
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
import ClimaUtilities.TimeManager: ITime, date
using Thermodynamics
using ClimaCore
using Dates
using Insolation
using SurfaceFluxes
import SurfaceFluxes.Parameters as SFP
using StaticArrays
import ClimaLand
import Interpolations
import ClimaUtilities

import ClimaLand.Parameters as LP

import ClimaCore.Geometry
import ClimaCore.Remapping

import ClimaAnalysis
import ClimaAnalysis.Template:
    TemplateVar,
    make_template_var,
    add_attribs,
    add_dim,
    add_time_dim,
    add_lon_dim,
    add_lat_dim,
    add_data,
    ones_data,
    zeros_data,
    one_to_n_data,
    initialize

import GeoMakie
import CairoMakie

import ClimaComms
ClimaComms.@import_required_backends
context = ClimaComms.context()

start_date = DateTime("2008-03-01")
stop_date = DateTime("2008-06-01")

nelements = (101, 15)
domain = ClimaLand.Domains.global_domain(
    FT;
    context,
    nelements,
    mask_threshold = FT(0.99),
)
surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
surface_space = domain.space.surface

era5_2008_data = "/net/sampo/data1/era5/forty_yrs_era5_land_forcing_data/forty_yrs_era5_land_forcing_data_artifact/era5_2008_1.0x1.0.nc"
era5_ncdata_path = [era5_2008_data]

regridder_type = :InterpolationsRegridder
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
interpolation_method = Interpolations.Linear()

T_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "t2m",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )

# TODO: Instead of this step, we can interpolate to an array whose coordinates
# are the same as the coordinate fields
# Then, use this to interpolate instead
# Instead of relying on ClimaCore for field -> array, use ClimaAnalysis for
# array -> array instead
dest = ClimaCore.Fields.zeros(surface_space)
coords = dest |> ClimaCore.Fields.coordinate_field

longs = coords.long |> ClimaCore.Fields.field2array
lats = coords.lat |> ClimaCore.Fields.field2array

longs = unique!(sort(longs))
lats = unique!(sort(lats))

t2m_var_from_nc = ClimaAnalysis.OutputVar(era5_2008_data, "t2m")
t2m_var_from_nc = ClimaAnalysis.select(t2m_var_from_nc; by = ClimaAnalysis.MatchValue(), time = DateTime("2008-03-01"))

# Like going from netcdf file to climacore field
fake_cc = ClimaAnalysis.resampled_as(t2m_var_from_nc, lat = collect(lats), lon = collect(longs))
# Like going from climacore field to netcdf file
longpts = range(-180.0, 180.0, 404)
latpts = range(-90.0, 90.0, 202)
fake_cc = ClimaAnalysis.resampled_as(fake_cc, lon = longpts, lat = latpts)


# Resample t2m_var_from_nc to match the grid of fake_cc
t2m_var_from_nc = ClimaAnalysis.resampled_as(t2m_var_from_nc, fake_cc)

# Make a plot here
fig = CairoMakie.Figure()
ClimaAnalysis.Visualize.plot_bias_on_globe!(fig, fake_cc, t2m_var_from_nc, mask = ClimaAnalysis.apply_oceanmask)

CairoMakie.save("bias_latest_climautilities_fake_cc_all_linear.png", fig)

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
interpolation_method = Interpolations.Constant()

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
fill!(ClimaCore.Fields.field_values(dest), NaN)
ClimaUtilities.TimeVaryingInputs.evaluate!(dest, T_atmos, FT(0.0))

longpts = range(-180.0, 179.0, 360)
latpts = range(-90.0, 90.0, 181)

hcoords = [Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]

# lon by lat
interpolated_array = Remapping.interpolate(dest, hcoords, nothing)

# To make it easier to compare, put it in a OutputVar
t2m_var_from_cc = TemplateVar() |> add_attribs(units = "K", long_name = "2 metre temperature") |> add_data(; data = interpolated_array) |>
        add_dim("lon", collect(longpts), units = "degrees_east") |> add_dim("lat", collect(latpts), units = "degrees_north") |> initialize


# Load data from OutputVar
t2m_var_from_nc = ClimaAnalysis.OutputVar(era5_2008_data, "t2m")
t2m_var_from_nc = ClimaAnalysis.select(t2m_var_from_nc; by = ClimaAnalysis.MatchValue(), time = DateTime("2008-03-01"))

t2m_var_from_nc = ClimaAnalysis.shift_longitude(t2m_var_from_nc, -180.0, 180.0)
# t2m_var_from_nc = ClimaAnalysis.resampled_as(t2m_var_from_nc, t2m_var_from_cc)


# Compare the bias
fig = CairoMakie.Figure()
ClimaAnalysis.Visualize.plot_bias_on_globe!(fig, t2m_var_from_cc, t2m_var_from_nc, cmap_extrema = (-1.0, 1.0))

CairoMakie.save("bias_latest_climautilities_no_regridding_constant.png", fig)

close(T_atmos)

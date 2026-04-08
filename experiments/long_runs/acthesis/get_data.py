import earthaccess
from ecmwf.datastores import Client
import numpy as np
import xarray as xr
import cfgrib
from pyhdf.SD import SD, SDC
from datetime import date, timedelta, datetime
import logging, os, calendar
import xesmf as xe
logging.basicConfig(level="INFO")

def get_MODIS_date_from_id(id):
    return date(int(id[1:5]), 1, 1) + timedelta(days = int(id[5:]) - 1)

def get_MODIS_date_from_result(result):
    id = result['umm']['GranuleUR'].split('.')[1]
    return get_MODIS_date_from_id(id)

def get_month_files(date_dict, year, month):
    _, nday = calendar.monthrange(year, month)
    dates = [date(year, month, day) for day in range(1, nday + 1)]
    return [date_dict.get(d, "NO_DATA_THIS_DATE") for d in dates], dates

def get_hdf_files(to_scrape_list, write_dir):
    earthaccess.download(to_scrape_list, local_path = write_dir)
    return {
        get_MODIS_date_from_id(d.split('.')[1]) : write_dir + "/" + d for d in os.listdir(write_dir)
    }

def parse_hdf_file(filename, fields):
    hdf = SD(filename, SDC.READ)
    data = extract_hdf_fields(hdf, fields)
    hdf.end()
    os.remove(filename)
    return data

def extract_hdf_fields(hdf, names):
    return {name : process_hdf_field(hdf.select(name)) for name in names}

def process_hdf_field(ds):
    attrs = ds.attributes()
    data = ds.get()
    if '_FillValue' in attrs:
        nanval = attrs['_FillValue']
        data = np.where(data == nanval, np.nan, data)
    if 'valid_range' in attrs:
        vmin, vmax = attrs['valid_range']
        data = np.where((data < vmin) | (data > vmax), np.nan, data)
    offset = attrs.get('add_offset', 0.0)
    scaling = attrs.get('scale_factor', 1.0)
    return (data - offset)*scaling

def downscale_mod_field_1(field_big, npix):
    # aggregate via mean to 1 degree by 1 degree pixels (0.05 right now)
    ny, nx = field_big.shape
    assert ny % npix == 0
    assert nx % npix == 0
    field_small_stack = field_big.reshape(ny // npix, npix, nx // npix, npix)
    with np.errstate(invalid ='ignore'):
        field_small = np.nanmean(field_small_stack, axis = (1, 3))
    return field_small

def downscale_mod_field_2(field_big, npix):
    #our points are at -90, 88, .... 89 while method 1 of downscaling interpolates to make +89.5, +88.5, ..., -89.5. Same for longitude.
    #we can interpolate better for the insides while using the half-cell to estimate the -90 and + 90 latitudes
    #1 degree is 20 pixels across - shift longitude by 10 pixels and average that way to get the right integer values
    #first dim is latitude, second dim is longitude
    ny, nx = field_big.shape
    assert ny % npix == 0
    assert nx % npix == 0
    assert npix % 2 == 0
    nshift = npix // 2
    field_shift = np.roll(field_big, shift = nshift, axis = -1) #positive shift bumps right-end to front, pushing 179.5-180 to the front to get averaged with -180 to -179.5
    field_big_inside = field_shift[nshift:-nshift, :]
    field_big_top = field_shift[:nshift, :]
    field_big_bottom = field_shift[-nshift:, :]
    field_inside_stack = field_big_inside.reshape((ny // npix) - 1, npix, nx // npix, npix)
    field_top_stack = field_big_top.reshape(1, nshift, nx // npix, npix)
    field_bottom_stack = field_big_bottom.reshape(1, nshift, nx // npix, npix)
    with np.errstate(invalid ='ignore'):
        field_small_inside = np.nanmean(field_inside_stack, axis = (1, 3)) #now -180, -179, ..., 179
        field_small_top = np.nanmean(field_top_stack, axis = (1, 3))
        field_small_bottom = np.nanmean(field_bottom_stack, axis = (1, 3))
    return np.vstack([field_small_top, field_small_inside, field_small_bottom])

def qa_and_scale(field, qa_field, npix):
    # remove variables based on albedo quality flag: for MCD43C3, keep if quality byte is 2 or less 
    field_big = np.where(qa_field < 3, field.astype(float), np.nan)
    return {"mode1": downscale_mod_field_1(field_big, npix), "mode2" : downscale_mod_field_2(field_big, npix)}

def reduce_hdf_fields(data):
    return {
        "bsa" : qa_and_scale(data["Albedo_BSA_shortwave"], data["Albedo_Quality"], 20),
        "wsa" : qa_and_scale(data["Albedo_WSA_shortwave"], data["Albedo_Quality"], 20),
        "scf" : qa_and_scale(data["Percent_Snow"], data["Albedo_Quality"], 20)
    }

def make_blank_data():
    return {
        "bsa" : {"mode1" : np.full((180, 360), np.nan), "mode2" : np.full((181, 360), np.nan)},
        "wsa" : {"mode1" : np.full((180, 360), np.nan), "mode2" : np.full((181, 360), np.nan)},
        "scf" : {"mode1" : np.full((180, 360), np.nan), "mode2" : np.full((181, 360), np.nan)}
    }

def get_modis_files(date_dict, year, month, dump_dir):
    to_scrape, dates = get_month_files(date_dict, year, month)
    #handle the no-data-dates (findall, return those dates separately as another output to pass to the make_netcdf func):
    no_data_dates = [dates[i] for i, x in enumerate(to_scrape) if x == "NO_DATA_THIS_DATE"]
    scrape_list = [x for x in to_scrape if x != "NO_DATA_THIS_DATE"]
    names = ['Albedo_BSA_shortwave', 'Albedo_WSA_shortwave', 'Albedo_Quality', 'Percent_Snow'] #, 'BRDF_Albedo_Uncertainty']
    fnames = get_hdf_files(scrape_list, dump_dir)
    data = {
        d : reduce_hdf_fields(parse_hdf_file(fnames[d], names)) for d in dates if d not in no_data_dates
    }
    #use no_data_dates and add a blank all-nan field into the data at those dates:
    for d in no_data_dates:
        data[d] = make_blank_data()
    return data

def build_ecmwf_request(year, month):
    return {
    "variable": [
        "snow_albedo",
        "snow_cover",
        #"snow_density",
        "snow_depth",
        "snow_depth_water_equivalent",
        "temperature_of_snow_layer",
        "forecast_albedo",
        #"surface_latent_heat_flux",
        #"surface_net_solar_radiation",
        #"surface_sensible_heat_flux",
        #"surface_solar_radiation_downwards",
    ],
    "year": f"{year}",
    "month": f"{month:02d}",
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15",
        "16", "17", "18",
        "19", "20", "21",
        "22", "23", "24",
        "25", "26", "27",
        "28", "29", "30",
        "31"
    ],
    #"daily_statistic": "daily_mean",
    #"time_zone": "utc+00:00",
    #"frequency": "3_hourly",
    "time": [
        "00:00", "03:00", "06:00",
        "09:00", "12:00", "15:00",
        "18:00", "21:00"
    ],
    "grid": [1.0, 1.0],
    "data_format": "grib",
    "download_format": "unarchived"
}

def get_grib_file(client, year, month, output_path): #you'll want to actually save these to /sampo or somewhere
    dataset_id = "reanalysis-era5-land"
    #dataset_id = "derived-era5-land-daily-statistics"
    request = build_ecmwf_request(year, month)
    client.retrieve(dataset_id, request, target = output_path)

def agg_time(field, axis, nsteps):
    ns = field.shape
    if ns[axis] == nsteps:
        with np.errstate(invalid ='ignore'):
            field_daily = np.nanmean(field, axis = axis)
        return field_daily
    else:
        #would arguably want a moveaxis() call here to do this via the axis arg, but in this case, I'm pretty sure its always axis = 0.
        assert ns[0] % nsteps == 0
        field_daily_stack = field.reshape(ns[0] // nsteps, nsteps, *ns[1:])
        with np.errstate(invalid ='ignore'):
            field_daily = np.nanmean(field_daily_stack, axis = 1)
        return field_daily

def parse_grib_file(path):
    datasets = cfgrib.open_datasets(path)
    grib_data = {}
    for ds in datasets:
        for var in ds.data_vars:
            ddims = list(ds[var].dims)
            data = ds[var].values
            times = ds[var]['time'].values.astype('datetime64[D]').astype(date)
            #collapse along time-dimension
            if len(data.shape) > 3 and 'step' in ddims:
                axis = ddims.index('step')
                data_daily = agg_time(data, axis, data.shape[axis])
                ddims.remove('step')
            else:
                axis = ddims.index('time')
                data_daily = agg_time(data, axis, 8)
                times = sorted(list(set(times)))
            #need to shift array since 0-360 longitude in era5, not -180 to 180
            lon_idx = ddims.index('longitude')
            data_daily = roll_lon(data_daily, lon_idx)
            grib_data[var] = {
                t : data_daily[i, ...] for (i, t) in enumerate(times)
            }
        ds.close()
    #do not delete these files this time, get rid of the idx file though
    return grib_data

def roll_lon(field, lon_idx):
    lon_dim = field.shape[lon_idx]
    return np.roll(field, shift = lon_dim//2, axis = lon_idx)

def collect_monthly_data(year, month, data):
    _, nday = calendar.monthrange(year, month)
    dates = [date(year, month, day) for day in range(1, nday + 1)]
    vars = list(data.keys())
    monthdata = {
        d : {
            v : data[v][d] for v in vars
        } for d in dates
    }
    return monthdata

def get_era5land_data(client, year, month, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    fname = f"{output_dir}/{year}-{month}.grib"
    get_grib_file(client, year, month, fname)
    data = parse_grib_file(fname)
    return collect_monthly_data(year, month, data)

def make_var(dates, data, tag, mode):
    t_data = np.empty((180, 360, len(dates)), dtype = float)
    if tag in ['asn', 'sd', 'tsn', 'fal', 'sde', 'snowc']:
        for (i, d) in enumerate(dates):
            t_data[:, :, i] = np.flip(data[d][tag][1:, :], axis = 0) #drop +90 degrees, make it go -90...89 and -180...179
    if tag in ['bsa', 'wsa', 'scf'] and mode == 1:
        for (i, d) in enumerate(dates):
            t_data[:, :, i] = np.flip(data[d][tag]["mode1"], axis = 0)
    if tag in ['bsa', 'wsa', 'scf'] and mode == 2:
        for (i, d) in enumerate(dates):
            t_data[:, :, i] = np.flip(data[d][tag]["mode2"][1:, :], axis = 0) #drop +90 degrees, make it go -90...89 and -180...179
    return t_data

def make_netcdf(mdata, edata, output_dir, output_name):
    lats = np.linspace(-90, 89, 180) #order in our output files, mark point values
    lons = np.linspace(-180, 179, 360) #order in our output files, mark point values
    dates_m = sorted(list(mdata.keys()))
    dates_e = sorted(list(edata.keys()))
    assert dates_m == dates_e
    #lat = 89.5...89.5, lon = -179.5...179.5 for MODIS mode1, mark the grid centers
    #lat = 90...-90, (need to drop +90), lon = -180...179 for MODIS mode2, aligns
    #lat[1] = 90...-90 (need to drop +90), lon = -180...179 (going to 359, then 0, then 179) for ERA5, aligns
    ds = xr.Dataset(
        data_vars={
            "snalb": (["lat", "lon", "date"], make_var(dates_m, edata, 'asn', 0)),
            "swe": (["lat", "lon", "date"], make_var(dates_m, edata, 'sd', 0)),
            "snd": (["lat", "lon", "date"], make_var(dates_m, edata, 'sde', 0)),
            "tsnow" : (["lat", "lon", "date"], make_var(dates_m, edata, 'tsn', 0)),
            "swa" : (["lat", "lon", "date"], make_var(dates_m, edata, 'fal', 0)),
            "snowc" : (["lat", "lon", "date"], make_var(dates_m, edata, 'snowc', 0)),
            "bsa_1" : (["lat", "lon", "date"], make_var(dates_m, mdata, 'bsa', 1)),
            "bsa_2" : (["lat", "lon", "date"], make_var(dates_m, mdata, 'bsa', 2)),
            "wsa_1" : (["lat", "lon", "date"], make_var(dates_m, mdata, 'wsa', 1)),
            "wsa_2" : (["lat", "lon", "date"], make_var(dates_m, mdata, 'wsa', 2)),
            "scf_1" : (["lat", "lon", "date"], make_var(dates_m, mdata, 'scf', 1)),
            "scf_2" : (["lat", "lon", "date"], make_var(dates_m, mdata, 'scf', 2)),
        },
        coords={
            "date": np.array(dates_m).astype('datetime64[ns]'),
            "lat": lats,
            "lon": lons,
        }
    )
    os.makedirs(output_dir, exist_ok=True)
    ds.to_netcdf(output_dir + f"/{output_name}.nc")

def gather_month_data(client, year, month, master_dict, era5dir, mod_dumpdir, nc_dir):
    output_tag = f"{year}" + "-" + f"{month:02d}"
    era5_data = get_era5land_data(client, year, month, era5dir)
    modis_data = get_modis_files(master_dict, year, month, mod_dumpdir)
    make_netcdf(modis_data, era5_data, nc_dir, output_tag)

def scrape_setup():
    client = Client()
    client.check_authentication()  # optional check
    earthaccess.login()
    all_mcd43c_files = earthaccess.search_data(
        short_name="MCD43C3",
        version="061",
        count=10000,
    ) #should be 9523 values
    master_date_dict = {
        get_MODIS_date_from_result(f) : f for f in all_mcd43c_files
    }
    return client, master_date_dict

def run_scrape():
    client, master_date_dict = scrape_setup()
    modis_dump_dir = "/home/acharbon/temp_modis"
    era5_dir = "/home/acharbon/sampo-data1/era5_data"
    nc_dir = "/home/acharbon/outputs/nc_files"
    missed_years = []
    progress = 0
    start_time = datetime.now()
    for year in range(2000, 2021):
        for month in range(1, 13):
            try:
                if date(year, month, 1) < date(2000, 3, 1):
                    continue
                if date(year, month, 1) > date(2020, 6, 1):
                    continue
                print(f"\n\n*** OBTAINING DATA FOR {year}-{month}: *****\n", flush = True)
                gather_month_data(client, year, month, master_date_dict, era5_dir, modis_dump_dir, nc_dir)
                progress += 1
                newtime = datetime.now()
                dt = newtime - start_time
                t_done = newtime + (dt*(244/progress) - dt)
                print(f"PROGRESS: {progress}/244")
                print(f"ESTIMATED COMPLETION: {t_done}")
            except:
                missed_years.append((year, month))
                print(f"****** ERROR WITH {year}-{month}*****", flush = True)
            else:
                print(f"{year}-{month} SUCCESSFULLY GATHERED.", flush = True)
    print("MISSED DATES: ")
    print(", ".join(missed_years))


#run_scrape()


# Suppress all the warnings
# get rid of leftover .idx files?

# 1 month of ERA5 land data -> ~6 min to grab it, 5-6 sec to process all of it
# 1 month of MODIS data -> ~30 sec to grab, 20-30 sec to process all of it

#earthaccess.login() #"persist = true" to setup once on cluster

#def get_hdf_file(file_data, write_dir):
#    earthaccess.download([file_data], local_path = write_dir)
#    return write_dir + "/" + os.listdir(write_dir)[0]

#def get_mod_file(filename, names, dump_dir):
#    #suppress warning and printing here? or do print something?
#    file = get_hdf_file(filename, dump_dir)
#    file_fields = parse_hdf_file(file, names)
#    return reduce_hdf_fields(file_fields)

#def agg_lat(field, axis):
#    #make 181 latitude points into 180:
#    arr_moved = np.moveaxis(field, axis, 0)
#    arr_temp = np.stack([arr_moved[:-1, ...], arr_moved[1:, ...]], axis = 0)
#    with np.errstate(invalid ='ignore'):
#        arr_downsampled = np.nanmean(arr_temp, axis = 0)
#    return np.moveaxis(arr_downsampled, 0, axis)

#def agg_lon(field, axis):
#    #lons run 0 - 359, while ours are 0.5, 1.5, ... 
#    field2 = np.roll(field, shift = -1, axis = axis)
#    arr_temp = np.stack([field, field2], axis = 0)
#    with np.errstate(invalid = 'ignore'):
#        arr_averaged = np.nanmean(arr_temp, axis = 0)
#    return arr_averaged

def plot_m(m):
    ny, nx = m.shape
    m_small = m.reshape(ny // 4, 4, nx // 4, 4)
    m_small = np.nanmean(m_small, axis = (1, 3))
    img = np.where(np.isnan(m_small), "  ", "**")
    for i in range(len(img)): #dim1
        for j in range(len(img[0])): #dim2
            print(img[i, j], end = "")
        print('\n')
    print("\n\n")

#request for the geopotential:
#month_request = {
#    "product_type": ["monthly_averaged_reanalysis"],
#    "variable": [
#        "surface_latent_heat_flux",
#        "surface_net_solar_radiation",
#        "surface_sensible_heat_flux",
#        "surface_solar_radiation_downwards"
#    ],
#    "year": [
#        "2000", "2001", "2002",
#        "2003", "2004", "2005",
#        "2006", "2007", "2008",
#        "2009", "2010", "2011",
#        "2012", "2013", "2014",
#        "2015", "2016", "2017",
#        "2018", "2019", "2020"
#    ],
#    "month": [
#        "01", "02", "03",
#        "04", "05", "06",
#        "07", "08", "09",
#        "10", "11", "12"
#    ],
#    "time": ["00:00"],
#    "grid": [1.0, 1.0],
#    "data_format": "netcdf",
#    "download_format": "unarchived"
#}

#request for geopotential (did not regrid):
#geo_request = {
#    "variable": ["geopotential"],
#    "grid": [1.0, 1.0],
#    "data_format": "grib",
#    "download_format": "unarchived"
#}

#request for glacier mask (did not regrid):
#geo_request = {
#    "variable": ["glacier_mask"],
#    "grid": [1.0, 1.0],
#    "data_format": "grib",
#    "download_format": "unarchived"
#}

def process_monthly_rad_data(path, output_dir):
    ds = xr.open_dataset(path)
    ds_data = {}
    for var in ['slhf', 'ssr', 'sshf', 'ssrd']:
        data = ds[var].values
        ddims = list(ds[var].dims)
        lon_idx = ddims.index('longitude')
        data_rolled = roll_lon(data, lon_idx)
        data_clima_aligned = np.flip(data_rolled[:, 1:, :], axis = 1) #flip latitude dimension (it's dim 1) and get rid of +90
        ds_data[var] = data_clima_aligned
    lats = np.linspace(-90, 89, 180) #order in our output files, mark point values
    lons = np.linspace(-180, 179, 360) #order in our output files, mark point values
    dates = ds["valid_time"].values
    new_ds = xr.Dataset(
        data_vars={
            "slhf": (["date", "lat", "lon"], ds_data['slhf']),
            "ssr": (["date", "lat", "lon"], ds_data['ssr']),
            "sshf": (["date", "lat", "lon",], ds_data['sshf']),
            "ssrd": (["date", "lat", "lon"], ds_data['ssrd']),
        },
        coords={
            "date": dates,
            "lat": lats,
            "lon": lons,
        }
    )
    os.makedirs(output_dir, exist_ok=True)
    new_ds.to_netcdf(output_dir + f"/era5_monthly_rad.nc")
    ds.close()

def process_month_folder(data_folder, output_folder):
    ds = xr.open_mfdataset(f"{data_folder}/*.nc")
    ds.to_netcdf(f"{output_folder}/all_snow_vars.nc")

def process_invariant(filename, npix, output_folder, filetype = "geopotential"):
    ds = xr.open_dataset(filename)
    lons = ds['longitude'].values
    invariant = (ds['z'].values / 9.80665) if filetype == "geopotential" else ds['glm'].values
    lon_idx = list(ds.dims).index('longitude')
    ny, nx = invariant.shape
    assert (ny-1) % npix == 0
    assert nx % npix == 0
    assert npix % 2 == 0
    nshift = npix // 2
    field_shift = np.roll(invariant, shift = len(lons)//2 + npix//2 - 1, axis = lon_idx) #positive shift bumps right-end to front, pushing 179.6-180 to the front to get averaged with -180 to -179.5
    field_big_inside = field_shift[nshift:-(nshift+1), :]
    field_big_top = field_shift[:nshift, :]
    field_big_bottom = field_shift[-(nshift+1):, :]
    field_inside_stack = field_big_inside.reshape(((ny-1) // npix) - 1, npix, nx // npix, npix)
    field_top_stack = field_big_top.reshape(1, nshift, nx // npix, npix)
    field_bottom_stack = field_big_bottom.reshape(1, nshift+1, nx // npix, npix)
    with np.errstate(invalid ='ignore'):
        field_small_inside = np.nanmean(field_inside_stack, axis = (1, 3)) #now -180, -179, ..., 179
        field_small_top = np.nanmean(field_top_stack, axis = (1, 3))
        field_small_bottom = np.nanmean(field_bottom_stack, axis = (1, 3))
    invar_avg = np.vstack([field_small_top, field_small_inside, field_small_bottom])
    invar_clima = np.flip(invar_avg[1:, :], axis = 0)
    lats = np.linspace(-90, 89, 180) #order in our output files, mark point values
    lons = np.linspace(-180, 179, 360) #order in our output files, mark point values
    varname = "elev" if filetype == "geopotential" else "glf"
    filename = "elevation" if filetype == "geopotential" else "glacialfrac"
    new_ds = xr.Dataset(
        data_vars={
            varname: (["lat", "lon"], invar_clima),
        },
        coords={
            "lat": lats,
            "lon": lons,
        }
    )
    os.makedirs(output_folder, exist_ok=True)
    new_ds.to_netcdf(output_folder + f"/{filename}.nc")
    ds.close()

def reshape_clm_field(var, ds, target_grid, regridder_cache, 
                      cons_field=[]):
    da = ds[var]
    # transpose
    other_dims = [d for d in da.dims if d not in ["lat", "lon"]]
    da = da.transpose(*other_dims, "lat", "lon")
    method = "conservative" if var in cons_field else "bilinear"
    key = (method, "lat", "lon")
    if key not in regridder_cache:
        regridder_cache[key] = xe.Regridder(
            da,
            target_grid,
            method,
        )
    regridder = regridder_cache[key]
    out = regridder(da)
    return out

def regrid_and_build_clm_outputs(file, year):
    ds = xr.open_dataset(file)
    target_grid = xr.Dataset({
        "lat": (["lat"], np.linspace(-90, 89, 180)),
        "lon": (["lon"], np.linspace(-180, 179, 360)),
    })
    regridder_cache = {}
    ds_out = xr.Dataset()
    for var in ds.data_vars:
        ds_out[var] = reshape_clm_field(
            var, ds, target_grid, regridder_cache
        )
    dates = np.array([date(year, i, 1) for i in range(1, 13)]).astype('datetime64[ns]')
    if "month" in ds_out.dims:
        ds_out = ds_out.rename({"month": "date"})
        ds_out = ds_out.assign_coords(date=dates)
    
    ds_out.to_netcdf(file.replace(".nc", "_clima.nc"))
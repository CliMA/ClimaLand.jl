module DataSourcesSMAP

import ClimaAnalysis
import NCDatasets
import Dates

export get_smap_obs_var_dict

"""
    get_smap_obs_var_dict() :: Dict{String, Function}

Return short_name => (start_date -> OutputVar) mapping for SMAP products.
"""
function get_smap_obs_var_dict()
    # EDIT these to your local files/varnames:
    smap_surface_path = "/data/SMAP/enhanced_L3/SMAP_L3_SM_P_E_2015_2024_daily.nc"
    sm_surface_var    = "sm_surface"    # e.g., SSMA or "soil_moisture"
    time_name, lat_name, lon_name = "time", "lat", "lon"

    surface_loader = function (_start_date::Dates.DateTime)
        # Load as ClimaAnalysis OutputVar
        # If you already have a helper like ClimaAnalysis.from_netcdf, use it.
        ds = NCDatasets.NCDataset(smap_surface_path)
        lon = vec(ds[lon_name][:]); lat = vec(ds[lat_name][:]); time = ds[time_name][:]
        data = ds[sm_surface_var][:]    # dims likely (time, lat, lon) or similar
        NCDatasets.close(ds)
        # Let ClimaAnalysis figure it out via a convenience helper:
        ov = ClimaAnalysis.wrap_netcdf(sm_surface_var; data=data, lon=lon, lat=lat, time=time)
        return ov
    end

    return Dict(
        "sm_surface" => surface_loader,
        # Add more entries if you want root-zone, etc.
        # "sm_root" => root_loader,
    )
end

end # module

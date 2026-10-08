# Write `ilamb_land_hist_site_rmse.csv`: the per-site FLUXNET2015 RMSE of the
# ILAMB land-hist models (https://www.ilamb.org/land-hist/), read from the
# per-model output files of the dashboard. `fluxnet_ilamb_rmse_plot.jl` averages
# these over the sites ClimaLand was run at. The CSV is checked in, so this only
# needs to be rerun if the dashboard changes.
#
# Usage:
#   julia --project=.buildkite experiments/integrated/generic_site/ilamb_land_hist_site_rmse.jl

import Downloads
import NCDatasets
import Statistics

const DASHBOARD = "https://www.ilamb.org/land-hist"
# ILAMB variable => dashboard page
const PAGES = (
    hfls = "HydrologyCycle/LatentHeat",
    hfss = "HydrologyCycle/SensibleHeat",
    rsus = "RadiationandEnergyCycle/SurfaceUpwardSWRadiation",
    rlus = "RadiationandEnergyCycle/SurfaceUpwardLWRadiation",
)
# The individual runs of the global table; the multi-model "Mean-*" rows are
# left out. land-hist has no JSBACH upward longwave.
const MODELS = [
    "$m-$f" for m in ("CLM", "ISBA-CTRIP", "JSBACH") for
    f in ("CRUJRA", "GSWP3", "Princeton")
]

function download_to(url, dir)
    path = joinpath(dir, basename(url))
    isfile(path) || Downloads.download(url, path)
    return path
end

# Site names are only in the benchmark data; the model files hold coordinates.
function benchmark_sites(var, dir)
    url = "https://www.ilamb.org/ILAMB-Data/DATA/$var/FLUXNET2015/$var.nc"
    path = download_to(url, mkpath(joinpath(dir, "DATA_$var")))
    return NCDatasets.NCDataset(path) do ds
        (string.(ds["site"][:]), Float64.(ds["lat"][:]), Float64.(ds["lon"][:]))
    end
end

function site_rmse(var, model, dir)
    url = "$DASHBOARD/$(PAGES[var])/FLUXNET2015/FLUXNET2015_$model.nc"
    path = try
        download_to(url, mkpath(joinpath(dir, String(var))))
    catch
        return nothing
    end
    return NCDatasets.NCDataset(path) do ds
        g = ds.group["MeanState"]
        # Sites without model output hold the default NetCDF fill value
        rmse =
            [v < 1e30 ? Float64(v) : missing for v in g["rmse_map_of_$var"][:]]
        # The dashboard's global value is the mean over sites, after a land
        # mask that drops one or two sites, so the two differ slightly
        @info "$var $model" mean_over_sites = Statistics.mean(skipmissing(rmse)) dashboard_global =
            g.group["scalars"]["RMSE global"][]
        (Float64.(g["lat"][:]), Float64.(g["lon"][:]), rmse)
    end
end

function main(outpath; dir = mktempdir())
    open(outpath, "w") do io
        println(io, join(["site", "variable", MODELS...], ","))
        for var in keys(PAGES)
            names, lats, lons = benchmark_sites(var, dir)
            values = fill("", length(names), length(MODELS))
            for (j, model) in enumerate(MODELS)
                r = site_rmse(var, model, dir)
                isnothing(r) && continue
                taken = falses(length(names))
                for (lat, lon, v) in zip(r...)
                    same_lon = .!taken .& (abs.(lons .- lon) .< 1e-3)
                    candidates =
                        findall(same_lon .& (abs.(lats .- lat) .< 1e-3))
                    # US-Myb: the benchmark latitude is 35.05 instead of 38.05
                    isempty(candidates) && (candidates = findall(same_lon))
                    isempty(candidates) &&
                        error("$var $model: no site at ($lat, $lon)")
                    # Co-located sites are in the same (alphabetical) order in
                    # both files
                    i = first(candidates)
                    taken[i] = true
                    ismissing(v) ||
                        (values[i, j] = string(round(v; sigdigits = 5)))
                end
            end
            for (i, name) in enumerate(names)
                all(isempty, values[i, :]) && continue
                println(io, join([name, String(var), values[i, :]...], ","))
            end
        end
    end
    return outpath
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(joinpath(@__DIR__, "ilamb_land_hist_site_rmse.csv"))
end

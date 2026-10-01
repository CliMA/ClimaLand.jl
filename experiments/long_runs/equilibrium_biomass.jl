# Equilibrium of the canopy carbon pools under the climate of a simulation, as initial
# conditions and compared with the XuSaatchi woody carbon.

import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import CairoMakie
import GeoMakie
import NCDatasets
import Statistics: mean
import ClimaLand
import ClimaLand.Canopy

const XUSAATCHI_PATH = joinpath(
    pkgdir(ClimaLand),
    "artifacts",
    "prognostic_carbon",
    "xusaatchi_woody_carbon_2deg.nc",
)

"""
    monthly(simdir, short_name; skip_months = 0)

The monthly means of the diagnostic `short_name` in the `ClimaAnalysis.SimDir`
`simdir` after its first `skip_months` months, as an array `(lon, lat, month)`.
"""
function monthly(simdir, short_name; skip_months = 0)
    var = get(simdir; short_name, reduction = "average", period = "1M")
    order = indexin(
        [
            ClimaAnalysis.Var.longitude_name(var),
            ClimaAnalysis.Var.latitude_name(var),
            ClimaAnalysis.Var.time_name(var),
        ],
        collect(keys(var.dims)),
    )
    return permutedims(var.data, order)[:, :, (skip_months + 1):end]
end

time_mean(x) = dropdims(mean(x; dims = 3); dims = 3)

"""
    equilibrium_pools(diagnostics_dir, parameters; skip_months = 0)

Equilibrium carbon pools with `parameters` (`Canopy.equilibrium_carbon_pools`) under
the climate of the simulation whose monthly diagnostics are in `diagnostics_dir`,
repeated, after its first `skip_months` months: the means of GPP, canopy leaf
respiration, air temperature, precipitation, C3 fraction, and the Q10 factor of the
monthly canopy temperature (air temperature where there is no canopy). Returns `lon`,
`lat` and a `NamedTuple` of `(lon, lat)` arrays: `C_leaf`, `C_stem`, `C_root`
(kg C m^-2), and the mean annual temperature `T_annual` (K) and precipitation
`P_annual` (m yr^-1).
"""
function equilibrium_pools(diagnostics_dir, parameters; skip_months = 0)
    simdir = ClimaAnalysis.SimDir(diagnostics_dir)
    kept(short_name) = monthly(simdir, short_name; skip_months)
    gpp = get(simdir; short_name = "gpp", reduction = "average", period = "1M")
    lon = ClimaAnalysis.longitudes(gpp)
    lat = ClimaAnalysis.latitudes(gpp)

    (; Q10, T_ref) = parameters
    T_air = kept("tair")
    T_canopy = kept("ct")
    T = @. ifelse(isnan(T_canopy), T_air, T_canopy)
    f_T = time_mean(@. Q10^((T - T_ref) / 10))
    T_annual = time_mean(T_air)
    # Precipitation is in kg m^-2 s^-1, negative downward
    P_annual = time_mean(kept("precip")) .* (-365 * 86400 / 1000)
    pools = Canopy.equilibrium_carbon_pools.(
        Ref(parameters),
        time_mean(kept("gpp")),
        time_mean(kept("crd")),
        f_T,
        T_annual,
        P_annual,
        time_mean(kept("fc3")),
    )
    return lon,
    lat,
    (;
        C_leaf = getproperty.(pools, :C_leaf),
        C_stem = getproperty.(pools, :C_stem),
        C_root = getproperty.(pools, :C_root),
        T_annual,
        P_annual,
    )
end

const IC_UNITS = Dict(
    :C_leaf => "kg m-2",
    :C_stem => "kg m-2",
    :C_root => "kg m-2",
    :T_annual => "K",
    :P_annual => "m yr-1",
    :LAI => "m2 m-2",
    :A0_daily => "mol m-2 day-1",
    :A0_annual => "mol m-2 yr-1",
    :precip_annual => "mol m-2 yr-1",
    :PET_annual => "mol m-2 yr-1",
    :VPDA0_annual => "Pa mol m-2 yr-1",
    :growing_days => "day",
    :A0c3_annual => "mol m-2 yr-1",
    :A0c4_annual => "mol m-2 yr-1",
    :GPPc3_annual => "mol m-2 yr-1",
)

"""
    write_initial_conditions(path, lon, lat, fields)

Writes the `(lon, lat)` arrays of the `NamedTuple` `fields`, named as the prognostic
variables they initialize, to a netCDF file at `path`.
"""
function write_initial_conditions(path, lon, lat, fields)
    NCDatasets.NCDataset(path, "c") do ds
        ds.attrib["title"] = "Initial conditions from a pre-industrial spin-up"
        NCDatasets.defVar(
            ds,
            "lon",
            lon,
            ("lon",);
            attrib = ["units" => "degrees_east"],
        )
        NCDatasets.defVar(
            ds,
            "lat",
            lat,
            ("lat",);
            attrib = ["units" => "degrees_north"],
        )
        for name in keys(fields)
            NCDatasets.defVar(
                ds,
                String(name),
                Float32.(fields[name]),
                ("lon", "lat");
                deflatelevel = 5,
                attrib = ["units" => IC_UNITS[name]],
            )
        end
    end
    return nothing
end

"""
    block_mean(lon, lat, x, block_lon, block_lat)

Mean of the finite values of `x[lon, lat]` in each block of the regular grid of block
centers `block_lon` and `block_lat`; NaN in blocks without any.
"""
function block_mean(lon, lat, x, block_lon, block_lat)
    Δlon = block_lon[2] - block_lon[1]
    Δlat = block_lat[2] - block_lat[1]
    total = zeros(length(block_lon), length(block_lat))
    count = zeros(Int, size(total))
    for (i, λ) in enumerate(lon), (j, φ) in enumerate(lat)
        isfinite(x[i, j]) || continue
        bi = mod(round(Int, (λ - block_lon[1]) / Δlon), length(block_lon)) + 1
        bj = clamp(
            round(Int, (φ - block_lat[1]) / Δlat) + 1,
            1,
            length(block_lat),
        )
        total[bi, bj] += x[i, j]
        count[bi, bj] += 1
    end
    return @. ifelse(count > 0, total / count, NaN)
end

"""
    lonlat_var(template, data, short_name, long_name, units)

An `OutputVar` with the `(lon, lat)` array `data` on the longitude-latitude grid of the
`OutputVar` `template`.
"""
function lonlat_var(template, data, short_name, long_name, units)
    order = indexin(
        collect(keys(template.dims)),
        [
            ClimaAnalysis.Var.longitude_name(template),
            ClimaAnalysis.Var.latitude_name(template),
        ],
    )
    return ClimaAnalysis.remake(
        template;
        data = permutedims(data, order),
        attributes = Dict(
            "short_name" => short_name,
            "long_name" => long_name,
            "units" => units,
        ),
    )
end

"""
    plot_against_benchmark(template, model, benchmark; short_name, titles,
        benchmark_name, units, colorrange, difference_range, path)

Maps the `(lon, lat)` arrays `model` and `benchmark`, titled `titles`, on the grid of
the `OutputVar` `template`, and their difference with its bias and RMSE over the
points both cover (weighted by area). Saves the figure at `path` and returns
`(; bias, rmse)`.
"""
function plot_against_benchmark(
    template,
    model,
    benchmark;
    short_name,
    titles,
    benchmark_name,
    units,
    colorrange,
    difference_range,
    path,
)
    difference = model .- benchmark
    both = isfinite.(difference)
    weights = [
        cosd(φ) for _ in ClimaAnalysis.longitudes(template),
        φ in ClimaAnalysis.latitudes(template)
    ][both]
    weighted_mean(x) = sum(weights .* x[both]) / sum(weights)
    bias = weighted_mean(difference)
    rmse = sqrt(weighted_mean(difference .^ 2))

    fig = CairoMakie.Figure(size = (1000, 1500))
    maps = Dict(:plot => Dict(:colormap => :viridis, :colorrange => colorrange))
    for (row, (data, long_name)) in enumerate(zip((model, benchmark), titles))
        viz.heatmap2D_on_globe!(
            fig,
            lonlat_var(template, data, short_name, long_name, units);
            p_loc = (row, 1),
            more_kwargs = maps,
        )
    end
    unit_label = isempty(units) ? "" : " ($units)"
    title = "Model - $benchmark_name$unit_label: bias $(round(bias; sigdigits = 2)), RMSE $(round(rmse; sigdigits = 2))"
    viz.heatmap2D_on_globe!(
        fig,
        lonlat_var(
            template,
            difference,
            short_name,
            "Model - $benchmark_name",
            units,
        );
        p_loc = (3, 1),
        more_kwargs = Dict(
            :plot => Dict(
                :colormap => CairoMakie.Reverse(:RdBu),
                :colorrange => difference_range,
            ),
            :axis => Dict(:title => title),
        ),
    )
    CairoMakie.save(path, fig)
    return (; bias, rmse)
end

"""
    plot_equilibrium_woody_carbon(lon, lat, C_stem; savedir, obs_path = XUSAATCHI_PATH)

Maps the equilibrium woody carbon `C_stem[lon, lat]` against the XuSaatchi woody carbon
in `obs_path`, in 2° blocks (`plot_against_benchmark`). Saves
`equilibrium_woody_carbon.png` in `savedir` and returns `(; bias, rmse)`.
"""
function plot_equilibrium_woody_carbon(
    lon,
    lat,
    C_stem;
    savedir,
    obs_path = XUSAATCHI_PATH,
)
    obs = ClimaAnalysis.OutputVar(obs_path, "woody_carbon")
    block_lon = ClimaAnalysis.longitudes(obs)
    block_lat = ClimaAnalysis.latitudes(obs)
    return plot_against_benchmark(
        obs,
        block_mean(lon, lat, C_stem, block_lon, block_lat),
        Float64.(coalesce.(obs.data, NaN));
        short_name = "cwood",
        titles = ("Equilibrium woody carbon", "XuSaatchi woody carbon"),
        benchmark_name = "XuSaatchi",
        units = "kg m^-2",
        colorrange = (0, 25),
        difference_range = (-15, 15),
        path = joinpath(savedir, "equilibrium_woody_carbon.png"),
    )
end

"""
    plot_spinup_state(template, fields; savedir)

Maps the `(lon, lat)` arrays of `fields`, a vector of
`(data, short_name, long_name, units, colorrange)`, on the grid of the `OutputVar`
`template`. Saves `spinup_state.png` in `savedir`.
"""
function plot_spinup_state(template, fields; savedir)
    fig = CairoMakie.Figure(size = (1000, 500 * length(fields)))
    for (row, (data, short_name, long_name, units, colorrange)) in
        enumerate(fields)
        viz.heatmap2D_on_globe!(
            fig,
            lonlat_var(template, data, short_name, long_name, units);
            p_loc = (row, 1),
            more_kwargs = Dict(
                :plot => Dict(:colormap => :viridis, :colorrange => colorrange),
            ),
        )
    end
    CairoMakie.save(joinpath(savedir, "spinup_state.png"), fig)
    return nothing
end

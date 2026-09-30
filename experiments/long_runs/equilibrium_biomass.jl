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
    plot_equilibrium_woody_carbon(lon, lat, C_stem; savedir, obs_path = XUSAATCHI_PATH)

Maps the equilibrium woody carbon `C_stem[lon, lat]`, the XuSaatchi woody carbon in
`obs_path`, and their difference, in 2° blocks, with the bias and RMSE over the blocks
both cover (weighted by area). Saves `equilibrium_woody_carbon.png` in `savedir` and
returns `(; bias, rmse)`.
"""
function plot_equilibrium_woody_carbon(
    lon,
    lat,
    C_stem;
    savedir,
    obs_path = XUSAATCHI_PATH,
)
    obs = ClimaAnalysis.OutputVar(obs_path, "woody_carbon")
    obs_data = Float64.(coalesce.(obs.data, NaN))
    block_lon = ClimaAnalysis.longitudes(obs)
    block_lat = ClimaAnalysis.latitudes(obs)
    model = block_mean(lon, lat, C_stem, block_lon, block_lat)
    difference = model .- obs_data
    both = isfinite.(difference)
    weights = [cosd(φ) for _ in block_lon, φ in block_lat][both]
    weighted_mean(x) = sum(weights .* x[both]) / sum(weights)
    bias = weighted_mean(difference)
    rmse = sqrt(weighted_mean(difference .^ 2))

    var(data, long_name) = ClimaAnalysis.remake(
        obs;
        data,
        attributes = Dict(
            "short_name" => "cwood",
            "long_name" => long_name,
            "units" => "kg m^-2",
        ),
    )
    maps = Dict(:plot => Dict(:colormap => :viridis, :colorrange => (0, 25)))
    fig = CairoMakie.Figure(size = (1000, 1500))
    viz.heatmap2D_on_globe!(
        fig,
        var(model, "Equilibrium woody carbon");
        p_loc = (1, 1),
        more_kwargs = maps,
    )
    viz.heatmap2D_on_globe!(
        fig,
        var(obs_data, "XuSaatchi woody carbon");
        p_loc = (2, 1),
        more_kwargs = maps,
    )
    title = "Model - XuSaatchi (kg m^-2): bias $(round(bias; digits = 2)), RMSE $(round(rmse; digits = 2))"
    viz.heatmap2D_on_globe!(
        fig,
        var(difference, "Model - XuSaatchi");
        p_loc = (3, 1),
        more_kwargs = Dict(
            :plot => Dict(
                :colormap => CairoMakie.Reverse(:RdBu),
                :colorrange => (-15, 15),
            ),
            :axis => Dict(:title => title),
        ),
    )
    CairoMakie.save(joinpath(savedir, "equilibrium_woody_carbon.png"), fig)
    return (; bias, rmse)
end

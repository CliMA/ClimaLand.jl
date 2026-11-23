"""
Depth bins (m) of the DIST panel of the water table leaderboard. The last one
holds every cell deeper than 10 m, including those without a water table in the
soil column.
"""
const _WATER_TABLE_BINS = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 10), (10, Inf)]

"""
    compute_water_table_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path;
        soil_depth,
        spin_up_months = 12,
        obs_path = ClimaLand.Artifacts.fan2017_water_table_depth_path(),
    )

Plot the simulated water table depth (`wtd`) against the estimate of Fan et al.
(2017), "Hydrologic regulation of plant rooting depth", PNAS 114 (40),
10572-10577, aggregated to the simulation grid: in each cell, the mean of the
depth at 30 arc-seconds capped at `soil_depth` (m), the depth of the simulated
soil column. The simulated water table cannot be deeper than the column, so that
mean is what the grid-mean water table of TOPMODEL stands for.

Columns are

- `SIM`, the simulated depth averaged over time,
- `ANN`, its bias against the estimate,
- `LAT`, zonal means of both,
- `DIST`, the land area by depth.

`obs_path` is a NetCDF file holding `wtd_mean_capped(lon, lat, cap)` in m, by
default the 1-degree aggregate of the estimate in the
`fan2017_water_table_depth_1degree` artifact, whose caps are 1, 2, ..., 100 m.
The leaderboard is skipped when the simulation has no `wtd`, or when
`soil_depth` is not one of the caps.
"""
function compute_water_table_leaderboard(
    leaderboard_base_path,
    diagnostics_folder_path;
    soil_depth,
    spin_up_months = 12,
    obs_path = ClimaLand.Artifacts.fan2017_water_table_depth_path(),
)
    sim_dir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    if !("wtd" in ClimaAnalysis.available_vars(sim_dir))
        @info "Skipping the water table leaderboard: no wtd in $diagnostics_folder_path"
        return nothing
    end
    sim_series = _load_partition_sim_var(sim_dir, "wtd", spin_up_months)
    year_range = _year_range(sim_series)
    sim = ClimaAnalysis.average_time(sim_series)
    obs_capped = ClimaAnalysis.OutputVar(obs_path, "wtd_mean_capped")
    caps = obs_capped.dims["cap"]
    cap = caps[argmin(abs.(caps .- soil_depth))]
    if !isapprox(cap, soil_depth; atol = 1e-3)
        @warn "Skipping the water table leaderboard: the observations are capped at $(first(caps)), ..., $(last(caps)) m, not at the soil depth, $soil_depth m"
        return nothing
    end
    obs = ClimaAnalysis.slice(obs_capped; cap)
    obs.attributes["short_name"] = "wtd"
    obs.attributes["long_name"] = "Water table depth (Fan et al. 2017)"
    obs = ClimaAnalysis.resampled_as(obs, sim; nan_threshold = 0.5)
    mask_fn = _intersect_mask(sim, obs)
    sim, obs = mask_fn(sim), mask_fn(obs)

    groups = ["SIM", "ANN", "LAT", "DIST"]
    fig = CairoMakie.Figure(; size = (600 * length(groups), 500))
    # The title and the column labels go in before any globe panel; see
    # `compute_partitioning_leaderboard`.
    titlelayout = CairoMakie.GridLayout(
        fig[-1, 1:length(groups)],
        halign = :center,
        tellwidth = false,
    )
    CairoMakie.Label(
        titlelayout[1, 1],
        "Water table depth\nexcluding spin up ($spin_up_months months)",
        fontsize = 40,
    )
    for (col_idx, group) in enumerate(groups)
        CairoMakie.Label(
            fig[0, col_idx],
            group,
            tellwidth = false,
            fontsize = 30,
        )
    end

    layout = fig[1, 1] = CairoMakie.GridLayout()
    ClimaAnalysis.Visualize.contour2D_on_globe!(
        layout,
        sim;
        more_kwargs = Dict(
            :plot => Dict(
                :levels => _water_table_levels(soil_depth),
                :extendhigh => :auto,
            ),
            :axis =>
                Dict(:title => "Water table depth (m), mean over $year_range"),
        ),
    )

    layout = fig[1, 2] = CairoMakie.GridLayout()
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        layout,
        sim,
        obs,
        cmap_extrema = (-10, 10),
    )

    lats, sim_zonal = _zonal_means(sim)
    _, obs_zonal = _zonal_means(obs)
    ax = CairoMakie.Axis(
        fig[1, 3],
        xlabel = "Water table depth (m)",
        ylabel = "Latitude (degrees)",
        yticks = -90:30:90,
        title = "Zonal mean",
    )
    CairoMakie.ylims!(ax, -90, 90)
    CairoMakie.lines!(
        ax,
        obs_zonal,
        lats;
        color = :black,
        linewidth = 4,
        label = "OBS",
    )
    CairoMakie.lines!(
        ax,
        sim_zonal,
        lats;
        color = :firebrick,
        linewidth = 4,
        label = "SIM",
    )
    CairoMakie.axislegend(ax, position = :rb, framevisible = false)

    sim_shares = _water_table_area_shares(sim)
    obs_shares = _water_table_area_shares(obs)
    nbins = length(_WATER_TABLE_BINS)
    ax = CairoMakie.Axis(
        fig[1, 4],
        xlabel = "Water table depth (m)",
        ylabel = "Land area (%)",
        xticks = (1:nbins, _water_table_bin_labels()),
        title = "Land area by depth",
    )
    CairoMakie.barplot!(
        ax,
        vcat(1:nbins, 1:nbins),
        vcat(obs_shares, sim_shares);
        dodge = vcat(fill(1, nbins), fill(2, nbins)),
        color = vcat(fill(:black, nbins), fill(:firebrick, nbins)),
    )
    CairoMakie.axislegend(
        ax,
        [
            CairoMakie.PolyElement(color = :black),
            CairoMakie.PolyElement(color = :firebrick),
        ],
        ["OBS", "SIM"],
        position = :lt,
        framevisible = false,
    )

    CairoMakie.save(
        joinpath(leaderboard_base_path, "WaterTable_leaderboard.png"),
        fig,
    )
    return nothing
end

"""
    _water_table_area_shares(var)

Return the percentage of the unmasked area of `var` in each depth bin of
`_WATER_TABLE_BINS`, weighting cells by the cosine of their latitude.
"""
function _water_table_area_shares(var)
    lats = ClimaAnalysis.latitudes(var)
    lat_dim = var.dim2index[ClimaAnalysis.latitude_name(var)]
    weights = reshape(
        cosd.(lats),
        ntuple(d -> d == lat_dim ? length(lats) : 1, ndims(var.data)),
    )
    weights = weights .* ones(size(var.data))
    valid = isfinite.(var.data)
    total = sum(weights[valid])
    return [
        100 * sum(weights[valid .& (lo .<= var.data .< hi)]) / total for
        (lo, hi) in _WATER_TABLE_BINS
    ]
end

_water_table_bin_labels() =
    [hi == Inf ? "> $(lo)" : "$(lo)–$(hi)" for (lo, hi) in _WATER_TABLE_BINS]

"""
    _water_table_levels(soil_depth)

Return the contour levels of the map of the water table depth, from 0 to
`soil_depth` with a round step, in at most 15 steps.
"""
function _water_table_levels(soil_depth)
    steps = (1, 2, 5, 10, 20, 50)
    step = steps[something(findfirst(s -> soil_depth / s <= 15, steps), 6)]
    return collect(0:step:soil_depth)
end

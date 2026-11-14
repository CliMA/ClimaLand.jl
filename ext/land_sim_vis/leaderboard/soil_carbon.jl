"""
    _column_weights(layer_edges, depth)

Thicknesses (m) of the parts of the layers bounded by `layer_edges` (m, from the
surface down) that lie above `depth` (m), from the top layer down.
"""
_column_weights(layer_edges, depth) = [
    max(min(layer_edges[k + 1], depth) - layer_edges[k], 0.0) for
    k in 1:(length(layer_edges) - 1)
]

"""
    _layer_edges(depths)

Edges (m, from the surface down) of the layers centred on `depths` (m, positive down),
each centre midway between the edges of its layer, as in the model's soil grid and in
SoilGrids.
"""
function _layer_edges(depths)
    edges = [0.0]
    for d in sort(depths)
        push!(edges, 2d - edges[end])
    end
    issorted(edges; lt = <=) ||
        error("Layer centres $depths are not midway between layer edges")
    return edges
end

"""
    _integrate_over_depth(var, depth)

The integral from the surface to `depth` (m) over the depth dimension of the
time-independent `var` (per m³), as a longitude-latitude `OutputVar` (per m²), over the
layers of `_layer_edges`.
"""
function _integrate_over_depth(var, depth)
    z_name = ClimaAnalysis.altitude_name(var)
    depths = abs.(ClimaAnalysis.altitudes(var))
    order = sortperm(depths)
    weights = zeros(length(depths))
    weights[order] .= _column_weights(_layer_edges(depths), depth)
    z_dim = findfirst(==(z_name), collect(keys(var.dims)))
    # Layers below `depth` are left out, not weighted by zero, as they may be NaN
    column = sum(
        selectdim(var.data, z_dim, k) .* weights[k] for
        k in eachindex(weights) if weights[k] > 0
    )
    surface = ClimaAnalysis.slice(
        var;
        Dict(Symbol(z_name) => first(ClimaAnalysis.altitudes(var)))...,
    )
    return ClimaAnalysis.remake(surface; data = column)
end

"""
    _land_total(var)

The land total of the longitude-latitude `var` (per m²) over its finite cells, in the
units of `var` times 1e12 m² (Pg C for kg C m⁻²).
"""
function _land_total(var)
    lons, lats = ClimaAnalysis.longitudes(var), ClimaAnalysis.latitudes(var)
    dlon, dlat =
        deg2rad(abs(lons[2] - lons[1])), deg2rad(abs(lats[2] - lats[1]))
    lat_dim =
        findfirst(==(ClimaAnalysis.latitude_name(var)), collect(keys(var.dims)))
    total = 0.0
    for I in CartesianIndices(var.data)
        x = var.data[I]
        isfinite(x) || continue
        total += x * 6.371e6^2 * cosd(lats[I[lat_dim]]) * dlon * dlat
    end
    return total * 1e-12
end

"""
    compute_soil_carbon_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path;
        spin_up_months = 12,
        max_ratio = 100.0,
    )

Compares the soil organic carbon a run would reach at steady state with SoilGrids, over
the top 1 m and 2 m, and saves `soil_carbon_steady_state.png` and `.csv`.

In each soil layer, the steady state of dSOC/dt = litter − Sm is estimated as
SOC × mean(litter)/mean(Sm), the means taken after the first `spin_up_months` months.
This is exact where Sm is proportional to SOC. Where mean(litter)/mean(Sm) exceeds
`max_ratio`, as in water-logged layers without O2, the layer has no steady state within
reach: it keeps its SOC, and the litter it receives is reported. The run needs the
monthly diagnostics `soc`, `scms` and `slit`, which a `LandModel` with carbon pools
writes by default; otherwise nothing is done.
"""
function compute_soil_carbon_leaderboard(
    leaderboard_base_path,
    diagnostics_folder_path;
    spin_up_months = 12,
    max_ratio = 100.0,
)
    sim_dir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    issubset(["soc", "scms", "slit"], ClimaAnalysis.available_vars(sim_dir)) ||
        return nothing
    @info "Steady-state soil carbon against SoilGrids"
    spinup_cutoff = spin_up_months * 30 * 86400.0
    function mean_after_spinup(short_name)
        var = get(sim_dir, short_name)
        ClimaAnalysis.times(var)[end] >= spinup_cutoff &&
            (var = ClimaAnalysis.window(var, "time", left = spinup_cutoff))
        return ClimaAnalysis.average_time(var)
    end
    obs = _preprocess_var(
        ClimaAnalysis.OutputVar(
            ClimaLand.Artifacts.soil_grids_ocd_artifact_path(),
            "ocd",
        ),
    )
    # SoilGrids covers fewer latitudes than the model
    obs_lats = extrema(ClimaAnalysis.latitudes(obs))
    soc, sm, litter = mean_after_spinup.(("soc", "scms", "slit"))
    covered = findall(
        lat -> obs_lats[1] <= lat <= obs_lats[2],
        ClimaAnalysis.latitudes(soc),
    )
    in_obs_lats(var) = ClimaAnalysis.window(
        var,
        "latitude";
        left = first(covered),
        right = last(covered),
        by = ClimaAnalysis.Index(),
    )
    soc, sm, litter = in_obs_lats.((soc, sm, litter))
    ratio = litter.data ./ sm.data
    no_steady_state = .!(ratio .<= max_ratio)
    steady = ifelse.(
        no_steady_state .| .!isfinite.(soc.data),
        soc.data,
        soc.data .* ratio,
    )
    soc_ss = ClimaAnalysis.remake(soc; data = steady)
    litter_trapped = ClimaAnalysis.remake(
        litter;
        data = ifelse.(
            no_steady_state .& isfinite.(litter.data),
            litter.data,
            0.0,
        ),
    )

    year = 365.25 * 86400
    rows =
        String["depth_m,soc_run_PgC,soc_steady_PgC,soilgrids_PgC,litter_PgC_yr,litter_no_steady_state_PgC_yr"]
    fig = CairoMakie.Figure(size = (2400, 900))
    for (row, depth) in enumerate((1.0, 2.0))
        run_col = _integrate_over_depth(soc, depth)
        sim_col = _integrate_over_depth(soc_ss, depth)
        # On the coordinates of `sim_col`: resampling keeps SoilGrids' Float32 ones
        obs_col = ClimaAnalysis.remake(
            sim_col;
            data = ClimaAnalysis.resampled_as(
                _integrate_over_depth(obs, depth),
                sim_col,
            ).data,
        )
        both = isfinite.(sim_col.data) .& isfinite.(obs_col.data)
        nan_out(v) = ClimaAnalysis.remake(v; data = ifelse.(both, v.data, NaN))
        run_col, sim_col, obs_col = nan_out.((run_col, sim_col, obs_col))
        litter_col = nan_out(_integrate_over_depth(litter, depth))
        trapped_col = nan_out(_integrate_over_depth(litter_trapped, depth))
        totals = _land_total.((run_col, sim_col, obs_col))
        litter_total = _land_total(litter_col) * year
        trapped_total = _land_total(trapped_col) * year
        push!(rows, join((depth, totals..., litter_total, trapped_total), ","))
        for (v, name) in ((sim_col, "SOC steady state"), (obs_col, "SoilGrids"))
            v.attributes["units"] = "kg C m^-2"
            v.attributes["short_name"] = name
            v.attributes["long_name"] = name
        end
        levels = collect(
            range(
                0,
                Statistics.quantile(filter(isfinite, vec(obs_col.data)), 0.98);
                length = 11,
            ),
        )
        for (col, v, title) in (
            (
                1,
                sim_col,
                Printf.@sprintf(
                    "Steady-state SOC, 0–%g m: %.0f Pg C (run mean %.0f)",
                    depth,
                    totals[2],
                    totals[1]
                )
            ),
            (
                2,
                obs_col,
                Printf.@sprintf(
                    "SoilGrids SOC, 0–%g m: %.0f Pg C",
                    depth,
                    totals[3]
                )
            ),
        )
            layout = fig[row, col] = CairoMakie.GridLayout()
            ClimaAnalysis.Visualize.contour2D_on_globe!(
                layout,
                v;
                more_kwargs = Dict(
                    :plot => Dict(:levels => levels, :extendhigh => :auto),
                    :axis => Dict(:title => title),
                ),
            )
        end
        layout = fig[row, 3] = CairoMakie.GridLayout()
        bias_extremum = Statistics.quantile(
            filter(isfinite, abs.(vec(sim_col.data .- obs_col.data))),
            0.95,
        )
        ClimaAnalysis.Visualize.plot_bias_on_globe!(
            layout,
            sim_col,
            obs_col;
            cmap_extrema = (-bias_extremum, bias_extremum),
        )
        lats, sim_zonal = _zonal_means(sim_col)
        _, obs_zonal = _zonal_means(obs_col)
        _, run_zonal = _zonal_means(run_col)
        ax = CairoMakie.Axis(
            fig[row, 4];
            xlabel = "SOC, 0–$(depth) m (kg C m^-2)",
            ylabel = "Latitude (degrees)",
            yticks = -90:30:90,
            title = Printf.@sprintf(
                "Zonal means; %.1f of %.1f Pg C/yr of litter\nin layers with no steady state",
                trapped_total,
                litter_total
            ),
        )
        CairoMakie.ylims!(ax, -90, 90)
        # A few cells far from steady state would flatten the other curves
        CairoMakie.xlims!(ax, 0, 3 * maximum(filter(isfinite, obs_zonal)))
        CairoMakie.lines!(
            ax,
            obs_zonal,
            lats;
            color = :black,
            linewidth = 3,
            label = "SoilGrids",
        )
        CairoMakie.lines!(
            ax,
            sim_zonal,
            lats;
            color = :firebrick,
            linewidth = 3,
            label = "steady state",
        )
        CairoMakie.lines!(
            ax,
            run_zonal,
            lats;
            color = :firebrick,
            linestyle = :dash,
            linewidth = 2,
            label = "run mean",
        )
        row == 1 && CairoMakie.axislegend(ax; position = :rb)
    end
    CairoMakie.Label(
        fig[0, 1:4],
        "Soil organic carbon at steady state (SOC × litter/Sm per layer, after $spin_up_months months) against SoilGrids",
        fontsize = 24,
        tellwidth = false,
    )
    CairoMakie.save(
        joinpath(leaderboard_base_path, "soil_carbon_steady_state.png"),
        fig,
    )
    write(
        joinpath(leaderboard_base_path, "soil_carbon_steady_state.csv"),
        join(rows, "\n") * "\n",
    )
    return nothing
end

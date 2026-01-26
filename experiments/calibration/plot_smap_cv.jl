"""
Plot SMAP Coefficient of Variation Results

This script reads the output from process_smap_cv.jl and creates
comprehensive visualization plots.

Usage:
    julia plot_smap_cv.jl /path/to/smap_cv_2015_2016.jld2
    
Or from Julia REPL:
    include("plot_smap_cv.jl")
    plot_smap_cv("/path/to/smap_cv_2015_2016.jld2")
"""

using JLD2
using Statistics
using Printf
import CairoMakie
import GeoMakie

"""
    plot_smap_cv(jld2_file::String; output_dir::String = dirname(jld2_file), downsample::Int = 2)

Create comprehensive visualization plots from SMAP CV analysis results.

# Arguments
- `jld2_file::String`: Path to the JLD2 file from process_smap_cv.jl
- `output_dir::String`: Directory to save plots (defaults to same directory as input file)
- `downsample::Int`: Downsampling factor for spatial plots (2 = every other pixel, 4 = every 4th pixel)

# Returns
- `fig`: CairoMakie Figure object with all plots
"""
function plot_smap_cv(jld2_file::String; output_dir::String = dirname(jld2_file), downsample::Int = 2)
    @info "Loading SMAP CV data from: $jld2_file"
    
    # Load data
    data = load(jld2_file)
    lons = data["lons"]
    lats = data["lats"]
    sm_cv_spatial = data["sm_cv_spatial"]
    sm_mean_spatial = data["sm_mean_spatial"]
    sm_std_spatial = data["sm_std_spatial"]
    summary_stats = data["summary_stats"]
    start_year = data["start_year"]
    stop_year = data["stop_year"]
    
    @info "Data loaded successfully"
    @info "Original grid size: $(size(sm_cv_spatial))"
    @info "Date range: $start_year - $stop_year"
    
    # Clean data: replace invalid values with NaN
    sm_cv_clean = copy(sm_cv_spatial)
    sm_mean_clean = copy(sm_mean_spatial)
    sm_std_clean = copy(sm_std_spatial)
    
    # Replace clearly invalid values
    sm_cv_clean[sm_cv_clean .< 0] .= NaN
    sm_cv_clean[sm_cv_clean .> 10] .= NaN
    sm_mean_clean[sm_mean_clean .< 0] .= NaN
    sm_mean_clean[sm_mean_clean .> 1] .= NaN
    
    # Clean coordinates - replace -9999 with NaN
    lons_clean = replace(lons, -9999.0f0 => NaN32)
    lats_clean = replace(lats, -9999.0f0 => NaN32)
    
    @info "Valid data points after cleaning: $(sum(.!isnan.(sm_cv_clean)))"
    @info "CV range: $(extrema(filter(!isnan, sm_cv_clean)))"
    
    # Downsample everything together
    if downsample > 1
        @info "Downsampling by factor of $downsample..."
        sm_cv_ds = sm_cv_clean[1:downsample:end, 1:downsample:end]
        sm_mean_ds = sm_mean_clean[1:downsample:end, 1:downsample:end]
        sm_std_ds = sm_std_clean[1:downsample:end, 1:downsample:end]
        lons_ds = lons_clean[1:downsample:end, 1:downsample:end]
        lats_ds = lats_clean[1:downsample:end, 1:downsample:end]
        
        @info "Downsampled grid: $(size(sm_cv_ds))"
        @info "Downsampled valid points: $(sum(.!isnan.(sm_cv_ds)))"
    else
        sm_cv_ds = sm_cv_clean
        sm_mean_ds = sm_mean_clean
        sm_std_ds = sm_std_clean
        lons_ds = lons_clean
        lats_ds = lats_clean
    end
    
    # Filter valid data for statistics
    valid_mask = .!isnan.(sm_cv_clean) .& .!isnan.(sm_mean_clean)
    valid_cv = sm_cv_clean[valid_mask]
    valid_mean = sm_mean_clean[valid_mask]
    valid_std = sm_std_clean[valid_mask]
    valid_lats = lats_clean[valid_mask]
    
    @info "Valid pixels: $(length(valid_cv)) / $(length(sm_cv_clean))"
    @info "Valid lats range: $(extrema(filter(!isnan, valid_lats)))"
    
    # Create output directory
    mkpath(output_dir)
    
    # Create multi-panel figure
    @info "Creating multi-panel figure..."
    fig = create_multipanel_figure(
        lons_ds, lats_ds,
        sm_cv_ds, sm_mean_ds, sm_std_ds,
        valid_cv, valid_mean, valid_std, valid_lats,
        summary_stats, start_year, stop_year
    )
    
    multipanel_file = joinpath(output_dir, "smap_cv_analysis_$(start_year)_$(stop_year).png")
    CairoMakie.save(multipanel_file, fig, px_per_unit = 1)
    @info "Saved multi-panel figure: $multipanel_file"
    
    # Create high-res map
    @info "Creating high-resolution global map..."
    fig_map = create_global_map(
        lons_clean, lats_clean, sm_cv_clean,
        valid_cv, start_year, stop_year
    )
    
    map_file = joinpath(output_dir, "smap_cv_map_$(start_year)_$(stop_year).png")
    CairoMakie.save(map_file, fig_map, px_per_unit = 1)
    @info "Saved high-res map: $map_file"
    
    # Create individual diagnostic plots
    @info "Creating diagnostic plots..."
    
    # Sample data for scatter plot if too many points
    max_scatter_points = 50000
    if length(valid_mean) > max_scatter_points
        @info "Sampling $(max_scatter_points) points for scatter plot (from $(length(valid_mean)))"
        sample_idx = rand(1:length(valid_mean), max_scatter_points)
        valid_mean_sample = valid_mean[sample_idx]
        valid_cv_sample = valid_cv[sample_idx]
    else
        valid_mean_sample = valid_mean
        valid_cv_sample = valid_cv
    end
    
    # 1. Histogram
    fig_hist = create_histogram(valid_cv, start_year, stop_year)
    hist_file = joinpath(output_dir, "smap_cv_histogram_$(start_year)_$(stop_year).png")
    CairoMakie.save(hist_file, fig_hist)
    @info "Saved histogram: $hist_file"
    
    # 2. Scatter plot (using sampled data)
    fig_scatter = create_scatter_plot(valid_mean_sample, valid_cv_sample, start_year, stop_year)
    scatter_file = joinpath(output_dir, "smap_cv_scatter_$(start_year)_$(stop_year).png")
    CairoMakie.save(scatter_file, fig_scatter)
    @info "Saved scatter plot: $scatter_file"
    
    # 3. Latitudinal profile
    fig_lat = create_latitudinal_profile(valid_lats, valid_cv, start_year, stop_year)
    lat_file = joinpath(output_dir, "smap_cv_latitude_$(start_year)_$(stop_year).png")
    CairoMakie.save(lat_file, fig_lat)
    @info "Saved latitudinal profile: $lat_file"
    
    @info "=" ^ 70
    @info "All plots created successfully!"
    @info "Output directory: $output_dir"
    @info "=" ^ 70
    
    return fig
end

"""
    create_multipanel_figure(...)

Create the main multi-panel summary figure.
"""
function create_multipanel_figure(
    lons, lats,
    sm_cv_spatial, sm_mean_spatial, sm_std_spatial,
    valid_cv, valid_mean, valid_std, valid_lats,
    summary_stats, start_year, stop_year
)
    fig = CairoMakie.Figure(size = (1200, 1000), fontsize = 12)
    
    # Panel 1: Match the standalone map settings exactly
    ax1 = GeoMakie.GeoAxis(
        fig[1, 1:2],
        dest = "+proj=robin",
        title = "Global Soil Moisture Coefficient of Variation ($start_year-$stop_year)",
    )
    
    hm1 = GeoMakie.surface!(
        ax1,
        lons, lats, sm_cv_spatial,
        colormap = :viridis,
        colorrange = (0, quantile(valid_cv, 0.95)),
    )
    
    # Use white coastlines like the standalone map (instead of black)
    GeoMakie.lines!(ax1, GeoMakie.coastlines(), color = :white, linewidth = 0.75)
    
    CairoMakie.Colorbar(
        fig[1, 3],
        label = "Coefficient of Variation (σ/μ)",
        height = CairoMakie.Relative(0.8),
    )
    
    # Panel 2: Histogram
    ax2 = CairoMakie.Axis(
        fig[2, 1],
        xlabel = "Coefficient of Variation",
        ylabel = "Frequency (log scale)",
        title = "Distribution of CV",
        yscale = log10,
    )
    
    CairoMakie.hist!(
        ax2,
        valid_cv,
        bins = 50,  # Reduced from 100
        color = (:steelblue, 0.7),
        strokewidth = 1,
        strokecolor = :black,
    )
    
    CairoMakie.vlines!(ax2, [median(valid_cv)], color = :red, linewidth = 2, linestyle = :dash, label = "Median")
    CairoMakie.vlines!(ax2, [mean(valid_cv)], color = :orange, linewidth = 2, linestyle = :dot, label = "Mean")
    CairoMakie.axislegend(ax2, position = :rt)
    
    # Panel 3: CV vs Mean SM scatter (use fewer bins)
    ax3 = CairoMakie.Axis(
        fig[2, 2],
        xlabel = "Mean Soil Moisture (m³/m³)",
        ylabel = "Coefficient of Variation",
        title = "CV vs Mean Soil Moisture",
    )
    
    CairoMakie.hexbin!(
        ax3,
        valid_mean,
        valid_cv,
        bins = 40,  # Reduced from 50
        colormap = :plasma,
        colorscale = log10,
    )
    
    # Panel 4: Latitudinal profile
    ax4 = CairoMakie.Axis(
        fig[3, 1:2],
        xlabel = "Latitude (°N)",
        ylabel = "Coefficient of Variation",
        title = "Latitudinal Profile",
    )
    
    # Filter out NaN latitudes first
    valid_lat_mask = .!isnan.(valid_lats)
    lats_clean = valid_lats[valid_lat_mask]
    cv_clean = valid_cv[valid_lat_mask]
    
    if length(lats_clean) > 0
        # Use actual latitude range
        lat_min = floor(minimum(lats_clean) / 5) * 5
        lat_max = ceil(maximum(lats_clean) / 5) * 5
        lat_bins = lat_min:5:lat_max
        lat_centers = (lat_bins[1:end-1] .+ lat_bins[2:end]) ./ 2
        cv_by_lat_mean = Float64[]
        cv_by_lat_std = Float64[]
        
        for i in 1:length(lat_bins)-1
            mask = (lats_clean .>= lat_bins[i]) .& (lats_clean .< lat_bins[i+1])
            if sum(mask) > 0
                push!(cv_by_lat_mean, mean(cv_clean[mask]))
                push!(cv_by_lat_std, std(cv_clean[mask]))
            else
                push!(cv_by_lat_mean, NaN)
                push!(cv_by_lat_std, NaN)
            end
        end
        
        valid_bins = .!isnan.(cv_by_lat_mean)
        
        if sum(valid_bins) > 0
            lat_centers_valid = lat_centers[valid_bins]
            cv_mean_valid = cv_by_lat_mean[valid_bins]
            cv_std_valid = cv_by_lat_std[valid_bins]
            
            CairoMakie.band!(
                ax4,
                lat_centers_valid,
                cv_mean_valid .- cv_std_valid,
                cv_mean_valid .+ cv_std_valid,
                color = (:steelblue, 0.3),
            )
            
            CairoMakie.lines!(
                ax4,
                lat_centers_valid,
                cv_mean_valid,
                color = :steelblue,
                linewidth = 3,
            )
        end
    end
    
    # Panel 5: Summary statistics
    ax5 = CairoMakie.Axis(fig[3, 3], aspect = 1)
    CairoMakie.hidedecorations!(ax5)
    CairoMakie.hidespines!(ax5)
    
    stats_text = """
    Summary Statistics
    ──────────────────
    Mean CV:      $(round(mean(valid_cv), digits=3))
    Median CV:    $(round(median(valid_cv), digits=3))
    Std Dev:      $(round(std(valid_cv), digits=3))
    Min CV:       $(round(minimum(valid_cv), digits=3))
    Max CV:       $(round(maximum(valid_cv), digits=3))
    
    10th %ile:    $(round(quantile(valid_cv, 0.10), digits=3))
    90th %ile:    $(round(quantile(valid_cv, 0.90), digits=3))
    
    Valid pixels: $(length(valid_cv))
    Coverage:     $(round(100*length(valid_cv)/length(sm_cv_spatial), digits=1))%
    """
    
    CairoMakie.text!(
        ax5,
        0.1, 0.9,
        text = stats_text,
        align = (:left, :top),
        fontsize = 11,
    )
    
    return fig
end

"""
    create_global_map(...)

Create high-resolution global map (using simple projection for speed).
"""
function create_global_map(lons, lats, sm_cv_spatial, valid_cv, start_year, stop_year)
    fig = CairoMakie.Figure(size = (1400, 800))
    
    ax = GeoMakie.GeoAxis(
        fig[1, 1],
        dest = "+proj=robin",
        title = "SMAP Soil Moisture Coefficient of Variation ($start_year-$stop_year)",
    )
    
    hm = GeoMakie.surface!(
        ax,
        lons, lats, sm_cv_spatial,
        colormap = :viridis,
        colorrange = (0, quantile(valid_cv, 0.95)),
        shading = false,
    )
    
    GeoMakie.lines!(ax, GeoMakie.coastlines(), color = :white, linewidth = 0.75)
    
    CairoMakie.Colorbar(
        fig[1, 2],
        hm,
        label = "Coefficient of Variation",
    )
    
    return fig
end

"""
    create_histogram(...)

Create CV distribution histogram.
"""
function create_histogram(valid_cv, start_year, stop_year)
    fig = CairoMakie.Figure(size = (800, 600))
    
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Coefficient of Variation",
        ylabel = "Frequency",
        title = "SMAP Soil Moisture CV Distribution ($start_year-$stop_year)",
    )
    
    CairoMakie.hist!(
        ax,
        valid_cv,
        bins = 100,
        color = (:steelblue, 0.7),
        strokewidth = 1,
        strokecolor = :black,
    )
    
    CairoMakie.vlines!(ax, [median(valid_cv)], color = :red, linewidth = 3, label = "Median")
    CairoMakie.vlines!(ax, [mean(valid_cv)], color = :orange, linewidth = 3, label = "Mean")
    
    CairoMakie.axislegend(ax, position = :rt)
    
    return fig
end

"""
    create_scatter_plot(...)

Create CV vs mean SM scatter plot.
"""
function create_scatter_plot(valid_mean, valid_cv, start_year, stop_year)
    fig = CairoMakie.Figure(size = (800, 600))
    
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Mean Soil Moisture (m³/m³)",
        ylabel = "Coefficient of Variation",
        title = "CV vs Mean Soil Moisture ($start_year-$stop_year)",
    )
    
    CairoMakie.hexbin!(
        ax,
        valid_mean,
        valid_cv,
        bins = 60,
        colormap = :plasma,
        colorscale = log10,
    )
    
    CairoMakie.Colorbar(fig[1, 2], label = "Count (log scale)")
    
    return fig
end

"""
    create_latitudinal_profile(...)

Create latitudinal profile plot.
"""
function create_latitudinal_profile(valid_lats, valid_cv, start_year, stop_year)
    fig = CairoMakie.Figure(size = (1000, 600))
    
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Latitude (°N)",
        ylabel = "Coefficient of Variation",
        title = "Latitudinal Profile of CV ($start_year-$stop_year)",
    )
    
    # Filter out NaN latitudes first
    valid_lat_mask = .!isnan.(valid_lats)
    lats_clean = valid_lats[valid_lat_mask]
    cv_clean = valid_cv[valid_lat_mask]
    
    if length(lats_clean) == 0
        @warn "No valid latitudes for latitudinal profile"
        return fig
    end
    
    # Use actual latitude range from cleaned data
    lat_min = floor(minimum(lats_clean) / 5) * 5
    lat_max = ceil(maximum(lats_clean) / 5) * 5
    lat_bins = lat_min:5:lat_max
    
    @info "Latitude bins: $lat_min to $lat_max ($(length(lat_bins)-1) bins)"
    
    lat_centers = (lat_bins[1:end-1] .+ lat_bins[2:end]) ./ 2
    cv_by_lat_mean = Float64[]
    cv_by_lat_std = Float64[]
    
    for i in 1:length(lat_bins)-1
        mask = (lats_clean .>= lat_bins[i]) .& (lats_clean .< lat_bins[i+1])
        n_in_bin = sum(mask)
        if n_in_bin > 0
            push!(cv_by_lat_mean, mean(cv_clean[mask]))
            push!(cv_by_lat_std, std(cv_clean[mask]))
        else
            push!(cv_by_lat_mean, NaN)
            push!(cv_by_lat_std, NaN)
        end
    end
    
    valid_bins = .!isnan.(cv_by_lat_mean)
    @info "Valid latitude bins: $(sum(valid_bins)) / $(length(cv_by_lat_mean))"
    
    if sum(valid_bins) == 0
        @warn "No valid latitude bins found!"
        return fig
    end
    
    lat_centers_valid = lat_centers[valid_bins]
    cv_mean_valid = cv_by_lat_mean[valid_bins]
    cv_std_valid = cv_by_lat_std[valid_bins]
    
    CairoMakie.band!(
        ax,
        lat_centers_valid,
        cv_mean_valid .- cv_std_valid,
        cv_mean_valid .+ cv_std_valid,
        color = (:steelblue, 0.3),
        label = "± 1 Std Dev",
    )
    
    CairoMakie.lines!(
        ax,
        lat_centers_valid,
        cv_mean_valid,
        color = :steelblue,
        linewidth = 3,
        label = "Mean CV",
    )
    
    CairoMakie.axislegend(ax, position = :rt)
    
    return fig
end

"""
    plot_mean_vs_cv_comparison(jld2_file::String; output_dir::String = dirname(jld2_file))

Create side-by-side comparison showing why CV is better than mean SM for heterogeneity.
"""
function plot_mean_vs_cv_comparison(jld2_file::String; output_dir::String = dirname(jld2_file))
    @info "Creating Mean vs CV comparison plots..."
    
    # Load data
    data = load(jld2_file)
    lons = replace(data["lons"], -9999.0f0 => NaN32)
    lats = replace(data["lats"], -9999.0f0 => NaN32)
    
    # Fix: Use array comprehension or direct assignment for conditional replacement
    sm_cv = copy(data["sm_cv_spatial"])
    sm_cv[(sm_cv .< 0) .| (sm_cv .> 10)] .= NaN
    
    sm_mean = copy(data["sm_mean_spatial"])
    sm_mean[(sm_mean .< 0) .| (sm_mean .> 1)] .= NaN
    
    sm_std = copy(data["sm_std_spatial"])
    sm_std[(sm_std .< 0) .| (sm_std .> 1)] .= NaN
    
    start_year = data["start_year"]
    stop_year = data["stop_year"]
    
    # Downsample for faster plotting
    ds = 2
    lons_ds = lons[1:ds:end, 1:ds:end]
    lats_ds = lats[1:ds:end, 1:ds:end]
    sm_cv_ds = sm_cv[1:ds:end, 1:ds:end]
    sm_mean_ds = sm_mean[1:ds:end, 1:ds:end]
    sm_std_ds = sm_std[1:ds:end, 1:ds:end]
    
    # Create figure with 3 panels
    fig = CairoMakie.Figure(size = (1600, 1400), fontsize = 14)
    
    # Panel 1: Mean Soil Moisture
    ax1 = GeoMakie.GeoAxis(
        fig[1, 1],
        dest = "+proj=robin",
        title = "A) Mean Soil Moisture ($start_year-$stop_year)",
    )
    
    hm1 = GeoMakie.surface!(
        ax1, lons_ds, lats_ds, sm_mean_ds,
        colormap = :Blues,
        colorrange = (0, 0.5),
    )
    GeoMakie.lines!(ax1, GeoMakie.coastlines(), color = :black, linewidth = 0.5)
    
    CairoMakie.Colorbar(fig[1, 2], hm1, label = "Mean SM (m³/m³)")
    
    # Panel 2: Standard Deviation
    ax2 = GeoMakie.GeoAxis(
        fig[2, 1],
        dest = "+proj=robin",
        title = "B) Standard Deviation of Soil Moisture",
    )
    
    hm2 = GeoMakie.surface!(
        ax2, lons_ds, lats_ds, sm_std_ds,
        colormap = :Oranges,
        colorrange = (0, 0.15),
    )
    GeoMakie.lines!(ax2, GeoMakie.coastlines(), color = :black, linewidth = 0.5)
    
    CairoMakie.Colorbar(fig[2, 2], hm2, label = "Std Dev SM (m³/m³)")
    
    # Panel 3: Coefficient of Variation
    ax3 = GeoMakie.GeoAxis(
        fig[3, 1],
        dest = "+proj=robin",
        title = "C) Coefficient of Variation (Normalized Variability)",
    )
    
    valid_cv = filter(!isnan, sm_cv_ds)
    hm3 = GeoMakie.surface!(
        ax3, lons_ds, lats_ds, sm_cv_ds,
        colormap = :viridis,
        colorrange = (0, quantile(valid_cv, 0.95)),
    )
    GeoMakie.lines!(ax3, GeoMakie.coastlines(), color = :white, linewidth = 0.75)
    
    CairoMakie.Colorbar(fig[3, 2], hm3, label = "CV (σ/μ)")
    
    # Add explanation text
    text_box = CairoMakie.Label(
        fig[4, 1:2],
        text = """
        Why CV is better than Mean or Std Dev for heterogeneity:
        
        • Panel A (Mean): Shows where soil is wet vs dry, but NOT how variable it is
        • Panel B (Std Dev): Shows absolute variability, but confounded by mean (wet regions have higher σ)
        • Panel C (CV): Shows RELATIVE variability, independent of mean moisture level
        
        Key insight: A desert with μ=0.05 and σ=0.03 has CV=0.6 (high variability)
                    A wetland with μ=0.50 and σ=0.03 has CV=0.06 (low variability)
        Both have same σ, but completely different temporal dynamics!
        """,
        fontsize = 12,
        halign = :left,
        valign = :top,
        tellwidth = false,
        tellheight = true,
    )
    
    output_file = joinpath(output_dir, "smap_mean_vs_cv_comparison_$(start_year)_$(stop_year).png")
    CairoMakie.save(output_file, fig, px_per_unit = 2)
    @info "Saved comparison plot: $output_file"
    
    return fig
end

"""
    plot_normalization_benefit(jld2_file::String; output_dir::String = dirname(jld2_file))

Show why normalization (CV) is necessary by comparing std dev vs CV.
"""
function plot_normalization_benefit(jld2_file::String; output_dir::String = dirname(jld2_file))
    @info "Creating normalization benefit plot..."
    
    # Load data
    data = load(jld2_file)
    sm_cv = copy(data["sm_cv_spatial"])
    sm_mean = copy(data["sm_mean_spatial"])
    sm_std = copy(data["sm_std_spatial"])
    start_year = data["start_year"]
    stop_year = data["stop_year"]
    
    # Clean data using boolean indexing
    valid_mask = .!isnan.(sm_cv) .& .!isnan.(sm_mean) .& .!isnan.(sm_std) .&
                 (sm_cv .> 0) .& (sm_cv .< 10) .&
                 (sm_mean .> 0) .& (sm_mean .< 1) .&
                 (sm_std .> 0) .& (sm_std .< 1)
    
    cv_valid = sm_cv[valid_mask]
    mean_valid = sm_mean[valid_mask]
    std_valid = sm_std[valid_mask]
    
    # Sample for plotting
    n_sample = min(100000, length(cv_valid))
    idx = rand(1:length(cv_valid), n_sample)
    
    # Create figure
    fig = CairoMakie.Figure(size = (1400, 600), fontsize = 14)
    
    # Panel A: Std Dev vs Mean (shows correlation)
    ax1 = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Mean Soil Moisture (m³/m³)",
        ylabel = "Standard Deviation (m³/m³)",
        title = "A) Problem: Std Dev correlates with Mean",
    )
    
    hb1 = CairoMakie.hexbin!(
        ax1,
        mean_valid[idx],
        std_valid[idx],
        bins = 50,
        colormap = :plasma,
        colorscale = log10,
    )
    
    # Add correlation coefficient
    r_std = cor(mean_valid, std_valid)
    CairoMakie.text!(
        ax1, 0.7, 0.95,
        text = @sprintf("Correlation: r = %.3f", r_std),
        align = (:left, :top),
        space = :relative,
        fontsize = 14,
        color = :red,
    )
    
    # Panel B: CV vs Mean (shows independence)
    ax2 = CairoMakie.Axis(
        fig[1, 2],
        xlabel = "Mean Soil Moisture (m³/m³)",
        ylabel = "Coefficient of Variation (σ/μ)",
        title = "B) Solution: CV is independent of Mean",
    )
    
    hb2 = CairoMakie.hexbin!(
        ax2,
        mean_valid[idx],
        cv_valid[idx],
        bins = 50,
        colormap = :plasma,
        colorscale = log10,
    )
    
    # Add correlation coefficient
    r_cv = cor(mean_valid, cv_valid)
    CairoMakie.text!(
        ax2, 0.7, 0.95,
        text = @sprintf("Correlation: r = %.3f", r_cv),
        align = (:left, :top),
        space = :relative,
        fontsize = 14,
        color = :green,
    )
    
    # Fixed: Remove colorscale from Colorbar, just use the hexbin's colormap
    CairoMakie.Colorbar(
        fig[1, 3], 
        hb1,  # Reference the hexbin plot
        label = "Count (log scale)",
    )
    
    # Add explanation
    CairoMakie.Label(
        fig[2, 1:2],
        text = """
        Key Finding: Std Dev has correlation r=$(round(r_std, digits=3)) with Mean, but CV has r=$(round(r_cv, digits=3))
        
        This means:
        • Left: Wetter regions automatically have higher std dev (artifact of measurement scale)
        • Right: CV removes this artifact - variability is now comparable across all moisture levels
        • CV allows fair comparison: Is Sahara Desert (dry, variable) more heterogeneous than Amazon (wet, stable)?
        """,
        fontsize = 12,
        halign = :left,
        tellheight = true,
    )
    
    output_file = joinpath(output_dir, "smap_normalization_benefit_$(start_year)_$(stop_year).png")
    CairoMakie.save(output_file, fig, px_per_unit = 2)
    @info "Saved normalization benefit plot: $output_file"
    
    return fig
end

"""
    plot_regional_examples(jld2_file::String; output_dir::String = dirname(jld2_file))

Show example regions where same std dev means different things.
"""
function plot_regional_examples(jld2_file::String; output_dir::String = dirname(jld2_file))
    @info "Creating regional example plots..."
    
    # Load data
    data = load(jld2_file)
    lons = data["lons"]
    lats = data["lats"]
    sm_cv = copy(data["sm_cv_spatial"])
    sm_mean = copy(data["sm_mean_spatial"])
    sm_std = copy(data["sm_std_spatial"])
    start_year = data["start_year"]
    stop_year = data["stop_year"]
    
    # Clean data
    sm_cv[(sm_cv .< 0) .| (sm_cv .> 10)] .= NaN
    sm_mean[(sm_mean .< 0) .| (sm_mean .> 1)] .= NaN
    sm_std[(sm_std .< 0) .| (sm_std .> 1)] .= NaN
    
    # Define regions of interest (approximate lat/lon boxes)
    regions = [
        ("Sahara Desert", -10.0, 30.0, 15.0, 35.0),
        ("Amazon Rainforest", -75.0, -50.0, -10.0, 5.0),
        ("US Great Plains", -105.0, -95.0, 35.0, 45.0),
        ("Southeast Asia", 95.0, 110.0, 10.0, 25.0),
    ]
    
    fig = CairoMakie.Figure(size = (1400, 1000), fontsize = 12)
    
    stats_text = String[]
    
    for (name, lon_min, lon_max, lat_min, lat_max) in regions
        # Create mask for region - handle 2D coordinate arrays and fill values
        mask = (lons .>= lon_min) .& (lons .<= lon_max) .& 
               (lats .>= lat_min) .& (lats .<= lat_max) .&
               (lons .!= -9999.0f0) .& (lats .!= -9999.0f0) .&  # Exclude fill values
               .!isnan.(sm_cv) .& .!isnan.(sm_mean) .& .!isnan.(sm_std)
        
        n_pixels = sum(mask)
        @info "Region '$name': found $n_pixels pixels"
        
        if n_pixels < 100
            @warn "Skipping '$name' - insufficient data ($n_pixels pixels)"
            continue
        end
        
        region_cv = sm_cv[mask]
        region_mean = sm_mean[mask]
        region_std = sm_std[mask]
        
        # Calculate statistics
        μ = mean(region_mean)
        σ = mean(region_std)
        cv = mean(region_cv)
        
        push!(stats_text, 
            @sprintf("%-20s: μ=%.3f, σ=%.3f, CV=%.3f (n=%d)", 
                     name, μ, σ, cv, n_pixels))
    end
    
    if isempty(stats_text)
        @warn "No valid regions found!"
        push!(stats_text, "No valid regions found - check coordinate ranges")
    end
    
    # Create summary table
    ax = CairoMakie.Axis(fig[1, 1], aspect = DataAspect())
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)
    CairoMakie.xlims!(ax, 0, 1)
    CairoMakie.ylims!(ax, 0, 1)
    
    summary = """
    Regional Comparison: Why CV matters
    ═══════════════════════════════════════════════════════════
    
    $(join(stats_text, "\n"))
    
    
    Interpretation:
    ──────────────
    • If two regions have similar σ but different μ, their CV will differ dramatically
    • High CV = high temporal variability relative to baseline moisture level
    • Low CV = stable moisture conditions (wetlands, rainforests)
    • CV allows fair comparison across different climate zones
    
    Example:
    • Desert: μ=0.05, σ=0.03 → CV=0.6 (highly variable, 60% fluctuation)
    • Wetland: μ=0.50, σ=0.03 → CV=0.06 (stable, only 6% fluctuation)
    • Both have SAME σ=0.03, but VERY different temporal dynamics!
    """
    
    CairoMakie.text!(
        ax, 0.05, 0.95,
        text = summary,
        align = (:left, :top),
        fontsize = 13,
        font = :regular,
    )
    
    output_file = joinpath(output_dir, "smap_regional_examples_$(start_year)_$(stop_year).png")
    CairoMakie.save(output_file, fig)
    @info "Saved regional examples: $output_file"
    
    return fig
end

# Command line execution
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia plot_smap_cv.jl <path_to_jld2_file> [downsample_factor]")
        println("Example: julia plot_smap_cv.jl /path/to/smap_cv_2015_2016.jld2 2")
        println("         downsample_factor: 1=full res (slow), 2=half res (default), 4=quarter res (fast)")
        exit(1)
    end
    
    jld2_file = ARGS[1]
    
    if !isfile(jld2_file)
        error("File not found: $jld2_file")
    end
    
    # Optional: downsample factor
    downsample = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
    
    # Optional: specify output directory
    output_dir = length(ARGS) >= 3 ? ARGS[3] : dirname(jld2_file)
    
    # Create plots
    #plot_smap_cv(jld2_file; output_dir = output_dir, downsample = downsample)
    
    # Create additional comparison plots
    plot_mean_vs_cv_comparison(jld2_file; output_dir = output_dir)
    plot_normalization_benefit(jld2_file; output_dir = output_dir)
    plot_regional_examples(jld2_file; output_dir = output_dir)
end
"""
Process SMAP soil moisture data to calculate coefficient of variation (CV)
and generate summary statistics and visualizations.

This script:
1. Reads SMAP L3_SM_P_E data using smap_utils.jl
2. Calculates temporal CV for each pixel using calculate_drydown.jl
3. Generates summary statistics (mean, median, percentiles)
4. Creates maps of average CV across the globe
5. Exports results to NetCDF and JLD2 formats
"""

using Dates
using Statistics
using JLD2
using Printf
import NCDatasets
import CairoMakie
import GeoMakie

# Load local utilities
include("smap_utils.jl")
include("calculate_drydown.jl")

"""
    process_smap_cv(start_year::Int, stop_year::Int; 
                    data_dir::String, 
                    output_dir::String)

Main processing function for SMAP CV calculation.

# Arguments
- `start_year::Int`: First year to process (SMAP starts April 2015)
- `stop_year::Int`: Last year to process
- `data_dir::String`: Directory containing SMAP HDF5 files
- `output_dir::String`: Directory for output files

# Returns
- `Dict`: Summary statistics including global mean, median, spatial patterns
"""
function process_smap_cv(
    start_year::Int,
    stop_year::Int;
    data_dir::String = get(ENV, "SMAP_DATA_PATH", "/resnick/scratch/egreich/SMAP/SMAP_L3_SM_P_E"),
    output_dir::String = joinpath(@__DIR__, "output"),
)
    @info "=" ^ 70
    @info "SMAP Soil Moisture Coefficient of Variation Analysis"
    @info "=" ^ 70
    @info "Processing years: $start_year - $stop_year"
    @info "Data directory: $data_dir"
    @info "Output directory: $output_dir"
    
    # Create output directory if needed
    mkpath(output_dir)
    
    # Set date range
    # SMAP mission started April 31, 2015
    start_date = DateTime(max(start_year, 2015), start_year == 2015 ? 4 : 1, 1)
    stop_date = DateTime(stop_year, 12, 31)
    
    @info "Date range: $(Dates.format(start_date, "yyyy-mm-dd")) to $(Dates.format(stop_date, "yyyy-mm-dd"))"
    
    # Step 1: Find SMAP files
    @info "\n" * "=" ^ 70
    @info "Step 1: Finding SMAP files..."
    @info "=" ^ 70
    
    smap_files = find_smap_files(start_date, stop_date; data_dir=data_dir)  # Remove smap_utils. prefix
    
    if isempty(smap_files)
        error("No SMAP files found in $data_dir for the specified date range")
    end
    
    @info "Found $(length(smap_files)) SMAP files"
    @info "First file: $(basename(smap_files[1]))"
    @info "Last file: $(basename(smap_files[end]))"
    
    # Step 2: Calculate CV using calculate_drydown.jl
    @info "\n" * "=" ^ 70
    @info "Step 2: Calculating coefficient of variation..."
    @info "=" ^ 70
    
    lons, lats, sm_std_spatial, sm_cv_spatial, sm_mean_spatial = 
        calculate_sm_variability_from_smap(smap_files, start_date, stop_date)
    
    # Step 3: Generate summary statistics
    @info "\n" * "=" ^ 70
    @info "Step 3: Computing summary statistics..."
    @info "=" ^ 70
    
    summary_stats = compute_cv_summary_stats(
        lons, lats, 
        sm_cv_spatial, 
        sm_std_spatial, 
        sm_mean_spatial
    )
    
    # Print summary
    print_summary_stats(summary_stats)
    
    # Step 4: Save results
    @info "\n" * "=" ^ 70
    @info "Step 4: Saving results..."
    @info "=" ^ 70
    
    # Save to JLD2 (for Julia)
    jld2_file = joinpath(
        output_dir, 
        "smap_cv_$(start_year)_$(stop_year).jld2"
    )
    
    jldsave(
        jld2_file;
        lons,
        lats,
        sm_cv_spatial,
        sm_std_spatial,
        sm_mean_spatial,
        summary_stats,
        start_year,
        stop_year,
        start_date,
        stop_date,
        n_files = length(smap_files),
    )
    
    @info "Saved JLD2 file: $jld2_file"
    
    # Save to NetCDF (for broader compatibility)
    nc_file = joinpath(
        output_dir,
        "smap_cv_$(start_year)_$(stop_year).nc"
    )
    
    save_to_netcdf(
        nc_file,
        lons, lats,
        sm_cv_spatial,
        sm_std_spatial,
        sm_mean_spatial,
        start_date, stop_date
    )
    
    @info "Saved NetCDF file: $nc_file"
    
    # Save summary statistics to text file
    txt_file = joinpath(
        output_dir,
        "smap_cv_summary_$(start_year)_$(stop_year).txt"
    )
    
    save_summary_to_txt(txt_file, summary_stats, start_year, stop_year)
    @info "Saved summary text file: $txt_file"
    
    # Step 5: Create visualization plots
    @info "\n" * "=" ^ 70
    @info "Step 5: Creating visualization plots..."
    @info "=" ^ 70
    
    try
        create_cv_plots(
            lons, lats,
            sm_cv_spatial,
            sm_mean_spatial,
            output_dir,
            start_year, stop_year
        )
    catch e
        @warn "Failed to create plots (continuing without visualization)" exception=(e, catch_backtrace())
    end
    
    @info "\n" * "=" ^ 70
    @info "Processing complete!"
    @info "=" ^ 70
    
    return summary_stats
end

"""
    compute_cv_summary_stats(lons, lats, sm_cv_spatial, sm_std_spatial, sm_mean_spatial)

Compute comprehensive summary statistics for the CV field.
"""
function compute_cv_summary_stats(lons, lats, sm_cv_spatial, sm_std_spatial, sm_mean_spatial)
    # Filter valid (non-NaN) values
    valid_cv = filter(!isnan, sm_cv_spatial)
    valid_std = filter(!isnan, sm_std_spatial)
    valid_mean = filter(!isnan, sm_mean_spatial)
    
    stats = Dict{String, Any}()
    
    # Global statistics for CV
    stats["cv_global_mean"] = mean(valid_cv)
    stats["cv_global_median"] = median(valid_cv)
    stats["cv_global_std"] = std(valid_cv)
    stats["cv_global_min"] = minimum(valid_cv)
    stats["cv_global_max"] = maximum(valid_cv)
    stats["cv_p10"] = quantile(valid_cv, 0.10)
    stats["cv_p25"] = quantile(valid_cv, 0.25)
    stats["cv_p75"] = quantile(valid_cv, 0.75)
    stats["cv_p90"] = quantile(valid_cv, 0.90)
    
    # Global statistics for mean SM
    stats["sm_mean_global_mean"] = mean(valid_mean)
    stats["sm_mean_global_median"] = median(valid_mean)
    stats["sm_mean_global_std"] = std(valid_mean)
    
    # Global statistics for std SM
    stats["sm_std_global_mean"] = mean(valid_std)
    stats["sm_std_global_median"] = median(valid_std)
    
    # Data availability
    total_pixels = length(sm_cv_spatial)
    valid_pixels = length(valid_cv)
    stats["n_total_pixels"] = total_pixels
    stats["n_valid_pixels"] = valid_pixels
    stats["data_coverage_pct"] = 100.0 * valid_pixels / total_pixels
    
    # Spatial heterogeneity metrics
    stats["spatial_cv_of_cv"] = std(valid_cv) / mean(valid_cv)
    stats["spatial_range_cv"] = maximum(valid_cv) - minimum(valid_cv)
    
    return stats
end

"""
    print_summary_stats(stats::Dict)

Print summary statistics in a formatted table.
"""
function print_summary_stats(stats::Dict)
    println("\n" * "=" ^ 70)
    println("SUMMARY STATISTICS - Coefficient of Variation (CV)")
    println("=" ^ 70)
    
    @printf("Global Mean CV:            %.4f\n", stats["cv_global_mean"])
    @printf("Global Median CV:          %.4f\n", stats["cv_global_median"])
    @printf("Global Std Dev CV:         %.4f\n", stats["cv_global_std"])
    @printf("Global Min CV:             %.4f\n", stats["cv_global_min"])
    @printf("Global Max CV:             %.4f\n", stats["cv_global_max"])
    println("-" ^ 70)
    @printf("10th Percentile:           %.4f\n", stats["cv_p10"])
    @printf("25th Percentile:           %.4f\n", stats["cv_p25"])
    @printf("75th Percentile:           %.4f\n", stats["cv_p75"])
    @printf("90th Percentile:           %.4f\n", stats["cv_p90"])
    
    println("\n" * "=" ^ 70)
    println("MEAN SOIL MOISTURE STATISTICS")
    println("=" ^ 70)
    @printf("Global Mean SM:            %.4f m³/m³\n", stats["sm_mean_global_mean"])
    @printf("Global Median SM:          %.4f m³/m³\n", stats["sm_mean_global_median"])
    @printf("Global Std Dev SM:         %.4f m³/m³\n", stats["sm_mean_global_std"])
    
    println("\n" * "=" ^ 70)
    println("DATA COVERAGE")
    println("=" ^ 70)
    @printf("Total Pixels:              %d\n", stats["n_total_pixels"])
    @printf("Valid Pixels:              %d\n", stats["n_valid_pixels"])
    @printf("Coverage:                  %.2f%%\n", stats["data_coverage_pct"])
    
    println("\n" * "=" ^ 70)
    println("SPATIAL VARIABILITY")
    println("=" ^ 70)
    @printf("Spatial CV of CV:          %.4f\n", stats["spatial_cv_of_cv"])
    @printf("Spatial Range (max-min):   %.4f\n", stats["spatial_range_cv"])
    
    println("=" ^ 70 * "\n")
end

"""
    save_to_netcdf(filename, lons, lats, sm_cv, sm_std, sm_mean, start_date, stop_date)

Save CV and related fields to NetCDF format for use in other tools.
"""
function save_to_netcdf(
    filename::String,
    lons::Array,
    lats::Array,
    sm_cv::Array,
    sm_std::Array,
    sm_mean::Array,
    start_date::DateTime,
    stop_date::DateTime,
)
    NCDatasets.Dataset(filename, "c") do ds
        # Define dimensions
        nlat, nlon = size(sm_cv)
        NCDatasets.defDim(ds, "lat", nlat)
        NCDatasets.defDim(ds, "lon", nlon)
        
        # Define variables
        nc_lat = NCDatasets.defVar(ds, "latitude", Float32, ("lat", "lon"))
        nc_lon = NCDatasets.defVar(ds, "longitude", Float32, ("lat", "lon"))
        nc_cv = NCDatasets.defVar(ds, "sm_cv", Float32, ("lat", "lon"))
        nc_std = NCDatasets.defVar(ds, "sm_std", Float32, ("lat", "lon"))
        nc_mean = NCDatasets.defVar(ds, "sm_mean", Float32, ("lat", "lon"))
        
        # Add attributes
        nc_lat.attrib["units"] = "degrees_north"
        nc_lat.attrib["long_name"] = "Latitude"
        
        nc_lon.attrib["units"] = "degrees_east"
        nc_lon.attrib["long_name"] = "Longitude"
        
        nc_cv.attrib["units"] = "1"
        nc_cv.attrib["long_name"] = "Coefficient of Variation of Soil Moisture"
        nc_cv.attrib["description"] = "Temporal CV = std(SM) / mean(SM) for each pixel"
        
        nc_std.attrib["units"] = "m3/m3"
        nc_std.attrib["long_name"] = "Standard Deviation of Soil Moisture"
        
        nc_mean.attrib["units"] = "m3/m3"
        nc_mean.attrib["long_name"] = "Mean Soil Moisture"
        
        # Global attributes
        ds.attrib["title"] = "SMAP Soil Moisture Coefficient of Variation"
        ds.attrib["source"] = "SMAP L3_SM_P_E Enhanced v6"
        ds.attrib["start_date"] = Dates.format(start_date, "yyyy-mm-dd")
        ds.attrib["stop_date"] = Dates.format(stop_date, "yyyy-mm-dd")
        ds.attrib["created"] = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        ds.attrib["grid"] = "EASE-Grid 2.0, 9 km resolution"
        
        # Write data
        nc_lat[:, :] = lats
        nc_lon[:, :] = lons
        nc_cv[:, :] = sm_cv
        nc_std[:, :] = sm_std
        nc_mean[:, :] = sm_mean
    end
end

"""
    save_summary_to_txt(filename, stats, start_year, stop_year)

Save summary statistics to a human-readable text file.
"""
function save_summary_to_txt(filename::String, stats::Dict, start_year::Int, stop_year::Int)
    open(filename, "w") do io
        println(io, "=" ^ 70)
        println(io, "SMAP Soil Moisture Coefficient of Variation Analysis")
        println(io, "Processing Period: $start_year - $stop_year")
        println(io, "Generated: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "=" ^ 70)
        println(io)
        
        println(io, "COEFFICIENT OF VARIATION (CV) STATISTICS")
        println(io, "-" ^ 70)
        @printf(io, "Global Mean:               %.6f\n", stats["cv_global_mean"])
        @printf(io, "Global Median:             %.6f\n", stats["cv_global_median"])
        @printf(io, "Global Std Dev:            %.6f\n", stats["cv_global_std"])
        @printf(io, "Global Minimum:            %.6f\n", stats["cv_global_min"])
        @printf(io, "Global Maximum:            %.6f\n", stats["cv_global_max"])
        println(io)
        
        println(io, "PERCENTILES")
        println(io, "-" ^ 70)
        @printf(io, "10th Percentile:           %.6f\n", stats["cv_p10"])
        @printf(io, "25th Percentile:           %.6f\n", stats["cv_p25"])
        @printf(io, "75th Percentile:           %.6f\n", stats["cv_p75"])
        @printf(io, "90th Percentile:           %.6f\n", stats["cv_p90"])
        println(io)
        
        println(io, "MEAN SOIL MOISTURE STATISTICS")
        println(io, "-" ^ 70)
        @printf(io, "Global Mean:               %.6f m³/m³\n", stats["sm_mean_global_mean"])
        @printf(io, "Global Median:             %.6f m³/m³\n", stats["sm_mean_global_median"])
        @printf(io, "Global Std Dev:            %.6f m³/m³\n", stats["sm_mean_global_std"])
        println(io)
        
        println(io, "DATA COVERAGE")
        println(io, "-" ^ 70)
        @printf(io, "Total Pixels:              %d\n", stats["n_total_pixels"])
        @printf(io, "Valid Pixels:              %d\n", stats["n_valid_pixels"])
        @printf(io, "Coverage Percentage:       %.2f%%\n", stats["data_coverage_pct"])
        println(io)
        
        println(io, "SPATIAL VARIABILITY METRICS")
        println(io, "-" ^ 70)
        @printf(io, "Spatial CV of CV:          %.6f\n", stats["spatial_cv_of_cv"])
        @printf(io, "Spatial Range:             %.6f\n", stats["spatial_range_cv"])
        println(io)
        println(io, "=" ^ 70)
    end
end

"""
    create_cv_plots(lons, lats, sm_cv_spatial, sm_mean_spatial, output_dir, start_year, stop_year)

Create visualization plots for soil moisture CV analysis.

Generates:
1. Global map of CV
2. Histogram of CV distribution
3. CV vs mean soil moisture scatter plot
4. Latitudinal profile of CV
"""
function create_cv_plots(
    lons::Array,
    lats::Array,
    sm_cv_spatial::Array,
    sm_mean_spatial::Array,
    output_dir::String,
    start_year::Int,
    stop_year::Int,
)
    @info "Creating visualization plots..."
    
    # Filter valid data
    valid_mask = .!isnan.(sm_cv_spatial) .& .!isnan.(sm_mean_spatial)
    valid_cv = sm_cv_spatial[valid_mask]
    valid_mean = sm_mean_spatial[valid_mask]
    valid_lats = lats[valid_mask]
    
    # Set up figure with multiple panels
    fig = CairoMakie.Figure(size = (1400, 1200), fontsize = 14)
    
    # Panel 1: Global map of CV
    ax1 = GeoMakie.GeoAxis(
        fig[1, 1:2],
        dest = "+proj=robin",  # Robinson projection for global view
        title = "Global Soil Moisture Coefficient of Variation ($start_year-$stop_year)",
        xlabel = "Longitude",
        ylabel = "Latitude",
    )
    
    # Create heatmap
    # Note: GeoMakie works best with 1D lon/lat vectors, so we need to handle EASE-Grid carefully
    # For simplicity, we'll use a regular heatmap with coastlines overlay
    hm1 = CairoMakie.heatmap!(
        ax1,
        lons[:, 1],  # Assuming lons are repeated across rows
        lats[1, :],  # Assuming lats are repeated across columns
        sm_cv_spatial',
        colormap = :viridis,
        colorrange = (0, quantile(valid_cv, 0.95)),  # Cap at 95th percentile for better visibility
        nan_color = :lightgray,
    )
    
    # Add coastlines
    GeoMakie.lines!(ax1, GeoMakie.coastlines(), color = :black, linewidth = 0.5)
    
    # Add colorbar
    CairoMakie.Colorbar(
        fig[1, 3],
        hm1,
        label = "Coefficient of Variation (σ/μ)",
        height = CairoMakie.Relative(0.8),
    )
    
    # Panel 2: Histogram of CV distribution
    ax2 = CairoMakie.Axis(
        fig[2, 1],
        xlabel = "Coefficient of Variation",
        ylabel = "Frequency (log scale)",
        title = "Distribution of CV across all pixels",
        yscale = log10,
    )
    
    CairoMakie.hist!(
        ax2,
        valid_cv,
        bins = 100,
        color = (:steelblue, 0.7),
        strokewidth = 1,
        strokecolor = :black,
    )
    
    # Add vertical lines for percentiles
    CairoMakie.vlines!(ax2, [median(valid_cv)], color = :red, linewidth = 2, linestyle = :dash, label = "Median")
    CairoMakie.vlines!(ax2, [mean(valid_cv)], color = :orange, linewidth = 2, linestyle = :dot, label = "Mean")
    CairoMakie.axislegend(ax2, position = :rt)
    
    # Panel 3: CV vs Mean SM scatter plot
    ax3 = CairoMakie.Axis(
        fig[2, 2],
        xlabel = "Mean Soil Moisture (m³/m³)",
        ylabel = "Coefficient of Variation",
        title = "CV vs Mean Soil Moisture",
    )
    
    # Density scatter (hexbin) for better visualization with many points
    CairoMakie.hexbin!(
        ax3,
        valid_mean,
        valid_cv,
        bins = 50,
        colormap = :plasma,
        colorscale = log10,
    )
    
    # Panel 4: Latitudinal profile
    ax4 = CairoMakie.Axis(
        fig[3, 1:2],
        xlabel = "Latitude (°N)",
        ylabel = "Coefficient of Variation",
        title = "Latitudinal Profile of CV (mean ± std)",
    )
    
    # Bin by latitude (5-degree bins)
    lat_bins = -90:5:90
    lat_centers = (lat_bins[1:end-1] .+ lat_bins[2:end]) ./ 2
    cv_by_lat_mean = Float64[]
    cv_by_lat_std = Float64[]
    
    for i in 1:length(lat_bins)-1
        mask = (valid_lats .>= lat_bins[i]) .& (valid_lats .< lat_bins[i+1])
        if sum(mask) > 0
            push!(cv_by_lat_mean, mean(valid_cv[mask]))
            push!(cv_by_lat_std, std(valid_cv[mask]))
        else
            push!(cv_by_lat_mean, NaN)
            push!(cv_by_lat_std, NaN)
        end
    end
    
    # Remove NaN bins
    valid_bins = .!isnan.(cv_by_lat_mean)
    lat_centers_valid = lat_centers[valid_bins]
    cv_mean_valid = cv_by_lat_mean[valid_bins]
    cv_std_valid = cv_by_lat_std[valid_bins]
    
    # Plot mean with error band
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
        label = "Mean CV ± Std Dev",
    )
    
    CairoMakie.axislegend(ax4, position = :rt)
    
    # Panel 5: Summary statistics box
    ax5 = CairoMakie.Axis(
        fig[3, 3],
        aspect = 1,
    )
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
        fontsize = 12,
        font = :regular,  # Changed from :mono to :regular
    )
    
    # Save figure
    plot_file = joinpath(output_dir, "smap_cv_analysis_$(start_year)_$(stop_year).png")
    CairoMakie.save(plot_file, fig, px_per_unit = 2)
    @info "Saved plot: $plot_file"
    
    # Also save individual high-res map
    fig_map = CairoMakie.Figure(size = (1600, 900))
    ax_map = GeoMakie.GeoAxis(
        fig_map[1, 1],
        dest = "+proj=robin",
        title = "SMAP Soil Moisture Coefficient of Variation ($start_year-$stop_year)",
    )
    
    hm_map = CairoMakie.heatmap!(
        ax_map,
        lons[:, 1],
        lats[1, :],
        sm_cv_spatial',
        colormap = :viridis,
        colorrange = (0, quantile(valid_cv, 0.95)),
        nan_color = :lightgray,
    )
    
    GeoMakie.lines!(ax_map, GeoMakie.coastlines(), color = :white, linewidth = 0.75)
    
    CairoMakie.Colorbar(
        fig_map[1, 2],
        hm_map,
        label = "Coefficient of Variation",
    )
    
    map_file = joinpath(output_dir, "smap_cv_map_$(start_year)_$(stop_year).png")
    CairoMakie.save(map_file, fig_map, px_per_unit = 2)
    @info "Saved high-res map: $map_file"
    
    return fig
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line arguments
    if length(ARGS) >= 2
        start_year = parse(Int, ARGS[1])
        stop_year = parse(Int, ARGS[2])
    else
        # Default: process 2015-2023 (full SMAP record as of 2024)
        start_year = 2015
        stop_year = 2023
    end
    
    # Optional: override data/output directories
    data_dir = length(ARGS) >= 3 ? ARGS[3] : get(ENV, "SMAP_DATA_PATH", "/resnick/scratch/egreich/SMAP/SMAP_L3_SM_P_E")
    output_dir = length(ARGS) >= 4 ? ARGS[4] : joinpath(@__DIR__, "output")
    
    # Run processing
    summary_stats = process_smap_cv(
        start_year, 
        stop_year;
        data_dir = data_dir,
        output_dir = output_dir
    )
end
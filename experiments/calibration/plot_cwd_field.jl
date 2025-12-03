using JLD2
using CairoMakie

"""
Plot the CWD field from the generated JLD2 file.
"""
function plot_cwd_field(cwd_file::String; output_file::String = "cwd_field.png")
    @info "Loading CWD field from" cwd_file
    
    # Load the data
    data = JLD2.load(cwd_file)
    lons = data["lons"]
    lats = data["lats"]
    CWD_spatial = data["CWD_spatial"]
    start_year = data["start_year"]
    stop_year = data["stop_year"]
    
    @info "Creating plot..."
    
    # Create figure
    fig = Figure(size = (1400, 700))
    ax = Axis(fig[1, 1],
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        title = "Climatic Water Deficit (CWD) $(start_year)-$(stop_year)",
        aspect = DataAspect()
    )
    
    # Reverse latitude for correct orientation (North at top)
    lats_reversed = reverse(lats)
    CWD_reversed = reverse(CWD_spatial, dims=2)
    
    # Plot the CWD field
    hm = heatmap!(ax, lons, lats_reversed, CWD_reversed,
        colormap = :turbo,  # or try :viridis, :plasma, :inferno
        colorrange = (0, maximum(CWD_spatial))
    )
    
    # Add coastlines/grid
    xlims!(ax, 0, 360)
    ylims!(ax, -90, 90)
    
    # Add colorbar
    Colorbar(fig[1, 2], hm, label = "CWD (mm/year)")
    
    # Save figure
    save(output_file, fig)
    @info "Saved plot to" output_file
    
    # Display some statistics
    @info "CWD Statistics:" minimum=minimum(CWD_spatial) maximum=maximum(CWD_spatial) mean=sum(CWD_spatial)/length(CWD_spatial)
    
    return fig
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    cwd_file = joinpath(@__DIR__, "cwd_era5_2015_2023.jld2")
    if !isfile(cwd_file)
        @error "CWD file not found" cwd_file
        exit(1)
    end
    plot_cwd_field(cwd_file)
end
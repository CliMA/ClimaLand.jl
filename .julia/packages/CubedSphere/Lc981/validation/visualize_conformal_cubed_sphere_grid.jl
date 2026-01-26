using Rotations
using CubedSphere
using CubedSphere: UniformSpacing, GeometricSpacing, ExponentialSpacing
using CubedSphere.SphericalGeometry
using DelimitedFiles
using CairoMakie

function write_output_to_file_1D(output_directory, x, y, file_name)
    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path)
        mkdir(path)
    end
    cd(path)

    file_name *= ".curve"
    outputfile = open(file_name, "w")

    write(outputfile, "#phi\n")
    for i in eachindex(x)
        write(outputfile, string(x[i], " ", y[i], "\n"))
        # The above line is equivalent to println(outputfile, "$(x[i]) $(y[i])")
    end

    close(outputfile)
    cd(cwd)
end

function read_output_from_file_1D(output_directory, file_name)
    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path)
        mkdir(path)
    end
    cd(path)

    data = []
    count = 1
    open(file_name, "r") do infile
        for line in eachline(infile)
            if count != 1
                push!(data, line)
            end
            count += 1
        end
    end
    data = readdlm(IOBuffer(join(data, "\n")))

    N = size(data, 1)
    x = zeros(N)
    y = zeros(N)
    for i in 1:N
        x[i] = data[i,1]
        y[i] = data[i,2]
    end

    cd(cwd)

    return (x, y)
end

output_directory = "visualize_conformal_cubed_sphere_grid"
path = joinpath(pwd(), output_directory)
isdir(path) || mkdir(path)

function create_single_line_or_scatter_plot(size, plot_type, x, y, axis_kwargs, title, plot_kwargs, file_name;
                                            figure_padding = (0, 5, 0, 0), # (left, right, bottom, top)
                                            specify_x_limits = false, x_limits = [0, 0], specify_y_limits = false,
                                            y_limits = [0, 0], tight_x_axis = false, tight_y_axis = false,
                                            format = ".png", output_directory = "visualize_conformal_cubed_sphere_grid")
    fig = Figure(size = size, figure_padding = figure_padding)
    ax = Axis(fig[1,1]; axis_kwargs...)

    if plot_type == "line_plot"
        lines!(ax, x, y, linewidth = plot_kwargs.linewidth, color = plot_kwargs.linecolor)
    elseif plot_type == "scatter_plot"
        scatter!(ax, x, y, marker = plot_kwargs.marker, markersize = plot_kwargs.markersize,
                 color = plot_kwargs.linecolor)
    elseif plot_type == "scatter_line_plot"
        scatterlines!(ax, x, y, linewidth = plot_kwargs.linewidth, marker = plot_kwargs.marker,
                      markersize = plot_kwargs.markersize, color = plot_kwargs.linecolor)
    end
    ax.title = title

    if specify_x_limits
        xlims!(ax, x_limits...)
    elseif tight_x_axis
        xlims!(ax, extrema(x)...)
    end

    if specify_y_limits
        ylims!(ax, y_limits...)
    elseif tight_y_axis
        ylims!(ax, extrema(y)...)
    end

    save(output_directory * "/" * file_name * format, fig)
end

function single_line_or_scatter_plot_example()
    x = range(0, 2π, length = 100)
    y = sin.(x)

    size = (750, 750)

    axis_kwargs = (xlabel = "x", ylabel = "sin(x)", xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5,
                   yticklabelsize = 17.5, xlabelpadding = 10, ylabelpadding = 10, aspect = 1, titlesize = 27.5,
                   titlegap = 15, titlefont = :bold)
    title = "sin(x) vs x"
    plot_kwargs = (linewidth = 2, linecolor = :black, marker = :rect, markersize = 10)

    plot_types = ["line_plot", "scatter_plot", "scatter_line_plot"]
    file_names = ["LinePlotExample", "ScatterPlotExample", "ScatterLinePlotExample"]

    for i in 1:3
        plot_type = plot_types[i]
        file_name = file_names[i]
        create_single_line_or_scatter_plot(size, plot_type, x, y, axis_kwargs, title, plot_kwargs, file_name;
                                           format = ".pdf")
    end
end

single_line_or_scatter_plot_example()

function create_multiple_line_or_scatter_plots(size, plot_type, x, ys, axis_kwargs, title, plot_kwargs, file_name;
                                               figure_padding = (0, 5, 0, 0), # (left, right, bottom, top)
                                               specify_x_limits = false, x_limits = [0, 0], specify_y_limits = false,
                                               y_limits = [0, 0], tight_x_axis = false, tight_y_axis = false,
                                               halign = :center, valign = :center, format = ".png",
                                               output_directory = "visualize_conformal_cubed_sphere_grid")
    fig = Figure(size = size, figure_padding = figure_padding)
    ax = Axis(fig[1,1]; axis_kwargs...)

    for i in axes(ys, 1)
        if plot_type == "line_plot"
            lines!(ax, x, ys[i, :], linewidth = plot_kwargs.linewidth, color = plot_kwargs.linecolors[i],
                   label = plot_kwargs.labels[i])
        elseif plot_type == "scatter_plot"
            scatter!(ax, x, ys[i, :], marker = plot_kwargs.markers[i], markersize = plot_kwargs.markersize,
                     color = plot_kwargs.linecolors[i], label = plot_kwargs.labels[i])
        elseif plot_type == "scatter_line_plot"
            scatterlines!(ax, x, ys[i, :], linewidth = plot_kwargs.linewidth, marker = plot_kwargs.markers[i],
                          markersize = plot_kwargs.markersize, color = plot_kwargs.linecolors[i],
                          label = plot_kwargs.labels[i])
        end
    end

    axislegend(ax; position = (halign, valign))

    ax.title = title

    if specify_x_limits
        xlims!(ax, x_limits...)
    elseif tight_x_axis
        xlims!(ax, extrema(x)...)
    end

    if specify_y_limits
        ylims!(ax, y_limits...)
    elseif tight_y_axis
        ylims!(ax, extrema(y)...)
    end

    save(output_directory * "/" * file_name * format, fig)
end

function multiple_line_or_scatter_plots_example()
    x = range(0, 2π, length = 100)
    y1 = sin.(x)
    y2 = cos.(x)
    y = zeros(2, length(x))
    y[1, :] = y1
    y[2, :] = y2

    size = (750, 750)

    axis_kwargs = (xlabel = "x", ylabel = "sin(x)", xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5,
                   yticklabelsize = 17.5, xlabelpadding = 10, ylabelpadding = 10, aspect = 1, titlesize = 27.5,
                   titlegap = 15, titlefont = :bold)
    title = "sin(x) vs x"
    plot_kwargs = (linewidth = 2, linecolors = [:red, :blue], markers = [:rect, :rect], markersize = 10,
                   labels = ["y1", "y2"])

    plot_types = ["line_plot", "scatter_plot", "scatter_line_plot"]
    file_names = ["LinePlotExample_2", "ScatterPlotExample_2", "ScatterLinePlotExample_2"]

    for i in 1:3
        plot_type = plot_types[i]
        file_name = file_names[i]
        create_multiple_line_or_scatter_plots(size, plot_type, x, y, axis_kwargs, title, plot_kwargs, file_name;
                                              format = ".pdf", halign = :left, valign = :bottom)
    end
end

multiple_line_or_scatter_plots_example()

function visualize_conformal_cubed_sphere_panel_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, color;
                                                   figure_padding = (0, 5, 0, 0)) # (left, right, bottom, top)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(size = (750, 750), figure_padding = figure_padding)

    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Conformal Cubed Sphere Panel", axis_kwargs_2D...)
    hide_decorations && hidedecorations!(ax2D)
    wireframe!(ax2D, X, Y, Z, color = color)

    return fig
end

function visualize_conformal_cubed_sphere_panel_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, color;
                                                   figure_padding = (25, 0, 35, 50)) # (left, right, bottom, top)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(size = (750, 750), figure_padding = figure_padding)

    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)),
                 title = "Conformal Cubed Sphere Panel", axis_kwargs_3D...)
    hide_decorations && hidedecorations!(ax3D)
    wireframe!(ax3D, X, Y, Z, color = color)

    return fig
end

function visualize_conformal_cubed_sphere_panel_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, color)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(size = (1500, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Conformal Cubed Sphere Panel", axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)),
                 title = "Conformal Cubed Sphere Panel", axis_kwargs_3D...)

    for ax in [ax2D, ax3D]
        hide_decorations && hidedecorations!(ax)
        wireframe!(ax, X, Y, Z, color = color)
    end

    colgap!(fig.layout, 60)

    return fig
end

function visualize_conformal_cubed_sphere_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, colors;
                                             figure_padding = (25, 100, 0, 0), # (left, right, bottom, top)
                                             title = "Conformal Cubed Sphere: 2D Projection")
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(size = (700, 750), figure_padding = figure_padding)

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title, axis_kwargs_2D...)
    hide_decorations && hidedecorations!(ax2D)
    wireframe!(ax2D, X, Y, Z, color = colors[1])

    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

    for (i, R) in enumerate(rotations)

        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)

        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end

        wireframe!(ax2D, X′, Y′, Z′, color = colors[i + 1])

    end

    return fig
end

axis_kwargs_2D = (xlabel = "x", ylabel = "y", xlabelsize = 22.5, ylabelsize = 22.5, xticksize = 8, yticksize = 8,
                  xticklabelsize = 17.5, yticklabelsize = 17.5, xlabelpadding = 10, ylabelpadding = 10,
                  xticklabelpad = 10, yticklabelpad = 10, aspect = 1, titlesize = 27.5, titlegap = 15,
                  titlefont = :bold)

axis_kwargs_3D = (xlabel = "x", ylabel = "y", zlabel = "z", xlabelsize = 22.5, ylabelsize = 22.5, zlabelsize = 22.5,
                  xticksize = 8, yticksize = 8, zticksize = 8, xticklabelsize = 17.5, yticklabelsize = 17.5,
                  zticklabelsize = 17.5, xticklabelpad = 10, yticklabelpad = 10, zticklabelpad = 10, titlesize = 27.5,
                  titlegap = 15, titlefont = :bold)

function visualize_conformal_cubed_sphere_2D(X, Y, Z, file_name;
                                             figure_padding = (25, 100, 0, 0), # (left, right, bottom, top)
                                             title = "Conformal Cubed Sphere: 2D Projection",
                                             output_directory = "visualize_conformal_cubed_sphere_grid")
    hide_decorations = false

    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]

    fig = Figure(size = (725, 750), figure_padding = figure_padding)

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title, axis_kwargs_2D...)
    hide_decorations && hidedecorations!(ax2D)
    wireframe!(ax2D, X, Y, Z, color = colors[1])

    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

    for (i, R) in enumerate(rotations)
        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)

        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end

        wireframe!(ax2D, X′, Y′, Z′, color = colors[i + 1])
    end

    save(output_directory * "/" * file_name, fig)
end

function visualize_conformal_cubed_sphere_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, colors, alphas;
                                             figure_padding = (25, 0, 0, 50), # (left, right, bottom, top)
                                             title = "Conformal Cubed Sphere: 3D View")
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(size = (750, 750), figure_padding = figure_padding)

    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title, axis_kwargs_3D...)
    hide_decorations && hidedecorations!(ax3D)

    wireframe!(ax3D, X, Y, Z, color = colors[1], alpha = alphas[1])

    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

    for (i, R) in enumerate(rotations)

        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)

        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end

        wireframe!(ax3D, X′, Y′, Z′, color = colors[i + 1], alpha = alphas[i + 1])

    end

    return fig
end

function visualize_conformal_cubed_sphere_3D(X, Y, Z, file_name;
                                             figure_padding = (25, 0, 0, 50), # (left, right, bottom, top)
                                             title = "Conformal Cubed Sphere: 3D View",
                                             output_directory = "visualize_conformal_cubed_sphere_grid")
    hide_decorations = false

    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

    fig = Figure(size = (750, 750), figure_padding = figure_padding)

    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title, axis_kwargs_3D...)
    hide_decorations && hidedecorations!(ax3D)

    wireframe!(ax3D, X, Y, Z, color = colors[1], alpha = alphas[1])

    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

    for (i, R) in enumerate(rotations)
        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)

        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end

        wireframe!(ax3D, X′, Y′, Z′, color = colors[i + 1], alpha = alphas[i + 1])
    end

    save(output_directory * "/" * file_name, fig)
end

function visualize_conformal_cubed_sphere_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, colors,
                                                alphas;
                                                title_2D = "Conformal Cubed Sphere: 2D Projection",
                                                title_3D = "Conformal Cubed Sphere: 3D View")
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(size = (1500, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title_2D, axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title_3D,
                 axis_kwargs_3D...)

    for ax in [ax2D, ax3D]
        hide_decorations && hidedecorations!(ax)
        wireframe!(ax, X, Y, Z, color = colors[1], alpha = alphas[1])
    end

    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

    for (i, R) in enumerate(rotations)

        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)

        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end

        wireframe!(ax2D, X′, Y′, Z′, color = colors[i+1])
        wireframe!(ax3D, X′, Y′, Z′, color = colors[i+1], alpha = alphas[i+1])

    end

    colgap!(fig.layout, 60)

    return fig
end

function visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, file_name;
                                                title_2D = "Conformal Cubed Sphere: 2D Projection",
                                                title_3D = "Conformal Cubed Sphere: 3D View",
                                                output_directory = "visualize_conformal_cubed_sphere_grid")
    hide_decorations = false

    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

    fig = Figure(size = (1500, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title_2D, axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title_3D,
                 axis_kwargs_3D...)

    for ax in [ax2D, ax3D]
        hide_decorations && hidedecorations!(ax)
        wireframe!(ax, X, Y, Z, color = colors[1], alpha = alphas[1])
    end

    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

    for (i, R) in enumerate(rotations)
        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)

        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end

        wireframe!(ax2D, X′, Y′, Z′, color = colors[i+1])
        wireframe!(ax3D, X′, Y′, Z′, color = colors[i+1], alpha = alphas[i+1])
    end

    colgap!(fig.layout, 60)

    save(output_directory * "/" * file_name, fig)
end

function conformal_cubed_sphere_2D_3D_visualization_example()
    Nx, Ny = 16, 16
    hide_decorations = false
    color = :blue
    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]
    output_directory = "visualize_conformal_cubed_sphere_grid"

    fig = visualize_conformal_cubed_sphere_panel_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, color)
    save(output_directory * "/conformal_cubed_sphere_panel_2D.pdf", fig)

    fig = visualize_conformal_cubed_sphere_panel_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, color)
    save(output_directory * "/conformal_cubed_sphere_panel_3D.pdf", fig)

    fig = visualize_conformal_cubed_sphere_panel_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, color)
    save(output_directory * "/conformal_cubed_sphere_panel_2D_3D.pdf", fig)

    fig = visualize_conformal_cubed_sphere_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, colors)
    save(output_directory * "/conformal_cubed_sphere_2D.pdf", fig)

    fig = visualize_conformal_cubed_sphere_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, colors, alphas;
                                              figure_padding = (25, 0, 35, 50))
    save(output_directory * "/conformal_cubed_sphere_3D.pdf", fig)

    fig = visualize_conformal_cubed_sphere_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, colors,
                                                 alphas)
    save(output_directory * "/conformal_cubed_sphere_2D_3D.pdf", fig)
end

conformal_cubed_sphere_2D_3D_visualization_example()

function visualize_optimized_non_uniform_conformal_cubed_sphere(Nx, Ny; spacing = GeometricSpacing(), optimized = false)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    visualize_conformal_cubed_sphere_2D(X, Y, Z, "conformal_cubed_sphere_2D.pdf")
    visualize_conformal_cubed_sphere_3D(X, Y, Z, "conformal_cubed_sphere_3D.pdf";
                                        figure_padding = (25, 0, 35, 50))
    visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, "conformal_cubed_sphere_2D_3D.pdf")
    cell_areas = compute_cell_areas(X, Y, Z)
    reference_minimum_cell_area = minimum(cell_areas)

    if optimized
        x, y, X, Y, Z = optimized_non_uniform_conformal_cubed_sphere_coordinates(Nx, Ny, spacing)
    else
        x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny; spacing)
    end

    if optimized
        file_name_suffix = "_optimized"
        title_2D = ("Conformal Cubed Sphere with $(string(typeof(spacing))) Spacing\nand Optimized with EKI: "
                    * "2D Projection")
        title_3D = ("Conformal Cubed Sphere with $(string(typeof(spacing))) Spacing\nand Optimized with EKI: "
                    * "3D View")
    else
        file_name_suffix = "_unoptimized"
        title_2D = "Conformal Cubed Sphere with $(string(typeof(spacing))) Spacing:\n2D Projection"
        title_3D = "Conformal Cubed Sphere with $(string(typeof(spacing))) Spacing:\n3D View"
    end

    visualize_conformal_cubed_sphere_2D(
    X, Y, Z, "non_uniform_conformal_cubed_sphere_2D_" * lowercase(string(typeof(spacing))) * file_name_suffix * ".pdf"; title = title_2D)
    visualize_conformal_cubed_sphere_3D(
    X, Y, Z, "non_uniform_conformal_cubed_sphere_3D_" * lowercase(string(typeof(spacing))) * file_name_suffix * ".pdf"; title = title_3D)
    visualize_conformal_cubed_sphere_2D_3D(
    X, Y, Z, "non_uniform_conformal_cubed_sphere_2D_3D_" * lowercase(string(typeof(spacing))) * file_name_suffix * ".pdf";
    title_2D = title_2D, title_3D = title_3D)

    cell_areas = compute_cell_areas(X, Y, Z)
    minimum_cell_area = minimum(cell_areas)

    spacing_type = optimized ? "optimized" : "non-optimized"
    spacing_type *= " $(string(typeof(spacing)))"
    print("The normalized minimum cell width of the non-uniform conformal cubed sphere for $spacing_type " *
          "is $(sqrt(minimum_cell_area/reference_minimum_cell_area))\n")
end

for optimized in [false, true]
    for spacing in [GeometricSpacing(), ExponentialSpacing()]
        spacing_type = optimized ? "optimized" : "non-optimized"
        spacing_type *= " $(string(typeof(spacing)))"
        @info "Visualizing non-uniform conformal cubed sphere for $spacing_type"
        N = 16
        Nx, Ny = N + 1, N + 1
        visualize_optimized_non_uniform_conformal_cubed_sphere(Nx, Ny; spacing, optimized)
    end
end

specify_title(spacing) = "Conformal Cubed Sphere with " * string(typeof(spacing))

function specify_file_name_suffixes(spacing, optimized)
    file_name_suffix_1 = "_" * string(typeof(spacing))[1:end-7]
    file_name_suffix_2 = optimized ? "_Optimized" : "_Unoptimized"
    return file_name_suffix_1, file_name_suffix_2
end

function minimum_cell_width_variation_with_resolution(spacing, optimized;
                                                      output_directory = "visualize_conformal_cubed_sphere_grid")
    resolutions = 100:50:1000
    normalized_minimum_cell_widths = zeros(length(resolutions))

    minimum_reference_cell_area = 0

    for (i, resolution) in enumerate(resolutions)
        Nx, Ny = resolution + 1, resolution + 1
        @info "  Computing minimum cell width for conformal cubed sphere grid with resolution Nx = Ny = $resolution"
        x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
        cell_areas = compute_cell_areas(X, Y, Z)
        minimum_reference_cell_area = minimum(cell_areas)
        if optimized
            x, y, X, Y, Z = optimized_non_uniform_conformal_cubed_sphere_coordinates(Nx, Ny, spacing)
        else
            x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny; spacing)
        end
        cell_areas = compute_cell_areas(X, Y, Z)
        minimum_cell_area = minimum(cell_areas)
        normalized_minimum_cell_widths[i] = sqrt(minimum_cell_area/minimum_reference_cell_area)
    end
    file_name_suffix_1, file_name_suffix_2 = specify_file_name_suffixes(spacing, optimized)
    file_name = "MinimumCellWidthVersusResolution" * file_name_suffix_1 * file_name_suffix_2
    write_output_to_file_1D(output_directory, resolutions, normalized_minimum_cell_widths, file_name)
end

compute_minimum_cell_width_variation_with_resolution = true

if compute_minimum_cell_width_variation_with_resolution
    for spacing in [GeometricSpacing(), ExponentialSpacing()]
        for optimized in [false, true]
            spacing_type = optimized ? "optimized" : "non-optimized"
            spacing_type *= " $(string(typeof(spacing)))"
            @info "Computing minimum cell width variation with resolution for $spacing_type"
            minimum_cell_width_variation_with_resolution(spacing, optimized)
        end
    end
end

plot_minimum_cell_width_variation_with_resolution = true

if plot_minimum_cell_width_variation_with_resolution
    plot_size = (750, 750)
    plot_type = "scatter_line_plot"
    xlabel = "Number of cells in each direction of a cubed sphere panel"
    ylabel = "Normalized minimum cell width"
    axis_kwargs = (xlabel = "Number of cells in each direction of a cubed sphere panel",
                   ylabel = "Normalized minimum cell width", xlabelsize = 22.5, ylabelsize = 22.5,
                   xticklabelsize = 17.5, yticklabelsize = 17.5, xlabelpadding = 10, ylabelpadding = 10, aspect = 1.0,
                   titlesize = 27.5, titlegap = 15, titlefont = :bold)
    plot_kwargs = (linewidth = 2, linecolor = :black, marker = :rect, markersize = 15)
    plot_kwargs_2 = (linewidth = 2, linecolors = [:red, :blue], markers = [:rect, :rect], markersize = 15,
                     labels = ["Unoptimized", "Optimized with EKI"])

    for spacing in [GeometricSpacing(), ExponentialSpacing()]
        title = specify_title(spacing) * ":\nNormalized Minimum Cell Width versus Resolution"

        resolutions_2 = Vector{Float64}()
        normalized_minimum_cell_widths_optimized = Vector{Float64}()
        normalized_minimum_cell_widths_unoptimized = Vector{Float64}()

        for optimized in [false, true]
            spacing_type = optimized ? "optimized" : "non-optimized"
            spacing_type *= " $(string(typeof(spacing)))"
            @info "Plotting minimum cell width variation with resolution for $spacing_type"

            file_name_suffix_1, file_name_suffix_2 = specify_file_name_suffixes(spacing, optimized)

            file_name = "MinimumCellWidthVersusResolution" * file_name_suffix_1 * file_name_suffix_2

            resolutions, normalized_minimum_cell_widths = read_output_from_file_1D(output_directory, file_name * ".curve")
            create_single_line_or_scatter_plot(plot_size, plot_type, resolutions[2:end],
                                               normalized_minimum_cell_widths[2:end], axis_kwargs, title, plot_kwargs,
                                               file_name; figure_padding = (-10, 20, 0, 0), format=".pdf")

            resolutions_2 = resolutions

            if optimized
                normalized_minimum_cell_widths_optimized = normalized_minimum_cell_widths
            else
                normalized_minimum_cell_widths_unoptimized = normalized_minimum_cell_widths
            end
        end
    end
end

using GLMakie
using Elliptic
using CubedSphere

"""
    visualize_conformal_mapping(w; x_min, y_min, x_max, y_max, n_lines, n_samples, filepath="mapping.png")

Visualize a complex mapping ``w(z)`` from the rectangle `x ∈ [x_min, x_max]`, `y ∈ [y_min, y_max]` using `n_lines`
equally spaced lines in ``x`` and ``y`` each evaluated at `n_samples` sample points.
"""
function visualize_conformal_mapping(w; x_min, y_min, x_max, y_max, n_lines, n_samples, filepath="mapping.png")

    z₁ = [[x + im * y for x in range(x_min, x_max, length=n_samples)] for y in range(y_min, y_max, length=n_lines)]
    z₂ = [[x + im * y for y in range(y_min, y_max, length=n_samples)] for x in range(x_min, x_max, length=n_lines)]

    zs = cat(z₁, z₂, dims=1)
    ws = [w.(z) for z in zs]

    fig = Figure(size=(1200, 1200), fontsize=30)
    ax = Axis(fig[1, 1];
              xlabel = "Re(w)",
              ylabel = "Im(w)",
              aspect = DataAspect())

    [lines!(ax, real(z), imag(z), color = :grey) for z in zs]
    [lines!(ax, real(w), imag(w), color = :blue) for w in ws]

    save(filepath, fig)

    display(fig)
end

# maps square (±Ke, ±Ke) to circle of unit radius
w(z) = (1 - cn(z, 1/2)) / sn(z, 1/2)

Ke = Elliptic.F(π/2, 1/2) # ≈ 1.854

visualize_conformal_mapping(w;
                            x_min = -Ke,
                            y_min = -Ke,
                            x_max =  Ke,
                            y_max =  Ke,
                            n_lines = 21,
                            n_samples = 1000,
                            filepath = "square_to_disk_conformal_mapping.png")

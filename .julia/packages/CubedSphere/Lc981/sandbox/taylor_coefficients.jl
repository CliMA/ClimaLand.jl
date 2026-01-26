# Compute Rančić coefficients with the method described in
# [Rancic-etal-1996](@citet)

using GLMakie
using Printf

include(joinpath(@__DIR__, "../examples/compute_taylor_coefficients.jl"))

function plot_transformation(A, r, Nφ; Lφ=π/2)
    dφ = Lφ / Nφ
    φ = range(-Lφ/2 + dφ/2, stop=Lφ/2 - dφ/2, length=Nφ)

    z = @. r * cis(φ)
    Z = @. z^4

    W  = zeros(eltype(z), size(z))
    W̃′ = zeros(eltype(z), size(z))

    for k = length(A):-1:1
        @. W  += A[k] * z^(4k)
        @. W̃′ += A[k] * (1 - z)^(4k)
    end

    w  = @. cbrt(W)
    w′ = @. (1 - w) / (1 + w/2)
    W′ = @. w′^3

    w̃′ = @. cbrt′(W̃′)
    w̃  = @. (1 - w̃′) / (1 + w̃′/2)
    W̃  = @. w̃^3

    fig = Figure(size=(1200, 1800), fontsize=30)

    axz  = Axis(fig[1, 1], title="z")
    axZ  = Axis(fig[1, 2], title="Z")
    axw  = Axis(fig[2, 1], title="w")
    axW  = Axis(fig[2, 2], title="W")
    axwp = Axis(fig[3, 1], title="w′")
    axWp = Axis(fig[3, 2], title="W′")

    lim = 2
    for ax in (axw, axW)
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)
    end

    lim = 1
    for ax in (axwp, axWp)
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)
    end

    lim = 1.5
    for ax in (axz, axZ)
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)
    end

    scatter!(axz, real.(z), imag.(z))
    scatter!(axZ, real.(Z), imag.(Z))

    scatter!(axw, real.(w), imag.(w), color=(:black, 0.5), label=L"w")
    scatter!(axw, real.(w̃), imag.(w̃), color=(:orange, 0.8), label=L"w̃")

    scatter!(axW, real.(W), imag.(W), color=(:black, 0.5), label=L"W")
    scatter!(axW, real.(W̃), imag.(W̃), color=(:orange, 0.8), label=L"W̃")

    scatter!(axwp, real.(w′), imag.(w′), color=(:black, 0.5), label=L"w′")
    scatter!(axwp, real.(w̃′), imag.(w̃′), color=(:orange, 0.8), label=L"w̃′")

    scatter!(axWp, real.(W′), imag.(W′), color=(:black, 0.5), label=L"W′")
    scatter!(axWp, real.(W̃′), imag.(W̃′), color=(:orange, 0.8), label=L"W̃′")

    for ax in [axw, axW, axwp, axWp]
        axislegend(ax)
    end

    save("consistency.png", fig)

    return fig
end

r = 1 - 1e-6

Nφ = find_N(r; decimals=10)

maximum_coefficients = 128

Ncoefficients = Int(Nφ/2) - 2 > maximum_coefficients ? maximum_coefficients : Int(Nφ/2) - 2

N_iterations = 30

A, B = find_taylor_coefficients(r; maximum_coefficients, N_iterations)

@info "After $N_iterations iterations we have:"

for (k, Aₖ) in enumerate(A[1:30])
    @printf("k = %2i, A ≈ %+.14f, A_Rancic = %+.14f, |A - A_Rancic| = %.2e \n", k, real(Aₖ), CubedSphere.A_Rancic[k+1], abs(CubedSphere.A_Rancic[k+1] - real(Aₖ)))
end

fig = plot_transformation(A, r, Nφ)
fig

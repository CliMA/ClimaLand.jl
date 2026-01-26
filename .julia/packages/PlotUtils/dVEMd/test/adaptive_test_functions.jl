using Plots

const test_funcs = [
    (x -> 2, 0, 1),
    (x -> x, -1, 1),
    (x -> 5x, -1, 1),
    (x -> 1 / x, 0, 1),
    (x -> 1 / x, -1, 1),
    (x -> tan(x), 0, 4),
    (x -> -1 / x, 0, 1),
    (x -> 1 / abs(x), -1, 1),
    (x -> 1e6x, -1, 1),
    (x -> 1e50x, -1, 1),
    (x -> log(1 + sin(cos(x))), -6, 6),
    (x -> sin(x^3) + cos(x^3), 0, 6.28),
    (x -> sin(x), -5, 200),
    (x -> sin(1 / x), -2, 2),
    (x -> sin(1 / x), 0, 2),
    (x -> sin(x^4), -4, 4),
    (x -> sin(300x), -4, 4),
    (x -> 1 + x^2 + 0.0125log(abs(1 - 3(x - 1))), -2, 2),
    (x -> sin(exp(x)), -6, 6),
    (x -> 1 / sin(x), -10, 10),
    (x -> sin(x) / x, -6, 6),
    (x -> tan(x^3 - x + 1) + 1 / (x + 3exp(x)), -2, 2),
]

main() = begin
    plots = [plot(func...) for func ∈ test_funcs]

    pitchfork = plot(x -> sqrt(x), 0, 10)
    plot!(pitchfork, x -> -sqrt(x), 0, 10)
    plot!(pitchfork, x -> 0, -10, 0)
    plot!(pitchfork, x -> 0, 0, 10, linestyle = :dash)
    push!(plots, pitchfork)

    np = length(plots)
    m = ceil(Int, √(np)) + 1
    n, r = divrem(np, m)
    r == 0 || (n += 1)
    append!(plots, [plot() for _ ∈ 1:(m * n - np)])

    @assert length(plots) == m * n

    png(plot(plots...; layout = (m, n), size = (m * 600, n * 400)), "grid")
    return
end

main()

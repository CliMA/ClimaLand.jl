using GaussQuadrature
using PyPlot

#x = collect(range(0; stop=1, length=201))
x = range(0; stop=1, length=201)
const n = 5

a, b = shifted_legendre_coefs(Float64, n)
p = orthonormal_poly(x, a, b)
grid(true)

figure(1)
plot(x, p)
const r = 0
α, β = logweight_coefs(Float64, n+1, r)
q = orthonormal_poly(x, α[1:n], β)
title("Orthonormal shifted Legendre polynomials")

figure(2)
plot(x, q)
grid(true)
s = "Orthonormal polymials with weight \$\\log x^{-1}\$"
title(latexstring(s))

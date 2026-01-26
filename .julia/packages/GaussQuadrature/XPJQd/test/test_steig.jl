using GaussQuadrature
using Printf

T = Float64

println("\nFloating point data type is ", T)
@printf("\teps = %0.2e\n\n", eps(T))

n = 10
d = rand(T, n)
e = rand(T, n)
z = zero(d)

A = SymTridiagonal(d, e[1:end-1])
D, V = eig(A)

steig!(d, e, z, 30)
idx = sortperm(d)
d = d[idx]
z = z[idx]

@printf("\nEigenvalues:\n\n")
@printf("%20s  %20s  %12s\n\n", "eig", "steig!", "difference")
for i=1:length(d)
    @printf("%20.15f  %20.15f  %12.2e\n", d[i], D[i], d[i]-D[i])
end

@printf("\nFirst component of eigenvectors:\n\n")
@printf("%20s  %20s  %12s\n\n", "eig", "steig!", "difference")
for i=1:length(d)
    @printf("%20.15f  %20.15f  %12.2e\n", z[i], V[1,i], 
            abs(z[i])-abs(V[1,i]))
end

# Orthogonality tests

using GaussQuadrature
using Printf

T = BigFloat
#T = Float64
println("\nFloating point type data type is ", T)
println("\teps = ", eps(T))

xlegendre(n, alpha, beta, endpt) = legendre(T, n, endpt)
xchebyshev_first(n, alpha, beta, endpt) = chebyshev(T, n, 1, endpt)
xchebyshev_second(n, alpha, beta, endpt) = chebyshev(T, n, 2, endpt)
xhermite(n, alpha, beta, endpt) = hermite(T, n)
xlaguerre(n, alpha, beta, endpt) = laguerre(n, alpha, endpt)
xlogweight(n, alpha, beta , endpt) = logweight(T, n, 
                                     round(Int64,alpha), endpt)

xlegendre_coefs(n, alpha, beta) = legendre_coefs(T, n)
xchebyshev_first_coefs(n, alpha, beta) = chebyshev_coefs(T, n, 1)
xchebyshev_second_coefs(n, alpha, beta) = chebyshev_coefs(T, n, 2)
xhermite_coefs(n, alpha, beta) = hermite_coefs(T, n)
xlaguerre_coefs(n, alpha, beta) = laguerre_coefs(n, alpha)
xlogweight_coefs(n, alpha, beta) = logweight_coefs(T, n, 
                                   round(Int64,alpha))

rule = [xlegendre, xchebyshev_first, xchebyshev_second, 
        xhermite, jacobi, xlaguerre, xlogweight]

coefs = [ xlegendre_coefs, xchebyshev_first_coefs, 
          xchebyshev_second_coefs, xhermite_coefs, jacobi_coefs, 
          xlaguerre_coefs, xlogweight_coefs ]

name = ["Legendre", "Chebyshev (first)", "Chebyshev (second)",
         "Hermite", "Jacobi", "Laguerre", "Log weight" ]

alpha = one(T)
beta  = one(T)/2
xargs = [(),(1),(2), (), (alpha, beta), (alpha)]

function discrepancy(n, dop, w, p)
    d = 0.0
    for k = 0:dop
        for j = 0:min(k-1, dop-k)
            s = 0.0
            for l = 1:n
                s += w[l] * p[l,j+1] * p[l,k+1]
            end
            d = max(d, abs(s))
        end 
        if 2*k <= dop
            s = 0.0
            for l = 1:n
                s += w[l] * p[l,k+1]^2
            end
            d = max(d, abs(s-1.0))
        end
    end 
    return d
end

function dop(n, endpt)
    if endpt == neither
        return 2*n - 1
    elseif endpt == left || endpt == right
        return 2*n - 2
    else
        return 2*n - 3
    end
end

function test_rule(descr, nmin, nmax, rule, coefs, name, endpt)
    println("\nTesting ", descr, " rules with ", nmin, " to ", 
            nmax, " points")
    for case in zip(rule, coefs, name)
        rulefunc, coeffunc, rulename = case
        maxd = 0.0
        for n = nmin:nmax
            x, w = rulefunc(n, alpha, beta, endpt)
            a, b = coeffunc(2n, alpha, beta)
            p = orthonormal_poly(x, a, b)
            d = discrepancy(n, dop(n, endpt), w, p)
            maxd = max(d, maxd)
        end
        @printf("\t%s: %12.4e\n", rpad(rulename, 20, ' '), maxd)
    end
end

test_rule("plain Gauss", 1, 5, rule, coefs, name, neither)

filter!(x->(x!=xhermite), rule)
filter!(x->(x!=xhermite_coefs), coefs)
filter!(x->(x!="Hermite"), name)

test_rule("left Radau", 1, 5, rule, coefs, name, left)

filter!(x->(x!=xlaguerre), rule)
filter!(x->(x!=xlaguerre_coefs), coefs)
filter!(x->(x!="Laguerre"), name)

test_rule("right Radau", 1, 5, rule, coefs, name, right)
test_rule("Lobatto", 2, 5, rule, coefs, name, both)

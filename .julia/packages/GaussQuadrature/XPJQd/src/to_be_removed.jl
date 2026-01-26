export legendre_coeff, chebyshev_coeff, jacobi_coeff
export laguerre_coeff, hermite_coeff

function depwarn(old, new)
    msg = """
    $old will be removed in GaussQuadrature 0.4.
    Use $new instead."""
    warn(msg)
end

function legendre_coeff{T<:AbstractFloat}(::Type{T},
                       n::Integer, endpt::EndPt)
    depwarn("legendre_coeff", "legendre_coefs")
    μ0 = convert(T, 2.0)
    a = zeros(T, n)
    b = zeros(T, n)
    for i = 1:n
        b[i] = i / sqrt(convert(T, 4*i^2-1))
    end
    return a, b, μ0
end

function chebyshev_coeff{T<:AbstractFloat}(::Type{T},
                        n::Integer, kind::Integer, endpt::EndPt)
    depwarn("chebyshev_coeff", "chebyshev_coefs")
    muzero = convert(T, pi)
    half = convert(T, 0.5)
    a = zeros(T, n)
    b = fill(half, n)
    if kind == 1
        b[1] = sqrt(half)
    elseif kind == 2
        muzero /= 2
    else
        error("Unsupported value for kind")
    end
    return a, b, muzero
end

function jacobi_coeff{T<:AbstractFloat}(n::Integer, alpha::T,
                                        beta::T, endpt::EndPt)
    depwarn("jacobi_coeff", "jacobi_coefs")
    ab = alpha + beta
    i = 2
    abi = ab + 2
    μ0 = 2^(ab+1) * exp(
             lgamma(alpha+1) + lgamma(beta+1) - lgamma(abi) )
    a = zeros(T, n)
    b = zeros(T, n)
    a[1] = ( beta - alpha ) / abi
    b[1] = sqrt( 4*(alpha+1)*(beta+1) / ( (ab+3)*abi*abi ) )
    a2b2 = beta*beta - alpha*alpha
    for i = 2:n
        abi = ab + 2i
        a[i] = a2b2 / ( (abi-2)*abi )
        b[i] = sqrt( 4i*(alpha+i)*(beta+i)*(ab+i) /
                     ( (abi*abi-1)*abi*abi ) )
    end
    return a, b, μ0
end

function laguerre_coeff{T<:AbstractFloat}(n::Integer, alpha::T,
                                          endpt::EndPt)
    depwarn("laguerre_coeff", "laguerre_coefs")
    @assert endpt in [neither, left]
    μ0 = gamma(alpha+1)
    a = zeros(T, n)
    b = zeros(T, n)
    for i = 1:n
        a[i] = 2i - 1 + alpha
        b[i] = sqrt( i*(alpha+i) )
    end
    return a, b, μ0
end

function hermite_coeff{T<:AbstractFloat}(::Type{T}, n::Integer)
    depwarn("hermite_coeff", "hermite_coefs")
    μ0 = sqrt(convert(T, pi))
    a = zeros(T, n)
    b = zeros(T, n)
    for i = 1:n
        iT = convert(T, i)
        b[i] = sqrt(iT/2)
    end
    return a, b, μ0
end

function custom_gauss_rule{T<:AbstractFloat}(lo::T, hi::T,
         a::Array{T,1}, b::Array{T,1}, μ0::T, endpt::EndPt,
         maxits::Integer=maxiterations[T])
    #
    # On entry:
    #
    # a, b hold the coefficients (as given, for instance, by
    # legendre_coeff) in the three-term recurrence relation
    # for the orthonormal polynomials p_0, p_1, p_2, ... , that is,
    #
    #    b[j] p (x) = (x-a[j]) p   (x) - b[j-1] p   (x).
    #          j                j-1              j-2
    #
    # μ0 holds the zeroth moment of the weight function, that is
    #
    #          / hi
    #         |
    #    μ0 = | w(x) dx.
    #         |
    #         / lo
    #
    # On return:
    #
    # x, w hold the points and weights.
    #
    warn("This method will be removed in GaussQuadrature 0.4")
    n = length(a)
    @assert length(b) == n
    if endpt == left
        if n == 1
            a[1] = lo
        else
            a[n] = solve(n, lo, a, b) * b[n-1]^2 + lo
        end
    elseif endpt == right
        if n == 1
            a[1] = hi
        else
            a[n] = solve(n, hi, a, b) * b[n-1]^2 + hi
        end
    elseif endpt == both
        if n == 1
            error("Must have at least two points for both ends.")
        end
        g = solve(n, lo, a, b)
        t1 = ( hi - lo ) / ( g - solve(n, hi, a, b) )
        b[n-1] = sqrt(t1)
        a[n] = lo + g * t1
    end
    w = zero(a)
    steig!(a, b, w, maxits)
    for i = 1:n
        w[i] = μ0 * w[i]^2
    end
    idx = sortperm(a)
    return a[idx], w[idx]
end

function steig!{T<:AbstractFloat}(d::Array{T,1}, e::Array{T,1}, 
                                  z::Array{T,1}, maxits::Integer)
    #
    # Finds the eigenvalues and first components of the normalised
    # eigenvectors of a symmetric tridiagonal matrix by the implicit
    # QL method.
    #
    # d[i]   On entry, holds the ith diagonal entry of the matrix. 
    #        On exit, holds the ith eigenvalue.
    #
    # e[i]   On entry, holds the [i+1,i] entry of the matrix for
    #        i = 1, 2, ..., n-1.  (The value of e[n] is not used.)
    #        On exit, e is overwritten.
    #
    # z[i]   On exit, holds the first component of the ith normalised
    #        eigenvector associated with d[i].
    #
    # maxits The maximum number of QL iterations.
    #
    # Martin and Wilkinson, Numer. Math. 12: 377-383 (1968).
    # Dubrulle, Numer. Math. 15: 450 (1970).
    # Handbook for Automatic Computation, Vol ii, Linear Algebra, 
    #        pp. 241-248, 1971.
    #
    # This is a modified version of the Eispack routine imtql2.
    #
    n = length(z)
    z[1] = 1
    z[2:n] = 0
    e[n] = 0

    if n == 1 # Nothing to do for a 1x1 matrix.
        return
    end
    for l = 1:n
        for j = 1:maxits
            # Look for small off-diagonal elements.
            m = n
            for i = l:n-1
                if abs(e[i]) <= eps(T) * ( abs(d[i]) + abs(d[i+1]) )
                    m = i
                    break   
                end
            end
            p = d[l]
            if m == l
                continue
            end
            if j == maxits
                msg = @sprintf("No convergence after %d iterations", j)
                msg *= " (try increasing maxits)"
                error(msg)
            end
            # Form shift
            g = ( d[l+1] - p ) / ( 2 * e[l] )
            r = hypot(g, one(T))
            g = d[m] - p + e[l] / ( g + copysign(r, g) )
            s = one(T)
            c = one(T)
            p = zero(T)
            for i = m-1:-1:l
                f = s * e[i]
                b = c * e[i]
                if abs(f) <  abs(g)
                    s = f / g
                    r = hypot(s, one(T))
                    e[i+1] = g * r
                    c = one(T) / r
                    s *= c
                else
                    c = g / f
                    r = hypot(c, one(T))
                    e[i+1] = f * r
                    s = one(T) / r
                    c *= s
                end 
                g = d[i+1] - p
                r = ( d[i] - g ) * s + 2 * c * b
                p = s * r
                d[i+1] = g + p
                g = c * r - b
                # Form first component of vector.
                f = z[i+1]
                z[i+1] = s * z[i] + c * f
                z[i]   = c * z[i] - s * f
            end # loop over i
            d[l] -= p
            e[l] = g
            e[m] = zero(T)
        end # loop over j
    end # loop over l
end

function orthonormal_poly{T<:AbstractFloat}(x::Array{T,1},
                         a::Array{T,1}, b::Array{T,1}, μ0::T)
    # p[i,j] = value at x[i] of orthonormal polynomial of degree j-1.
    m = length(x)
    n = length(a)
    p = zeros(T, m, n+1)
    c = one(T) / sqrt(μ0)
    rb = one(T) / b[1]
    for i = 1:m
        p[i,1] = c
        p[i,2] = rb * ( x[i] - a[1] ) * c
    end
    for j = 2:n
       rb = one(T) / b[j]
       for i = 1:m
           p[i,j+1] = rb * ( (x[i]-a[j]) * p[i,j]
                                - b[j-1] * p[i,j-1] )
       end
    end
    return p
end


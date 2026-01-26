function lanczosIterations(A, stepsize, αs, βs, vs, reorthogonalize)

    iter = length(αs)
    v1 = vs[iter]
    v = vs[iter + 1]
    β = βs[iter]

    for i = iter .+ (1:stepsize)
        w = A*v
        α = real(dot(w,v))

        axpy!(-α, v, w)
        axpy!(-β, v1, w)

        # Orthogonlize (FixMe! consider a partial reorthogonal strategy here)
        if reorthogonalize
            for j = 1:i - 1
                axpy!(-dot(vs[j], w), vs[j], w)
            end
        end

        v1 = v
        v, β = qr!(w)

        push!(αs, α)
        push!(βs, β)
        push!(vs, v)
    end
    return αs, βs, vs
end

function _teig(A,
    nvals = 1;
    maxiter = size(A, 1),
    initvec = convert(Vector{eltype(A)}, randn(size(A, 1))),
    tolconv = sqrt(eps(real(eltype(A)))),
    reorthogonalize = true,
    stepsize = max(1, div(nvals, 10)),
    debug = false)

    # LinAlg.chksquare(A)

    ### initialize αs, βs and vs ###
    # run a single iteration of the Lanczos procedure. Note that this step handles promotion Lanczos vectors as well. They will in in general be different from `initvec`.
    v1 = A*initvec
    v0 = convert(typeof(v1), initvec)
    v0, nm0 = qr!(v0)
    rmul!(v1, inv(nm0)) # because we don't assume that initvec is normalized
    α = real(dot(v1, v0))
    axpy!(-α, v0, v1)
    v1, β = qr!(v1)
    αs = [α]
    βs = [β]
    vs = typeof(v1)[v0, v1]

    ### initial iterations ###
    lanczosIterations(A, nvals - 1, αs, βs, vs, reorthogonalize)

    ### Iteration count ###
    iter = nvals

    # values0 = sort!(abs(eigvals(SymTridiagonal(αs, βs[1:end - 1]))), rev = true)
    # values1 = values0

    z = zeros(eltype(αs), 1)

    ### run itertions ###
    # continue until convergence or the maximum number of iterations has been reached.
    converged = UInt[]
    while iter <= maxiter
        lanczosIterations(A, stepsize, αs, βs, vs, reorthogonalize)
        # values, vectors = LinearAlgebra.EigenSelfAdjoint.eigQL!(SymTridiagonal(copy(αs), βs[1:end - 1]), [zeros(length(αs) - 1); ones(eltype(αs), 1)]')
        values, vectors = eigen(SymTridiagonal(αs, βs[1:end-1]))

        ### Use criteria (9.18) of Peter Arbenz' lecture notes ###
        βj = βs[end]
        converged = Int[]
        for i = 1:length(αs)
            if βj*abs(vectors[end,i]) < abs(values[i])*tolconv
                push!(converged, i)
            end
        end

        debug && @show converged

        # if length(converged) >= nvals
        #     break
        # end

        # For now we hardcode convergence for largest positive value. This should be a choice later on
        if all([length(αs) - i + 1 in converged for i = 1:nvals])
            break
        end

        iter += stepsize
        debug && @show iter
    end
    if length(converged) < nvals
        error("no convergence")
    end

    S = SymTridiagonal(αs, βs[1:end-1])
    values, vectors = eigen(S)
    return values[converged], vectors[:, converged], S, vs, βs[end]
end

teig(A,
    nvals = 1;
    maxiter = size(A, 1),
    initvec = convert(Vector{eltype(A)}, randn(size(A, 1))),
    tolconv = 1e-12,
    reorthogonalize = true,
    stepsize = max(1, div(nvals, 10)),
    debug = false) = _teig(A, nvals, maxiter = maxiter, initvec = initvec,
        tolconv = tolconv, reorthogonalize = reorthogonalize, stepsize = stepsize, debug = debug)

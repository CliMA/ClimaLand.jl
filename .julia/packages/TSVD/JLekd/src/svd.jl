function _biLanczosIterations!(A, stepsize, αs, βs, U, V, μs, νs, maxνs, maxμs, τ, reorth_in, tolreorth, debug)

    m, n = size(A)
    reorth_μ = reorth_in
    nReorth = 0
    nReorthVecs = 0

    T = eltype(eltype(U))
    Tr = real(T)

    iter = length(αs)

    u = U[iter + 1]
    v = V[iter]
    β = βs[iter]

    for j = iter .+ (1:stepsize)
        # The v step
        vOld = v
        ## apply operator
        v = A'u
        axpy!(T(-β), vOld, v)
        α = norm(v)

        ## update norm(A) estimate. FixMe! Use tighter bounds, see Larsen's thesis page 33
        τ = max(τ, eps(Tr) * (α + β))
        debug && @show τ

        ## run ω recurrence
        reorth_ν = false
        for i = 1:j - 1
            ν = βs[i]*μs[i + 1] + αs[i]*μs[i] - β*νs[i]
            ν = (ν + copysign(τ, ν))/α
            if abs(ν) > tolreorth
                reorth_ν |= true
            end
            νs[i] = ν
        end
        if j > 1
            push!(maxνs, maximum(abs, νs))
        end
        push!(νs, 1)

        ## reorthogonalize if necessary
        if reorth_ν || reorth_μ
            debug && println("Reorth v")
            for i in 1:j - 1
                axpy!(-dot(V[i], v), V[i], v)
                νs[i] = eps(Tr)
                nReorthVecs += 1
            end
            α = norm(v)
        end

        ## update the result vectors
        push!(αs, α)
        rmul!(v, inv(α))
        push!(V, v)

        # The u step
        uOld = u
        ## apply operator
        u = A*v
        axpy!(T(-α), uOld, u)
        β = norm(u)

        ## update norm(A) estimate. FixMe! Use tighter bounds, see Larsen's thesis page 33
        τ = max(τ, eps(Tr) * (α + β))
        debug && @show τ

        ## run ω recurrence
        reorth_μ = false
        for i = 1:j
            μ = αs[i]*νs[i] - α*μs[i]
            if i > 1
                μ += βs[i - 1]*νs[i-1]
            end
            μ = (μ + copysign(τ, μ))/β
            if abs(μ) > tolreorth
                reorth_μ |= true
            end
            μs[i] = μ
        end
        push!(maxμs, maximum(μs))
        push!(μs, 1)

        ## reorthogonalize if necessary
        if reorth_ν || reorth_μ
            debug && println("Reorth u")
            for i in 1:j
                axpy!(-dot(U[i], u), U[i], u)
                μs[i] = eps(Tr)
                nReorthVecs += 1
            end
            β = norm(u)
            nReorth += 1
        end

        ## update the result vectors
        push!(βs, β)
        rmul!(u, inv(β))
        push!(U, u)
    end

    return reorth_μ
end

function biLanczos(A,
    nvals = 1;
    maxiter = 1000,
    initvec = convert(Vector{float(eltype(A))}, randn(size(A,1))),
    tolconv = sqrt(eps(real(eltype(initvec)))),
    tolreorth = sqrt(eps(real(eltype(initvec)))),
    stepsize = max(1, div(nvals, 10)),
    debug = false)

    Tv = eltype(initvec)
    Tr = real(Tv)

    # I need to append βs with a zero at each iteration. Tt is much easier for type inference if it is a vector with the right element type
    z = zeros(Tr, 1)

    # initialize the αs, βs, U and V. Use result of first matvec to infer the correct types.
    # So the first iteration is run here, but slightly differently from the rest of the iterations
    nrmInit = norm(initvec)
    v = A'initvec
    rmul!(v, inv(nrmInit))
    α = norm(v)
    if isone(α)
        # We need v to change when scaled with α.
        # Otherwise the basis generation might break down.
        α = nextfloat(α)
    end
    rmul!(v, inv(α))
    V = fill(v, 1)
    αs = fill(α, 1)

    u = A*v
    uOld = similar(u)
    copyto!(uOld, initvec)
    rmul!(uOld, inv(nrmInit))
    axpy!(eltype(u)(-α), uOld, u)
    β = norm(u)
    rmul!(u, inv(β))
    U = typeof(u)[uOld, u]
    βs = fill(β, 1)

    # error estimate used in ω recurrence
    τ = eps(Tr)*(α + β)
    ν = 1 + τ/α
    μ = τ/β

    # Arrays for saving the estimates of the maximum angles between Lanczos vectors
    maxμs = Tr[]
    maxνs = Tr[]

    μs = Tr[μ, 1]
    νs = Tr[one(μ)]

    reorth_μ = _biLanczosIterations!(A, nvals - 1, αs, βs, U, V, μs, νs, maxμs, maxνs, τ, false, tolreorth, debug)

    # Iteration count
    iter = nvals

    # Save the estimates of the maximum angles between Lanczos vectors
    # append!(maxμs, maxμ)

    hasConv = false
    while iter <= maxiter

        reorth_μ = _biLanczosIterations!(A, stepsize, αs, βs, U, V, μs, νs, maxμs, maxνs, τ, reorth_μ, tolreorth, debug)
        iter += stepsize

        # This is more expensive than necessary because we only need the last components. However, LAPACK doesn't support this.
        UU, ss, VV = svd(Bidiagonal([αs;z], βs, :L))

        debug && @show βs[end]

        # Test for convergence. A Ritzvalue is considered converged if
        # either the last component of the corresponding vector is (relatively)
        # small or if the last component in βs is small (or both)
        if all(abs.(UU[end, 1:nvals])*βs[end] .< tolconv*ss[1:nvals]) &&
           all(abs.(VV[end, 1:nvals])*βs[end] .< tolconv*ss[1:nvals])
            hasConv = true
            break
        end

        τ = eps(eltype(ss))*ss[1]

        debug && @show iter
        debug && @show τ
    end
    if !hasConv
        error("no convergence")
    end

    # Form upper bidiagonal square matrix
    # m = length(U[1])
    # for j = 1:length(αs)
    #     # Calculate Givens rotation
    #     c, s, αs[j] = givensAlgorithm(αs[j], βs[j])

    #     # Update left vector
    #     # for i = 1:m
    #     #     uij       = U[j][i]
    #     #     uij1      = U[j+1][i]
    #     #     U[j][i]   = uij*c + uij1*s'
    #     #     U[j+1][i] = uij1*c - uij*s
    #     # end

    #     # Update bidiagonal matrix
    #     if j < length(αs)
    #         αj1 = αs[j + 1]
    #         αs[j + 1] = c*αj1
    #         βs[j] = s*αj1
    #     end
    # end

    return U, Bidiagonal(αs, βs[1:end-1], :L), V, maxμs, maxνs
end

function _tsvd(A,
    nvals = 1;
    maxiter = 1000,
    # The initial vector is critical in determining the output type.
    # We use the result of A*initvec to detemine the storage type for
    # the Lanczos vectors. Hence, the user would need to either
    # have an appropriate multiplication method defined or supply
    # an appropriate initial vector.
    initvec = convert(Vector{float(eltype(A))}, randn(size(A,1))),
    tolconv = sqrt(eps(real(eltype(initvec)))),
    tolreorth = sqrt(eps(real(eltype(initvec)))),
    stepsize = max(1, div(nvals, 10)),
    debug = false)

    U, B, V, maxμs, maxνs = biLanczos(A,
        nvals;
        maxiter = maxiter,
        initvec = initvec,
        tolconv = tolconv,
        tolreorth = tolreorth,
        stepsize = stepsize,
        debug = debug)

    # Calculate the bidiagonal SVD
    smU, sms, smV = svd(B)

    # Create matrices from the Vectors of vectors and update U and V
    mU = hcat(U[1:end-1])
    mV = hcat(V)

    # Adapt the matrix type of the (LAPACK) SVD matrices to the type of mU and mV
    # E.g. if the vectors have been stored on a GPU device or on a distributed
    # system, this step will ensure that the result ends up where A/initvec is stored
    mUall = mU*adapt(typeof(initvec), smU)
    mVall = mV*adapt(typeof(initvec), smV)
    msall = adapt(typeof(initvec), sms)

    return mUall[:,1:nvals],
        sms[1:nvals],
        mVall[:,1:nvals],
        B,
        maxμs,
        maxνs
end

"""
    tsvd(A, nvals = 1; [maxiter, initvec, tolconv, tolreorth, debug])

Computes the truncated singular value decomposition (TSVD) by Lanczos bidiagonalization of the operator `A`. The Lanczos vectors are partially orthogonalized as described in

R. M. Larsen, *Lanczos bidiagonalization with partial reorthogonalization*, Department of Computer Science, Aarhus University, Technical report, DAIMI PB-357, September 1998.



# Positional arguments:

- `A`: Anything that supports the in place update operations


    mul!(y::AbstractVector, A, x::AbstractVector, α::Number, β::Number)

and

    mul!(y::AbstractVector, A::Adjoint, x::AbstractVector, α::Number, β::Number)

corresponding to the operations `y := α*op(A)*x + β*y` where `op` can be either the identity or the conjugate transpose of `A`. If the `initvec` argument is not supplied then it is furthermore required that `A` supports `eltype` and `size`.

- `nvals`: The number of singular values and vectors to compute. Default is one (the largest).



# Keyword arguments:

- `maxiter`: The maximum number of iterations of the Lanczos bidiagonalization. Default is 1000, but usually much fewer iterations are needed.

- `initvec`: Initial `U` vector for the Lanczos procesdure. Default is a vector of Gaussian random variates. The `length` and `eltype` of the `initvec` will control the size and element types of the basis vectors in `U` and `V`.

- `tolconv`: Relative convergence criterion for the singular values. Default is `sqrt(eps(real(eltype(A))))`.

- `tolreorth`: Absolute tolerance for the inner product of the Lanczos vectors as measured by the ω recurrence. Default is `sqrt(eps(real(eltype(initvec))))`. `0.0` and `Inf` correspond to complete and no reorthogonalization respectively.

- `debug`: Boolean flag for printing debug information



# Output:

The output of the procesure it the truple tuple `(U, s, V)`

- `U`: `size(A, 1)` times `nvals` matrix of left singular vectors.
- `s`: Vector of length `nvals` of the singular values of `A`.
- `V`: `size(A, 2)` times `nvals` matrix of right singular vectors.



# Examples

```jldoctest
julia> A = matrixdepot("LPnetlib/lp_osa_30")
4350×104374 SparseArrays.SparseMatrixCSC{Float64, Int64} with 604488 stored entries:
⠙⠮⠷⠶⠽⠶⠽⠶⠮⠷⠮⠷⠶⠽⠶⠽⠶⠬⠷⠮⠷⠦⠽⠶⠽⠶⠽⠶⠮⠷⠮⠷⠶⠽⠶⠽⠶⠭⠷⠦

julia> U, s, V = tsvd(A, 5);

julia> round.(s, digits=7)
5-element Vector{Float64}:
 1365.8944098
 1033.2125634
  601.3524529
  554.107656
  506.0414587
```
"""
function tsvd(
    A,
    nvals = 1;
    maxiter = 1000,
    initvec = convert(Vector{float(eltype(A))}, randn(size(A,1))),
    tolconv = sqrt(eps(real(eltype(initvec)))),
    tolreorth = sqrt(eps(real(eltype(initvec)))),
    stepsize = max(1, div(nvals, 10)),
    debug = false)

    result = _tsvd(A, nvals, maxiter = maxiter, initvec = initvec, tolconv = tolconv,
            tolreorth = tolreorth, debug = debug)
    return (result[1], result[2], result[3])
end


### SVD by Lanczos on A'A

mutable struct AtA{T,S<:AbstractMatrix,V<:AbstractVecOrMat} <: AbstractMatrix{T}
    matrix::S
    vector::V
end
function AtA(A::AbstractMatrix, x::AbstractVecOrMat)
    y = A*x
    return AtA{eltype(y),typeof(A),typeof(y)}(A,y)
end

function size(A::AtA, i::Integer)
    if i < 1
        error("arraysize: dimension out of range")
    elseif i < 3
        return size(A.matrix, 2)
    else
        return 1
    end
end
size(A::AtA) = (size(A, 1), size(A, 2))

# Split Vector and Matrix to avoid ambiguity
function mul!(y::AbstractVector{T},
              A::AtA{T,S,V},
              x::AbstractVector{T},
              α::T = one(T),
              β::T = zero(T)) where {T<:Number,S,V}
    mul!(A.vector, A.matrix, x, one(T), zero(T))
    mul!(y, A.matrix', A.vector, α, β)
    return y
end
function mul!(y::AbstractMatrix{T},
              A::AtA{T,S,V},
              x::AbstractMatrix{T},
              α::T = one(T),
              β::T = zero(T)) where {T<:Number,S,V}
    mul!(A.vector, A.matrix, x, one(T), zero(T))
    mul!(y, A.matrix', A.vector, α, β)
    return y
end
(*)(A::AtA, x::AbstractVector) = mul!(similar(A.vector, size(x)), A, x)
(*)(A::AtA, x::AbstractMatrix) = mul!(similar(A.vector, size(x)), A, x)

function tsvd2(A,
    nvals = 1;
    maxiter = min(size(A)...),
    initvec = convert(Vector{eltype(A)}, randn(size(A,2))),
    tolconv = sqrt(eps(real(eltype(A)))),
    stepsize = max(1, div(nvals, 10)),
    debug = false)
    values, vectors, S, lanczosVecs = _teig(AtA(A, initvec), nvals, maxiter = maxiter,
        initvec = initvec, tolconv = tolconv, stepsize = stepsize, debug = debug)
    mV = hcat(lanczosVecs[1:end-1])*vectors
    return sqrt.(reverse(values)[1:nvals]), mV[:,end:-1:1][:,1:nvals]
end

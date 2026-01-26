#
"""
    BatchedDiagonalOperator(diag; update_func, update_func!, accepted_kwargs)

Represents a time-dependent elementwise scaling (diagonal-scaling) operation.
Acts on `AbstractArray`s of the same size as `diag`. The update function is
called by `update_coefficients!` and is assumed to have the following signature:

    update_func(diag::AbstractArray, u, p, t; <accepted kwarg fields>) -> [modifies diag]
"""
struct BatchedDiagonalOperator{T, D, F, F!} <: AbstractSciMLOperator{T}
    diag::D
    update_func::F
    update_func!::F!

    function BatchedDiagonalOperator(diag::AbstractArray, update_func, update_func!)
        new{
            eltype(diag),
            typeof(diag),
            typeof(update_func),
            typeof(update_func!)
        }(diag,
            update_func,
            update_func!)
    end
end

function DiagonalOperator(u::AbstractArray;
        update_func = nothing,
        update_func! = nothing,
        accepted_kwargs = nothing)
    update_func = preprocess_update_func(update_func, accepted_kwargs)
    update_func! = preprocess_update_func(update_func!, accepted_kwargs)

    BatchedDiagonalOperator(u, update_func, update_func!)
end

# traits
function Base.show(io::IO, L::BatchedDiagonalOperator)
    n, k = size(L.diag)
    print(io, "BatchedDiagonalOperator($n × $n, k = $k)")
end
Base.size(L::BatchedDiagonalOperator) = (N = size(L.diag, 1); (N, N))
Base.iszero(L::BatchedDiagonalOperator) = iszero(L.diag)
Base.transpose(L::BatchedDiagonalOperator) = L
Base.adjoint(L::BatchedDiagonalOperator) = conj(L)
function Base.conj(L::BatchedDiagonalOperator) # TODO - test this thoroughly
    update_func,
    update_func! = if isreal(L)
        L.update_func, L.update_func!
    else
        uf = L.update_func === nothing ? nothing :
             (
            L, u, p, t; kwargs...) -> conj(L.update_func(conj(L.diag),
            u,
            p,
            t;
            kwargs...))
        uf! = L.update_func! === nothing ? nothing :
              (
            L, u, p, t; kwargs...) -> begin
            L.update_func!(conj!(L.diag), u, p, t; kwargs...)
            conj!(L.diag)
        end
        uf, uf!
    end

    DiagonalOperator(conj(L.diag);
        update_func = update_func,
        update_func! = update_func!,
        accepted_kwargs = NoKwargFilter())
end

function Base.convert(::Type{AbstractMatrix}, L::BatchedDiagonalOperator)
    m, n = size(L)
    msg = """$L cannot be represented by an $m × $n AbstractMatrix"""
    throw(ArgumentError(msg))
end

LinearAlgebra.issymmetric(L::BatchedDiagonalOperator) = true
function LinearAlgebra.ishermitian(L::BatchedDiagonalOperator)
    if isreal(L)
        true
    else
        vec(L.diag) |> Diagonal |> ishermitian
    end
end
LinearAlgebra.isposdef(L::BatchedDiagonalOperator) = isposdef(Diagonal(vec(L.diag)))

function update_coefficients(L::BatchedDiagonalOperator, u, p, t; kwargs...)
    if !isnothingfunc(L.update_func)
        return @reset L.diag = L.update_func(L.diag, u, p, t; kwargs...)
    elseif !isnothingfunc(L.update_func!)
        L.update_func!(L.diag, u, p, t; kwargs...)
        return L
    end
end

function update_coefficients!(L::BatchedDiagonalOperator, u, p, t; kwargs...)
    if !isnothingfunc(L.update_func!)
        L.update_func!(L.diag, u, p, t; kwargs...)
    elseif !isnothingfunc(L.update_func)
        L.diag = L.update_func(L.diag, u, p, t; kwargs...)
    end
    nothing
end

getops(L::BatchedDiagonalOperator) = (L.diag,)

# Copy method to avoid aliasing
function Base.copy(L::BatchedDiagonalOperator)
    BatchedDiagonalOperator(
        copy(L.diag),
        L.update_func,
        L.update_func!
    )
end

function isconstant(L::BatchedDiagonalOperator)
    update_func_isconstant(L.update_func) & update_func_isconstant(L.update_func!)
end
islinear(::BatchedDiagonalOperator) = true
isconvertible(::BatchedDiagonalOperator) = false
has_adjoint(L::BatchedDiagonalOperator) = true
has_ldiv(L::BatchedDiagonalOperator) = all(x -> !iszero(x), L.diag)
has_ldiv!(L::BatchedDiagonalOperator) = has_ldiv(L)

# operator application
Base.:*(L::BatchedDiagonalOperator, u::AbstractVecOrMat) = L.diag .* u
Base.:\(L::BatchedDiagonalOperator, u::AbstractVecOrMat) = L.diag .\ u

function LinearAlgebra.mul!(v::AbstractVecOrMat,
        L::BatchedDiagonalOperator,
        u::AbstractVecOrMat)
    V = vec(v)
    U = vec(u)
    d = vec(L.diag)
    D = Diagonal(d)
    mul!(V, D, U)

    v
end

function LinearAlgebra.mul!(v::AbstractVecOrMat,
        L::BatchedDiagonalOperator,
        u::AbstractVecOrMat,
        α,
        β)
    V = vec(v)
    U = vec(u)
    d = vec(L.diag)
    D = Diagonal(d)
    mul!(V, D, U, α, β)

    v
end

function LinearAlgebra.ldiv!(v::AbstractVecOrMat,
        L::BatchedDiagonalOperator,
        u::AbstractVecOrMat)
    V = vec(v)
    U = vec(u)
    d = vec(L.diag)
    D = Diagonal(d)
    ldiv!(V, D, U)

    v
end

function LinearAlgebra.ldiv!(L::BatchedDiagonalOperator, u::AbstractVecOrMat)
    U = vec(u)
    d = vec(L.diag)
    D = Diagonal(d)
    ldiv!(D, U)
    u
end

function (L::BatchedDiagonalOperator)(v::AbstractVecOrMat, u, p, t; kwargs...)
    if !isconstant(L)
        L = update_coefficients(L, u, p, t; kwargs...)
    end
    L.diag .* v
end

function (L::BatchedDiagonalOperator)(
        w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t; kwargs...)
    update_coefficients!(L, u, p, t; kwargs...)
    w .= L.diag .* v
    return w
end

function (L::BatchedDiagonalOperator)(
        w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t, α, β; kwargs...)
    update_coefficients!(L, u, p, t; kwargs...)
    if β == 0
        w .= α .* (L.diag .* v)
    else
        w .= α .* (L.diag .* v) .+ β .* w
    end
    return w
end
#

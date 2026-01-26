using Test

function __rprint_diff(
    io::IO,
    x::T,
    y::T;
    pc,
    xname,
    yname,
) where {T <: NamedTuple}
    for pn in propertynames(x)
        pc_full = (pc..., ".", pn)
        xi = getproperty(x, pn)
        yi = getproperty(y, pn)
        __rprint_diff(io, xi, yi; pc = pc_full, xname, yname)
    end
end;

function __rprint_diff(io::IO, xi, yi; pc, xname, yname) # assume we can compute difference here
    if !(xi == yi)
        xs = xname * string(join(pc))
        ys = yname * string(join(pc))
        println(io, "==================== Difference found:")
        println(io, "maximum(abs.(Δ)) = $(maximum(abs.(xi .- yi)))")
        println(io, "maximum(abs.(xi)) = $(maximum(abs.(xi)))")
        println(io, "maximum(abs.(yi)) = $(maximum(abs.(yi)))")
        println(io, "extrema(xi) = $(extrema(xi))")
        println(io, "extrema(yi) = $(extrema(yi))")
        # println(io, "$xs: ", xi)
        # println(io, "$ys: ", yi)
        # println(io, "($xs .- $ys): ", (xi .- yi))
    end
    return nothing
end

"""
    rprint_diff(io::IO, ::T, ::T) where {T <: NamedTuple}
    rprint_diff(::T, ::T) where {T <: NamedTuple}

Recursively print differences in given `NamedTuple`.
"""
_rprint_diff(io::IO, x::T, y::T, xname, yname) where {T <: NamedTuple} =
    __rprint_diff(io, x, y; pc = (), xname, yname)
_rprint_diff(x::T, y::T, xname, yname) where {T <: NamedTuple} =
    _rprint_diff(stdout, x, y, xname, yname)

"""
    @rprint_diff(::T, ::T) where {T <: NamedTuple}

Recursively print differences in given `NamedTuple`.
"""
macro rprint_diff(x, y)
    return :(_rprint_diff(
        stdout,
        $(esc(x)),
        $(esc(y)),
        $(string(x)),
        $(string(y)),
    ))
end


# Recursively compare contents of similar types
function _rcompare(pass, x::T, y::T; use_cuda) where {T}
    if use_cuda
        return pass && (x ≈ y) # CUDA doesn't always satisfy ==
    else
        return pass && (x == y)
    end
end

function _rcompare(pass, x::T, y::T; use_cuda) where {T <: NamedTuple}
    for pn in propertynames(x)
        pass &=
            _rcompare(pass, getproperty(x, pn), getproperty(y, pn); use_cuda)
    end
    return pass
end

"""
    rcompare(x::T, y::T) where {T <: NamedTuple}

Recursively compare given types via `==`.
Returns `true` if `x == y` recursively.
"""
rcompare(x::T, y::T; use_cuda) where {T <: NamedTuple} =
    _rcompare(true, x, y; use_cuda)
rcompare(x, y; use_cuda) = false

function test_compare(x, y; use_cuda)
    rcompare(x, y; use_cuda) || @rprint_diff(x, y)
    @test rcompare(x, y; use_cuda)
end

function test_kernel!(use_cuda; fused!, unfused!, X, Y)
    for x in X
        x .= map(_ -> rand(), x)
    end
    for y in Y
        y .= map(_ -> rand(), y)
    end
    X_fused = deepcopy(X)
    X_unfused = X
    Y_fused = deepcopy(Y)
    Y_unfused = Y
    fused!(X_fused, Y_fused)
    unfused!(X_unfused, Y_unfused)
    @testset "Test correctness of $(nameof(typeof(fused!)))" begin
        test_compare(X_fused, X_unfused; use_cuda)
        test_compare(Y_fused, Y_unfused; use_cuda)
    end
end

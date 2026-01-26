### Constraints
# 
# Constraints are specified by the user as
#    lx_i ≤   x[i]  ≤ ux_i  # variable (box) constraints
#    lc_i ≤ c(x)[i] ≤ uc_i  # linear/nonlinear constraints
# and become equality constraints with l_i = u_i. ±∞ are allowed for l
# and u, in which case the relevant side(s) are unbounded.
# 
# The user supplies functions to calculate c(x) and its derivatives.
# 
# Of course we could unify the box-constraints into the
# linear/nonlinear constraints, but that would force the user to
# provide the variable-derivatives manually, which would be silly.
# 
# This parametrization of the constraints gets "parsed" into a form
# that speeds and simplifies the IPNewton algorithm, at the cost of many
# additional variables. See `parse_constraints` for details.

struct ConstraintBounds{T}
    nc::Int            # Number of linear/nonlinear constraints supplied by user
    # Box-constraints on variables (i.e., directly on x)
    eqx::Vector{Int}   # index-vector of equality-constrained x (not actually variable...)
    valx::Vector{T}    # value of equality-constrained x
    ineqx::Vector{Int} # index-vector of other inequality-constrained variables
    σx::Vector{Int8}   # ±1, in constraints σ(v-b) ≥ 0 (sign depends on whether v>b or v<b)
    bx::Vector{T}      # bound (upper or lower) on variable
    # Linear/nonlinear constraint functions and bounds
    eqc::Vector{Int}   # index-vector equality-constrained entries in c
    valc::Vector{T}    # value of the equality-constraint
    ineqc::Vector{Int} # index-vector of inequality-constraints
    σc::Vector{Int8}   # same as σx, bx except for the nonlinear constraints
    bc::Vector{T}
end
function ConstraintBounds(lx, ux, lc, uc)
    _cb(symmetrize(lx, ux)..., symmetrize(lc, uc)...)
end
function _cb(lx::AbstractArray{Tx}, ux::AbstractArray{Tx}, lc::AbstractVector{Tc}, uc::AbstractVector{Tc}) where {Tx,Tc}
    T = promote_type(Tx, Tc)
    ConstraintBounds{T}(length(lc), parse_constraints(T, lx, ux)..., parse_constraints(T, lc, uc)...)
end

Base.eltype(::Type{ConstraintBounds{T}}) where T = T
Base.eltype(cb::ConstraintBounds) = eltype(typeof(cb))

Base.convert(::Type{ConstraintBounds{T}}, cb::ConstraintBounds{S}) where {T,S} =
    ConstraintBounds(cb.nc, cb.eqx, convert(Vector{T}, cb.valx),
                     cb.ineqx, cb.σx, convert(Vector{T}, cb.bx),
                     cb.eqc, convert(Vector{T}, cb.valc), cb.ineqc,
                     cb.σc, convert(Vector{T}, cb.bc))


"""
    nconstraints(bounds) -> nc

The number of linear/nonlinear constraint functions supplied by the
user. This does not include bounds-constraints on variables.

See also: nconstraints_x.
"""
nconstraints(cb::ConstraintBounds) = cb.nc

"""
    nconstraints_x(bounds) -> nx

The number of "meaningful" constraints (not `±Inf`) on the x coordinates.

See also: nconstraints.
"""
function nconstraints_x(cb::ConstraintBounds)
    mi = isempty(cb.ineqx) ? 0 : maximum(cb.ineqx)
    me = isempty(cb.eqx) ? 0 : maximum(cb.eqx)
    nmax = max(mi, me)
    hasconstraint = falses(nmax)
    hasconstraint[cb.ineqx] .= true
    hasconstraint[cb.eqx] .= true
    sum(hasconstraint)
end

function Base.show(io::IO, cb::ConstraintBounds)
    indent = "    "
    print(io, "ConstraintBounds:")
    print(io, "\n  Variables:")
    showeq(io, indent, cb.eqx, cb.valx, 'x', :bracket)
    showineq(io, indent, cb.ineqx, cb.σx, cb.bx, 'x', :bracket)
    print(io, "\n  Linear/nonlinear constraints:")
    showeq(io, indent, cb.eqc, cb.valc, 'c', :subscript)
    showineq(io, indent, cb.ineqc, cb.σc, cb.bc, 'c', :subscript)
    nothing
end

# Synonym constructor for ConstraintBounds with no c(x)
BoxConstraints(lx, ux) = ConstraintBounds(lx, ux, [], [])

abstract type AbstractConstraints end

nconstraints(constraints::AbstractConstraints) = nconstraints(constraints.bounds)

struct OnceDifferentiableConstraints{F,J,T} <: AbstractConstraints
    c!::F         # c!(storage, x) stores the value of the constraint-functions at x
    jacobian!::J  # jacobian!(storage, x) stores the Jacobian of the constraint-functions
    bounds::ConstraintBounds{T}
end

function OnceDifferentiableConstraints(c!, jacobian!,
                                       lx::AbstractArray, ux::AbstractArray,
                                       lc::AbstractArray, uc::AbstractArray)
    b = ConstraintBounds(lx, ux, lc, uc)
    OnceDifferentiableConstraints(c!, jacobian!, b)
end

function OnceDifferentiableConstraints(lx::AbstractArray, ux::AbstractArray)
    bounds = ConstraintBounds(lx, ux, [], [])
    OnceDifferentiableConstraints(bounds)
end

function OnceDifferentiableConstraints(bounds::ConstraintBounds)
    c! = (c, x)->nothing
    J! = (J, x)->nothing
    OnceDifferentiableConstraints(c!, J!, bounds)
end

function checked_chunk(lx)
    if isempty(lx)
        throw(ArgumentError("autodiff on constraints require the full lower bound vector `lx`."))
    end
    ForwardDiff.Chunk(lx)
end

function OnceDifferentiableConstraints(c!, lx::AbstractVector, ux::AbstractVector,
                                       lc::AbstractVector, uc::AbstractVector,
                                       autodiff::Symbol = :central,
                                       chunk::ForwardDiff.Chunk = checked_chunk(lx))
                                       
    bounds = ConstraintBounds(lx, ux, lc, uc)
    T = eltype(bounds)
    sizex = size(lx)
    sizec = size(lc)

    xcache = zeros(T, sizex)
    ccache = zeros(T, sizec)

    backend = get_adtype(autodiff, chunk)
    jac_prep = DI.prepare_jacobian(c!, ccache, backend, xcache)
    function j!(_j, _x)
        DI.jacobian!(c!, ccache, _j, jac_prep, backend, _x)
        return _j
    end
    return OnceDifferentiableConstraints(c!, j!, bounds)
end


struct TwiceDifferentiableConstraints{F,J,H,T} <: AbstractConstraints
    c!::F # c!(storage, x) stores the value of the constraint-functions at x
    jacobian!::J # jacobian!(storage, x) stores the Jacobian of the constraint-functions
    h!::H   # h!(storage, x) stores the hessian of the constraint functions
    bounds::ConstraintBounds{T}
end

function TwiceDifferentiableConstraints(c!, jacobian!, h!, lx, ux, lc, uc)
    b = ConstraintBounds(lx, ux, lc, uc)
    TwiceDifferentiableConstraints(c!, jacobian!, h!, b)
end

function TwiceDifferentiableConstraints(c!, lx::AbstractVector, ux::AbstractVector,
    lc::AbstractVector, uc::AbstractVector,
    autodiff::Symbol = :central,
    chunk::ForwardDiff.Chunk = checked_chunk(lx))
    bounds = ConstraintBounds(lx, ux, lc, uc)
    T = eltype(bounds)
    nc = length(lc)
    nx = length(lx)
    x_example = zeros(T, nx)
    λ_example = zeros(T, nc)
    ccache = zeros(T, nc)
   
    function sum_constraints(_x, _λ)
        # TODO: get rid of this allocation with DI.Cache
        ccache_righttype = zeros(promote_type(T, eltype(_x)), nc)
        c!(ccache_righttype, _x)
        return sum(_λ[i] * ccache[i] for i in eachindex(_λ, ccache))
    end

    backend = get_adtype(autodiff, chunk)


    jac_prep = DI.prepare_jacobian(c!, ccache, backend, x_example)
    function con_jac!(_j, _x)
        DI.jacobian!(c!, ccache, _j, jac_prep, backend, _x) 
        return _j
    end

    hess_prep = DI.prepare_hessian(sum_constraints, backend, x_example, DI.Constant(λ_example))
    function con_hess!(_h, _x, _λ)
        DI.hessian!(sum_constraints, _h, hess_prep, backend, _x, DI.Constant(_λ)) 
        return _h
    end
    
    return TwiceDifferentiableConstraints(c!, con_jac!, con_hess!, bounds) 
end

function TwiceDifferentiableConstraints(c!, con_jac!,lx::AbstractVector, ux::AbstractVector,
    lc::AbstractVector, uc::AbstractVector,
    autodiff::Symbol = :central,
    chunk::ForwardDiff.Chunk = checked_chunk(lx))
    # TODO: is con_jac! still useful? we ignore it here
    
    return TwiceDifferentiableConstraints(c!, lx, ux, lc, uc, autodiff, chunk)
end



function TwiceDifferentiableConstraints(lx::AbstractArray, ux::AbstractArray)
    bounds = ConstraintBounds(lx, ux, [], [])
    TwiceDifferentiableConstraints(bounds)
end

function TwiceDifferentiableConstraints(bounds::ConstraintBounds)
    c! = (x, c)->nothing
    J! = (x, J)->nothing
    h! = (x, λ, h)->nothing
    TwiceDifferentiableConstraints(c!, J!, h!, bounds)
end


## Utilities

function symmetrize(l, u)
    if isempty(l) && !isempty(u)
        l = fill!(similar(u), -Inf)
    elseif !isempty(l) && isempty(u)
        u = fill!(similar(l), Inf)
    end
    # TODO: change to indices?
    size(l) == size(u) || throw(DimensionMismatch("bounds arrays must be consistent, got sizes $(size(l)) and $(size(u))"))
    _symmetrize(l, u)
end
_symmetrize(l::AbstractArray{T,N}, u::AbstractArray{T,N}) where {T,N} = l, u
_symmetrize(l::Vector{Any}, u::Vector{Any}) = _symm(l, u)
_symmetrize(l, u) = _symm(l, u)

# Designed to ensure that bounds supplied as [] don't cause
# unnecessary broadening of the eltype. Note this isn't type-stable; if
# the user cares, it can be avoided by supplying the same concrete
# type for both l and u.
function _symm(l, u)
    if isempty(l) && isempty(u)
        if eltype(l) == Any
            # prevent promotion from returning eltype Any
            l = Array{Union{}}(undef, 0)
        end
        if eltype(u) == Any
            u = Array{Union{}}(undef, 0)
        end
    end
    promote(l, u)
end

"""
    parse_constraints(T, l, u) -> eq, val, ineq, σ, b

From user-supplied constraints of the form

    l_i ≤  v_i  ≤ u_i

(which include both inequality and equality constraints, the latter
when `l_i == u_i`), convert into the following representation:

    - `eq`, a vector of the indices for which `l[eq] == u[eq]`
    - `val = l[eq] = u[eq]`
    - `ineq`, `σ`, and `b` such that the inequality constraints can be written as
             σ[k]*(v[ineq[k]] - b[k]) ≥ 0
       where `σ[k] = ±1`.

Note that since the same `v_i` might have both lower and upper bounds,
`ineq` might have the same index twice (once with `σ`=-1 and once with `σ`=1).

Supplying `±Inf` for elements of `l` and/or `u` implies that `v_i` is
unbounded in the corresponding direction. In such cases there is no
corresponding entry in `ineq`/`σ`/`b`.

T is the element-type of the non-Int outputs
"""
function parse_constraints(::Type{T}, l, u) where T
    size(l) == size(u) || throw(DimensionMismatch("l and u must be the same size, got $(size(l)) and $(size(u))"))
    eq, ineq = Int[], Int[]
    val, b = T[], T[]
    σ = Array{Int8}(undef, 0)
    for i = 1:length(l)
        li, ui = l[i], u[i]
        li <= ui || throw(ArgumentError("l must be smaller than u, got $li, $ui"))
        if li == ui
            push!(eq, i)
            push!(val, ui)
        else
            if isfinite(li)
                push!(ineq, i)
                push!(σ, 1)
                push!(b, li)
            end
            if isfinite(ui)
                push!(ineq, i)
                push!(σ, -1)
                push!(b, ui)
            end
        end
    end
    eq, val, ineq, σ, b
end

### Compact printing of constraints

function showeq(io, indent, eq, val, chr, style)
    if !isempty(eq)
        print(io, '\n', indent)
        if style == :bracket
            eqstrs = map((i, v)->"$chr[$i]=$v", eq, val)
        else
            eqstrs = map((i, v)->"$(chr)_$i=$v", eq, val)
        end
        foreach(s->print(io, s * ", "), eqstrs[1:end - 1])
        print(io, eqstrs[end])
    end
end

function showineq(io, indent, ineqs, σs, bs, chr, style)
    if !isempty(ineqs)
        print(io, '\n', indent)
        if style == :bracket
            ineqstrs = map((i, σ, b)->"$chr[$i]" * ineqstr(σ, b), ineqs, σs, bs)
        else
            ineqstrs = map((i, σ, b)->"$(chr)_$i" * ineqstr(σ, b), ineqs, σs, bs)
        end
        foreach(s->print(io, s * ", "), ineqstrs[1:end - 1])
        print(io, ineqstrs[end])
    end
end
ineqstr(σ, b) = σ > 0 ? "≥$b" : "≤$b"

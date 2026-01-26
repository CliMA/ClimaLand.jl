module MuladdMacro

"""
    postwalk(f, expr)

Apply the function `f` to each subexpression within expr.

Example:
```
julia> double(x::Int)=2x;  # Int Literals

julia> double(x::Any)=x;   # Everything else (e.g. Expressions and Symbols)

julia> MuladdMacro.postwalk(double, :(1+2*3));

julia> MuladdMacro.postwalk(double, :(1+2*3))
:(2 + 4 * 6)
```
"""
postwalk(f, x) = f(x)
postwalk(f, expr::Expr) = f(Expr(expr.head, (postwalk(f, arg) for arg in expr.args)...))

"""
    @muladd

Convert every combination of addition/subtraction and multiplication to a call of `muladd`.

If both of the involved operators are dotted, `muladd` is applied as a dot call.
The order of summation might be changed.
"""
macro muladd(ex)
    esc(to_muladd(ex))
end

"""
    to_muladd(ex)

Convert every combination of addition/subtraction and multiplication in expression `ex` to a call of `muladd`.

If both of the involved operators are dotted, `muladd` is applied as a dot call
The order of summation might be changed.
"""
function to_muladd(ex)
    postwalk(ex) do x
        # Modify summations
        issum(x) && return sum_to_muladd(x)

        # Modify subtractions
        issub(x) && return sub_to_muladd(x)

        return x
    end
end

"""
    sum_to_muladd(ex)

Replace sum `ex` by sequence of `muladd` if possible. Hereby the order of summation might be changed.
"""
function sum_to_muladd(ex)
    # Retrieve summands
    summands = args(ex)

    # Skip further calculations if no summand is a multiplication
    dotcall = isdotcall(ex)
    any(x -> ismul(x, dotcall), summands) || return ex

    # Split summands into two groups, one with expressions
    # of multiplications and one with other expressions
    mulsummands = filter(x -> ismul(x, dotcall), summands)
    oddsummands = filter(x -> !ismul(x, dotcall), summands)

    # If all summands are multiplications the first one is not reduced
    isempty(oddsummands) && push!(oddsummands, popfirst!(mulsummands))

    # Reduce sum to a composition of muladd
    foldl(mulsummands; init = newargs(ex, oddsummands...)) do s₁, s₂
        newmuladd(splitargs(s₂)..., s₁, dotcall)
    end
end

"""
    sub_to_muladd(ex)

Replace subtraction `ex` by `muladd` if possible.
"""
function sub_to_muladd(ex)
    # Retrieve operands
    x, y = args(ex)

    # Modify subtraction if possible
    dotcall = isdotcall(ex)
    if ismul(y, dotcall)
        y₁, y₂ = splitargs(y)
        return newmuladd(:(-$y₁), y₂, x, dotcall)
    elseif ismul(x, dotcall)
        return newmuladd(splitargs(x)..., :(-$y), dotcall)
    end

    return ex
end

"""
    iscall(ex, op)

Determine whether `ex` is a call of operation `op` with at least two arguments.
"""
iscall(ex::Expr, op) = ex.head == :call && length(ex.args) > 2 && ex.args[1] == op
iscall(ex, op) = false

"""
    isdotcall(ex)

Determine whether `ex` is a dot call.
"""
function isdotcall(ex::Expr)
    (ex.head == :. && length(ex.args) == 2 && Meta.isexpr(ex.args[2], :tuple)) ||
        (ex.head == :call && !isempty(ex.args) && startswith(string(ex.args[1]), '.'))
end
isdotcall(ex) = false

"""
    isdotcall(ex, op)

Determine whether `ex` is a dot call of operation `op` with at least two arguments.
"""
function isdotcall(ex::Expr, op)
    (ex.head == :. && length(ex.args) == 2 && ex.args[1] == op &&
     Meta.isexpr(ex.args[2], :tuple) && length(ex.args[2].args) > 1) ||
        (ex.head == :call && length(ex.args) > 2 && ex.args[1] == Symbol('.', op))
end
isdotcall(ex, op) = false

"""
    issum(ex)

Determine whether `ex` is a sum.
"""
issum(ex) = iscall(ex, :+) || isdotcall(ex, :+)

"""
    issub(ex)

Determine whether `ex` is a subtraction.
"""
issub(ex) = iscall(ex, :-) || isdotcall(ex, :-)

"""
    ismul(ex, dot::Bool)

Determine whether expression `ex` is a multiplication that is dotted if
`dot` is `true` and not dotted otherwise.
"""
ismul(ex, dot::Bool) = dot ? isdotcall(ex, :*) : iscall(ex, :*)

"""
    newmuladd(x, y, z, dot::Bool)

Return expression `(muladd).(x, y, z)` if `dot` is `true` and
`muladd(x, y, z)` otherwise.
"""
function newmuladd(x, y, z, dot::Bool)
    # Quoting seems to be required for the @. macro to work
    if dot
        :(($(Meta.quot(Base.muladd))).($x, $y, $z))
    else
        :($(Meta.quot(Base.muladd))($x, $y, $z))
    end
end

"""
    args(ex)

Return arguments of function call in `ex`.
"""
function args(ex::Expr)
    if ex.head == :call && length(ex.args) > 1
        return ex.args[2:end]
    end

    if ex.head == :. && length(ex.args) == 2 && Meta.isexpr(ex.args[2], :tuple) &&
       !isempty(ex.args[2].args)
        return ex.args[2].args
    end

    error("expression is not a function call with arguments")
end

"""
    splitargs(ex)

Split arguments of function call `ex` before last argument and combine first
arguments to one expression if possible.
"""
function splitargs(ex)
    if ex.head == :call && length(ex.args) > 2
        x = ex.args[2:(end - 1)]
        y = ex.args[end]
    elseif ex.head == :. && length(ex.args) == 2 && Meta.isexpr(ex.args[2], :tuple) &&
           length(ex.args[2].args) > 1
        x = ex.args[2].args[1:(end - 1)]
        y = ex.args[2].args[end]
    else
        error("cannot split arguments")
    end

    return newargs(ex, x...), y
end

"""
    newargs(ex, args)

Create new expression of function call `ex` with arguments `args`.

Unary function calls are not considered, i.e. if only one function
argument is provided it is returned.
"""
function newargs(ex::Expr, args...)
    # Return single argument
    length(args) == 1 && return args[1]

    # Create function calls with new arguments
    if ex.head == :call && !isempty(ex.args)
        return Expr(:call, ex.args[1], args...)
    end
    if ex.head == :. && length(ex.args) == 2 && Meta.isexpr(ex.args[2], :tuple)
        return Expr(:., ex.args[1], Expr(:tuple, args...))
    end

    error("expression is not a function call")
end

export @muladd

end # module

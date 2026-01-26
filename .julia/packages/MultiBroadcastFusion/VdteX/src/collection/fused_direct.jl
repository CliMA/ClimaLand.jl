#####
##### Simple version
#####

# General case: do nothing (identity)
transform(x) = x
transform(s::Symbol) = s
# Expression: recursively transform for Expr
function transform(e::Expr)
    if e.head == :macrocall && e.args[1] == Symbol("@__dot__")
        se = code_lowered_single_expression(e)
        margs = materialize_args(se)
        subexpr = :(Pair($(margs[1]), $(margs[2])))
        subexpr
    else
        Expr(transform(e.head), transform.(e.args)...)
    end
end

"""
    fused_direct(expr::Expr)

Directly transforms the input expression
into a tuple of `Pair`s, containing (firsts)
the destination of broadcast expressions and
(seconds) the broadcasted objects.

For example:
```julia
import MultiBroadcastFusion as MBF
using Test
expr_in = quote
    @. y1 = x1 + x2 + x3 + x4
    @. y2 = x2 + x3 + x4 + x5
end

expr_out = :(tuple(
    Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)),
    Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)),
))
@test MBF.fused_direct(expr_in) == expr_out
```

This can be used to make a custom kernel fusion macro:
```
import MultiBroadcastFusion as MBF
import MultiBroadcastFusion: fused_direct
MBF.@make_type MyFusedBroadcast
MBF.@make_fused fused_direct MyFusedBroadcast my_fused

Base.copyto!(fmb::MyFusedBroadcast) = println("You're ready to fuse!")

x1 = rand(3,3)
y1 = rand(3,3)
y2 = rand(3,3)

# 4 reads, 2 writes
@my_fused begin
  @. y1 = x1
  @. y2 = x1
end
```

Also see [`fused_assemble`](@ref)
"""
function fused_direct(expr::Expr)
    check_restrictions(expr)
    e = transform(expr)
    @assert e.head == :block
    ex = Expr(:call, :tuple, e.args...)
    # Filter out LineNumberNode, as this will not be valid due to prepending `tup = ()`
    linefilter!(ex)
    ex
end

function check_restrictions(expr::Expr)
    for _expr in expr.args
        _expr isa LineNumberNode && continue
        s_error = if _expr isa QuoteNode
            "Dangling symbols are not allowed inside fused blocks"
        elseif _expr.head == :for
            "Loops are not allowed inside fused blocks"
        elseif _expr.head == :if
            "If-statements are not allowed inside fused blocks"
        elseif _expr.head == :call
            "Function calls are not allowed inside fused blocks"
        elseif _expr.head == :(=)
            "Non-broadcast assignments are not allowed inside fused blocks"
        elseif _expr.head == :let
            "Let-blocks are not allowed inside fused blocks"
        elseif _expr.head == :quote
            "Quotes are not allowed inside fused blocks"
        else
            ""
        end
        isempty(s_error) || error(s_error)
        if _expr.head == :macrocall && _expr.args[1] == Symbol("@__dot__")
        else
            @show dump(_expr)
            error("Uncaught edge case")
        end
    end
end

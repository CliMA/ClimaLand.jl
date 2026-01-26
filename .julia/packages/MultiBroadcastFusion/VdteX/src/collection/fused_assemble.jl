#####
##### Fused assemble
#####

# General case: do nothing (identity)
transform_assemble(x, sym) = x
transform_assemble(s::Symbol, sym) = s
# Expression: recursively transform_assemble for Expr
function transform_assemble(e::Expr, sym)
    if e.head == :macrocall && e.args[1] == Symbol("@__dot__")
        se = code_lowered_single_expression(e)
        margs = materialize_args(se)
        subexpr = :($sym = ($sym..., Pair($(margs[1]), $(margs[2]))))
        subexpr
    else
        Expr(
            transform_assemble(e.head, sym),
            transform_assemble.(e.args, sym)...,
        )
    end
end

"""
    fused_assemble(expr::Expr)

Transforms the input expressions
into a runtime assembly of a tuple
of `Pair`s, containing (firsts)
the destination of broadcast expressions
and (seconds) the broadcasted objects.

For example:
```julia
import MultiBroadcastFusion as MBF
using Test
expr_in = quote
    @. y1 = x1 + x2 + x3 + x4
    @. y2 = x2 + x3 + x4 + x5
end

expr_out = quote
    tup = ()
    tup = (tup..., Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)))
    tup = (tup..., Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)))
    tup
end

@test MBF.linefilter!(MBF.fused_assemble(expr_in, :tup)) ==
      MBF.linefilter!(expr_out)
@test MBF.fused_assemble(expr_in, :tup) == expr_out
```

This can be used to make a custom kernel fusion macro:
```
import MultiBroadcastFusion as MBF
import MultiBroadcastFusion: fused_assemble
MBF.@make_type MyFusedBroadcast
MBF.@make_fused fused_assemble MyFusedBroadcast my_fused

Base.copyto!(fmb::MyFusedBroadcast) = println("You're ready to fuse!")

x1 = rand(3,3)
y1 = rand(3,3)
y2 = rand(3,3)

# 4 reads, 2 writes
@my_fused begin
    for i in 1:3
        @. y1 = x1
        @. y2 = x1
    end
end
```

Also see [`fused_direct`](@ref)
"""
fused_assemble(expr::Expr) = fused_assemble(expr, gensym())
function fused_assemble(expr::Expr, sym::Symbol)
    check_restrictions_assemble(expr)
    e = transform_assemble(expr, sym)
    @assert e.head == :block
    ex = Expr(:block, :($sym = ()), e.args..., sym)
    # Filter out LineNumberNode, as this will not be valid due to prepending `tup = ()`
    linefilter!(ex)
    ex
end

function check_restrictions_assemble(expr::Expr)
    for arg in expr.args
        arg isa LineNumberNode && continue
        s_error = if arg isa QuoteNode
            "Dangling symbols are not allowed inside fused blocks"
        elseif arg.head == :call
            "Function calls are not allowed inside fused blocks"
        elseif arg.head == :(=)
            "Non-broadcast assignments are not allowed inside fused blocks"
        elseif arg.head == :let
            "Let-blocks are not allowed inside fused blocks"
        elseif arg.head == :quote
            "Quotes are not allowed inside fused blocks"
        else
            ""
        end
        isempty(s_error) || error(s_error)

        if arg.head == :macrocall && arg.args[1] == Symbol("@__dot__")
        elseif arg.head == :for
            check_restrictions(arg.args[2])
        elseif arg.head == :if
            check_restrictions(arg.args[2])
        elseif arg.head == :macrocall && arg.args[1] == Symbol("@inbounds")
        else
            @show dump(arg)
            error("Uncaught edge case")
        end
    end
    return nothing
end

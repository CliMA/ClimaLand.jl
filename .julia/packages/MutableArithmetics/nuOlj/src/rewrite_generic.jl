# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# We need these two methods because we're changing how * is re-written.
operate!(::typeof(*), x::AbstractArray{T}, y::T) where {T} = (x .*= y)
operate!(::typeof(*), x::AbstractArray, y) = (x .= operate(*, x, y))

"""
    _rewrite_generic(stack::Expr, expr::T)::Tuple{Any,Bool}

This method is the heart of the rewrite logic. It converts `expr` into a mutable
equivalent, places any intermediate calculations onto `stack`, and returns a
tuple containing the return value---which is either `expr` or a `gensym`ed
variable equivalent to `expr`---and a boolean flag that indicates whether the
return value can be mutated by future callers.
"""
function _rewrite_generic end

"""
    _rewrite_generic(::Expr, x)

A generic fallback. Given a type `x` we return it without mutation. In addition,
this type should not be mutated by future callers.
"""
_rewrite_generic(::Expr, x) = esc(x), false

"""
    _rewrite_generic(::Expr, x::Number)

If `x` is a `Number` at macro expansion time, it _must_ be a constant literal.
We return `x` without mutation, but we return `true` because other callers may
mutate the value without fear. Put aother way, they don't need to wrap the value
in `copy_if_mutable(x)` before using it as the first argument to `operate!!`.

This most commonly happens in situations like `x^2`.
"""
_rewrite_generic(::Expr, x::Number) = x, true

function _is_generator(expr)
    return Meta.isexpr(expr, :call, 2) && Meta.isexpr(expr.args[2], :generator)
end

function _is_flatten(expr)
    return Meta.isexpr(expr, :call, 2) && Meta.isexpr(expr.args[2], :flatten)
end

function _is_parameters(expr)
    return Meta.isexpr(expr, :call, 3) && Meta.isexpr(expr.args[2], :parameters)
end

function _is_kwarg(expr, kwarg::Symbol)
    return Meta.isexpr(expr, :kw) && expr.args[1] == kwarg
end

function _rewrite_elseif!(if_expr, expr::Any)
    if expr isa Expr && Meta.isexpr(expr, :elseif)
        new_ifelse_expr = Expr(:elseif, esc(expr.args[1]))
        push!(if_expr.args, new_ifelse_expr)
        @assert 2 <= length(expr.args) <= 3
        return mapreduce(&, 2:length(expr.args)) do i
            return _rewrite_elseif!(new_ifelse_expr, expr.args[i])
        end
    end
    stack = quote end
    root, is_mutable = _rewrite_generic(stack, expr)
    push!(stack.args, root)
    push!(if_expr.args, stack)
    return is_mutable
end

_rewrite_generic(stack::Expr, expr::LineNumberNode) = expr, false

"""
    _rewrite_generic(stack::Expr, expr::Expr)

This method is the heart of the rewrite logic. It converts `expr` into a mutable
equivalent.
"""
function _rewrite_generic(stack::Expr, expr::Expr)
    if Meta.isexpr(expr, :block)
        new_stack = quote end
        for arg in expr.args
            root, _ = _rewrite_generic(new_stack, arg)
            push!(new_stack.args, root)
        end
        root = gensym()
        push!(stack.args, :($root = $new_stack))
        return root, false
    elseif Meta.isexpr(expr, :if)
        # `if` blocks are special, because we can't lift the computation inside
        # them into the stack; the values might be defined only if the branch is
        # true.
        if_expr = Expr(:if, esc(expr.args[1]))
        @assert 2 <= length(expr.args) <= 3
        is_mutable = mapreduce(&, 2:length(expr.args)) do i
            return _rewrite_elseif!(if_expr, expr.args[i])
        end
        root = gensym()
        push!(stack.args, :($root = $if_expr))
        return root, is_mutable
    elseif !Meta.isexpr(expr, :call)
        # In situations like `x[i]`, we do not attempt to rewrite. Return `expr`
        # and don't let future callers mutate.
        return esc(expr), false
    elseif Meta.isexpr(expr, :call, 1)
        # A zero-argument function
        return esc(expr), false
    elseif Meta.isexpr(expr.args[2], :(...))
        # If the first argument is a splat.
        return esc(expr), false
    elseif _is_generator(expr) || _is_flatten(expr) || _is_parameters(expr)
        if !(expr.args[1] in (:sum, :Σ, :∑))
            # We don't know what this is. Return the expression and don't let
            # future callers mutate.
            return esc(expr), false
        end
        # This is a generator expression like `sum(i for i in args)`. Generators
        # come in two forms: `sum(i for i=I, j=J)` or `sum(i for i=I for j=J)`.
        # The latter is a `:flatten` expression and needs additional handling,
        # but we delay this complexity for _rewrite_generic_generator.
        if Meta.isexpr(expr.args[2], :parameters)
            # The summation has keyword arguments. We can deal with `init`, but
            # not any of the others.
            p = expr.args[2]
            is_init = length(p.args) == 1 && _is_kwarg(p.args[1], :init)
            if is_init && expr.args[3] isa Expr
                # sum(iter ; init) form!
                # We rewrite only if `iter` is an Expr; if it's just a Symbol,
                # we don't enter this branch.
                root = gensym()
                init, _ = _rewrite_generic(stack, p.args[1].args[2])
                push!(stack.args, :($root = $init))
                return _rewrite_generic_generator(stack, :+, expr.args[3], root)
            end
            # We don't know how to deal with this
            return esc(expr), false
        else
            # Summations use :+ as the reduction operator.
            init_expr = expr.args[2].args[end]
            if Meta.isexpr(init_expr, :(=)) && init_expr.args[1] == :init
                # sum(iter, init) form!
                root = gensym()
                init, _ = _rewrite_generic(stack, init_expr.args[2])
                push!(stack.args, :($root = $init))
                new_expr = copy(expr.args[2])
                pop!(new_expr.args)
                return _rewrite_generic_generator(stack, :+, new_expr, root)
            elseif Meta.isexpr(expr.args[2], :flatten)
                # sum(iter for iter, init) form!
                first_generator = expr.args[2].args[1].args[1]
                init_expr = first_generator.args[end]
                if Meta.isexpr(init_expr, :(=)) && init_expr.args[1] == :init
                    root = gensym()
                    init, _ = _rewrite_generic(stack, init_expr.args[2])
                    push!(stack.args, :($root = $init))
                    new_expr = copy(expr.args[2])
                    pop!(new_expr.args[1].args[1].args)
                    return _rewrite_generic_generator(stack, :+, new_expr, root)
                end
            end
            return _rewrite_generic_generator(stack, :+, expr.args[2])
        end
    end
    # At this point, we have an expression like `op(args...)`. We can either
    # choose to convert the operation to it's mutable equivalent, or return the
    # non-mutating operation.
    if expr.args[1] == :+
        # +(args...) => add_mul(add_mul(arg1, arg2), arg3)
        @assert length(expr.args) > 1
        if length(expr.args) == 2  # +(arg)
            return _rewrite_generic(stack, expr.args[2])
        elseif length(expr.args) == 3 && _is_call(expr.args[3], :*)
            # +(x, *(y...)) => add_mul(x, y...)
            x, is_mutable = _rewrite_generic(stack, expr.args[2])
            rhs = if is_mutable
                Expr(:call, operate!!, add_mul, x)
            else
                Expr(:call, operate, add_mul, x)
            end
            for i in 2:length(expr.args[3].args)
                yi, _ = _rewrite_generic(stack, expr.args[3].args[i])
                push!(rhs.args, yi)
            end
            root = gensym()
            push!(stack.args, :($root = $rhs))
            return root, true
        end
        return _rewrite_generic_to_nested_op(stack, expr, add_mul)
    elseif expr.args[1] == :-
        # -(args...) => sub_mul(sub_mul(arg1, arg2), arg3)
        @assert length(expr.args) > 1
        if length(expr.args) == 2  # -(arg)
            value, is_mutable = _rewrite_generic(stack, expr.args[2])
            root = gensym()
            rhs = if is_mutable
                Expr(:call, operate!!, *, value, -1)
            else
                Expr(:call, operate, *, -1, value)
            end
            push!(stack.args, :($root = $rhs))
            return root, true
        end
        return _rewrite_generic_to_nested_op(stack, expr, sub_mul)
    elseif expr.args[1] == :*
        # *(args...) => *(*(arg1, arg2), arg3)
        @assert length(expr.args) > 2
        arg1, is_mutable = _rewrite_generic(stack, expr.args[2])
        arg2, _ = _rewrite_generic(stack, expr.args[3])
        rhs = if is_mutable
            Expr(:call, operate!!, *, arg1, arg2)
        else
            Expr(:call, operate, *, arg1, arg2)
        end
        root = gensym()
        push!(stack.args, :($root = $rhs))
        for i in 4:length(expr.args)
            arg, _ = _rewrite_generic(stack, expr.args[i])
            # It is always safe to modify `root` here
            rhs = Expr(:call, operate!!, *, root, arg)
            root = gensym()
            push!(stack.args, :($root = $rhs))
        end
        return root, is_mutable
    elseif expr.args[1] == :.+
        # .+(args...) => add_mul.(add_mul.(arg1, arg2), arg3)
        @assert length(expr.args) > 1
        if length(expr.args) == 2  # +(arg)
            return _rewrite_generic(stack, expr.args[2])
        end
        return _rewrite_generic_to_nested_op(
            stack,
            expr,
            add_mul;
            broadcast = true,
        )
    elseif expr.args[1] == :.-
        # .-(args...) => sub_mul.(sub_mul.(arg1, arg2), arg3)
        @assert length(expr.args) > 1
        if length(expr.args) == 2  # .-(arg)
            return _rewrite_generic(stack, Expr(:call, :.*, -1, expr.args[2]))
        end
        return _rewrite_generic_to_nested_op(
            stack,
            expr,
            sub_mul;
            broadcast = true,
        )
    else
        # Use the non-mutating call.
        result = Expr(:call, esc(expr.args[1]))
        for i in 2:length(expr.args)
            arg, _ = _rewrite_generic(stack, expr.args[i])
            push!(result.args, arg)
        end
        root = gensym()
        push!(stack.args, Expr(:(=), root, result))
        # This value isn't safe to mutate, because it might be a reference to
        # another object.
        return root, false
    end
end

function _rewrite_generic_to_nested_op(stack, expr, op; broadcast::Bool = false)
    root, is_mutable = _rewrite_generic(stack, expr.args[2])
    if !is_mutable
        # The first argument isn't mutable, so we need to make a copy.
        arg = Expr(:call, copy_if_mutable, root)
        root = gensym()
        push!(stack.args, Expr(:(=), root, arg))
    end
    for i in 3:length(expr.args)
        arg, _ = _rewrite_generic(stack, expr.args[i])
        rhs = if broadcast
            Expr(:call, broadcast!!, op, root, arg)
        else
            Expr(:call, operate!!, op, root, arg)
        end
        root = gensym()
        push!(stack.args, Expr(:(=), root, rhs))
    end
    return root, true
end

_is_call(expr, op) = Meta.isexpr(expr, :call) && expr.args[1] == op
_is_call(expr, op, n) = Meta.isexpr(expr, :call, 1 + n) && expr.args[1] == op

"""
    _rewrite_generic_generator(stack::Expr, op::Symbol, expr::Expr)

Special handling for generator expressions.

`op` is `:+` and `expr` is a `:generator` or `:flatten` expression.
"""
function _rewrite_generic_generator(
    stack::Expr,
    op::Symbol,
    expr::Expr,
    root = nothing,
)
    @assert op == :+
    is_flatten = Meta.isexpr(expr, :flatten)
    if is_flatten
        expr = expr.args[1]
    end
    # The value we're going to mutate. Start it off at `Zero`.
    if root === nothing
        root = gensym()
        push!(stack.args, Expr(:(=), root, Zero()))
    end
    # We need a new stack to go inside our for-loops since we want to
    # recursively rewrite the inner part as well.
    new_stack = quote end
    if _is_call(expr.args[1], op)
        # Optimization time! Instead of operate!!(op, root, op(args...)),
        # rewrite as operate!!(op, root, arg) for arg in args
        for arg in expr.args[1].args[2:end]
            value, _ = _rewrite_generic(new_stack, arg)
            rhs = Expr(:call, operate!!, add_mul, root, value)
            push!(new_stack.args, :($root = $rhs))
        end
    elseif op == :+ && _is_call(expr.args[1], :*)
        # Optimization time! Instead of operate!!(+, root, *(args...)), rewrite
        # this as operate!!(add_mul, root, args...)
        rhs = Expr(:call, operate!!, add_mul, root)
        for arg in expr.args[1].args[2:end]
            value, _ = _rewrite_generic(new_stack, arg)
            push!(rhs.args, value)
        end
        push!(new_stack.args, :($root = $rhs))
    elseif op == :+ && _is_call(expr.args[1], :-, 1)
        # Optimization time! Instead of operate!!(+, root, -(arg)), rewrite this
        # as operate!!(sub_mul, root, arg)
        value, _ = _rewrite_generic(new_stack, expr.args[1].args[2])
        rhs = Expr(:call, operate!!, sub_mul, root, value)
        push!(new_stack.args, :($root = $rhs))
    elseif is_flatten
        # The first argument is itself a generator
        _rewrite_generic_generator(new_stack, op, expr.args[1], root)
    else
        # expr.args[1] is the inner part of the loop. Rewrite it. We don't care
        # if it is mutable because we need a new value every iteration.
        inner, _ = _rewrite_generic(new_stack, expr.args[1])
        # Now build up the summation or product part of the inner loop. It's
        # always safe to mutate because we're going to start with `root=Zero()`.
        rhs = Expr(:call, operate!!, add_mul, root, inner)
        push!(new_stack.args, :($root = $rhs))
    end
    # This is a little complicated: walk back out of the generator statements
    # wrapping each level in a for loop and the over-writing the `new_stack`
    # variable.
    #
    # !!! warning
    #     The Julia syntax sum(i for i in 1:2, j in 1:i) is incorrect, but we
    #     handle it anyway! Because the user will write dependencies from left
    #     to right, we need to wrap from right to left.
    for i in length(expr.args):-1:2
        new_stack = _iterable_condition(new_stack, expr.args[i])
    end
    # Finally, push our new_stack onto the old `stack`...
    push!(stack.args, new_stack)
    # and return the `root`. We can mutate this in future because it started off
    # as `Zero`.
    return root, true
end

function _iterable_condition(new_stack, expr)
    if !Meta.isexpr(expr, :filter)
        return Expr(:for, esc(expr), new_stack)
    end
    body = quote
        if $(esc(expr.args[1]))
            $new_stack
        end
    end
    # A filter might be over multiple index sets
    for i in length(expr.args):-1:2
        body = Expr(:for, esc(expr.args[i]), body)
    end
    return body
end

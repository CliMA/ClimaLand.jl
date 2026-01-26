"""
    Atomix.@atomic

A superset of `Base.@atomic` supporting atomic operations on array elements.
Atomic operations on fields dispatches to `Base.@atomic`.
"""
:(@atomic)

macro atomic(ex)
    ans = atomic_impl(QuoteNode(:sequentially_consistent), ex)
    ans === nothing || return ans
    esc(:($Base.@atomic($ex)))
end

macro atomic(order, ex)
    ans = atomic_impl(order, ex)
    ans === nothing || return ans
    esc(:($Base.@atomic($order, $ex)))
end

macro atomic(a1, op, a2)
    ans = atomic_impl(QuoteNode(:sequentially_consistent), a1, op, a2)
    ans === nothing || return ans
    esc(:($Base.@atomic($a1, $op, $a2)))
end

macro atomic(order, a1, op, a2)
    ans = atomic_impl(order, a1, op, a2)
    ans === nothing || return ans
    esc(:($Base.@atomic($order, $a1, $op, $a2)))
end

function ref_expr(ex)
    @nospecialize
    isexpr(ex, :ref) || return nothing
    m = esc(ex.args[1])
    indices = map(esc, ex.args[2:end])
    return :(referenceable($m)[$(indices...)])
end

function order_expr(order)
    @nospecialize
    if order isa QuoteNode
        llvm_ordering_from_juila(order.value)
    else
        :(llvm_ordering_from_juila($(esc(order))))
    end
end

function atomic_impl(order, ex)
    @nospecialize
    if ex isa Expr
        if (ref = ref_expr(ex)) !== nothing
            return :(Atomix.get($ref, $(order_expr(order))))
        elseif isexpr(ex, :call, 3)
            return atomic_impl(order, ex.args[2], ex.args[1], ex.args[3])
        elseif ex.head === :(=)
            l, r = ex.args[1], esc(ex.args[2])
            if (ref = ref_expr(l)) !== nothing
                return :(Atomix.set!($ref, $r, $(order_expr(order))))
            end
        elseif length(ex.args) == 2
            shead = string(ex.head)
            if endswith(shead, '=')
                op = Symbol(shead[1:prevind(shead, end)])
                ans = atomic_impl(order, ex.args[1], op, ex.args[2])
                ans === nothing || return :($ans[2])
            end
        end
    end
    return nothing
end

function atomic_impl(order, a1, op, a2)
    @nospecialize
    ref = ref_expr(a1)
    ref === nothing && return nothing
    :(Atomix.modify!($ref, $(esc(op)), $(esc(a2)), $(order_expr(order))))
end

"""
    Atomix.@atomicswap

A superset of `Base.@atomicswap` supporting atomic operations on array elements.
Atomic operations on fields dispatches to `Base.@atomicswap`.
"""
:(@atomicswap)

macro atomicswap(ex)
    ans = atomicswap_impl(QuoteNode(:sequentially_consistent), ex)
    ans === nothing || return ans
    esc(:($Base.@atomicswap($ex)))
end

macro atomicswap(order, ex)
    ans = atomicswap_impl(order, ex)
    ans === nothing || return ans
    esc(:($Base.@atomicswap($order, $ex)))
end

function atomicswap_impl(order, ex)
    @nospecialize
    isexpr(ex, :(=), 2) || return nothing
    ref = ref_expr(ex.args[1])
    ref === nothing && return nothing
    val = esc(ex.args[2])
    return :(first(Atomix.swap!($ref, $val, $(order_expr(order)))))
end

"""
    Atomix.@atomicreplace

A superset of `Base.@atomicreplace` supporting atomic operations on array
elements.  Atomic operations on fields dispatches to `Base.@atomicreplace`.
"""
:(@atomicreplace)

macro atomicreplace(ex, old_new)
    order = QuoteNode(:sequentially_consistent)
    ans = atomicreplace_impl(order, order, ex, old_new)
    ans === nothing || return ans
    esc(:($Base.@atomicreplace($ex, $old_new)))
end

macro atomicreplace(order, ex, old_new)
    ans = atomicreplace_impl(order, order, ex, old_new)
    ans === nothing || return ans
    esc(:($Base.@atomicreplace($order, $ex, $old_new)))
end

macro atomicreplace(success_order, fail_order, ex, old_new)
    ans = atomicreplace_impl(success_order, fail_order, ex, old_new)
    ans === nothing || return ans
    esc(:($Base.@atomicreplace($success_order, $fail_order, $ex, $old_new)))
end

function atomicreplace_impl(success_order, fail_order, ex, old_new)
    @nospecialize
    ref = ref_expr(ex)
    ref === nothing && return nothing
    so = order_expr(success_order)
    fo = order_expr(fail_order)
    if isexpr(old_new, :call, 3) && old_new.args[1] === :(=>)
        exp, rep = esc(old_new.args[2]), esc(old_new.args[3])
        return :(Atomix.replace!($ref, $exp, $rep, $so, $fo))
    else
        old_new = esc(old_new)
        return :(Atomix.replace!($ref, $old_new::Pair..., $so, $fo))
    end
end

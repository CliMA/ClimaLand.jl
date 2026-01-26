module FastClosures

using Base.Meta

export @closure

struct Var
    name::Symbol
    num_esc::Int
end

function makelet(v::Var)
    ex = v.name
    for i=1:v.num_esc
        ex = esc(ex)
    end
    :($ex=$ex)
end

# Wrap `closure_expression` in a `let` block to improve efficiency.
function wrap_closure(module_, closure_expression)
    ex = macroexpand(module_, closure_expression)
    if isexpr(ex, :error)
        # There was an error in macroexpand - just return the original
        # expression so that the user gets a comprehensible error message.
        return closure_expression
    end
    if isexpr(ex, :do) && length(ex.args) >= 2 && isexpr(ex.args[2], :(->)) # do syntax
        ex = ex.args[2]
    end
    if isexpr(ex, :(->))
        args1 = ex.args[1]
        if isa(args1, Symbol)
            funcargs = [args1]
        else
            @assert isexpr(args1, :tuple)
            funcargs = args1.args
        end
    elseif isexpr(ex, :function)
        @assert isexpr(ex.args[1], :call)
        funcargs = ex.args[1].args[2:end]
    else
        throw(ArgumentError("Argument to @closure must be a closure!  (Got $closure_expression)"))
    end
    # FIXME support type assertions and kw args
    bound_vars = Var[Var(v,0) for v in funcargs]
    @assert isa(ex.args[2], Expr) && ex.args[2].head == :block
    captured_vars = Var[]
    find_var_uses!(captured_vars, bound_vars, ex.args[2], 0)
    quote
        let $(map(makelet, captured_vars)...)
            $closure_expression
        end
    end
end

"""
    @closure closure_expression

Wrap the closure definition `closure_expression` in a let block to encourage
the julia compiler to generate improved type information.  For example:

```julia
callfunc(f) = f()

function foo(n)
   for i=1:n
       if i >= n
           # Unlikely event - should be fast.  However, capture of `i` inside
           # the closure confuses the julia-0.6 compiler and causes it to box
           # the variable `i`, leading to a 100x performance hit if you remove
           # the `@closure`.
           callfunc(@closure ()->println("Hello \$i"))
       end
   end
end
```

There's nothing nice about this - it's a heuristic workaround for some
inefficiencies in the type information inferred by the julia 0.6 compiler.
However, it can result in large speedups in many cases, without the need to
restructure the code to avoid the closure.
"""
macro closure(ex)
    esc(wrap_closure(__module__, ex))
end

# Find arguments in closure arg list
#FIXME function find_closure_args(ex)
#end

# Utility function - fill `varlist` with all accesses to variables inside `ex`
# which are not bound before being accessed.  Variables which were bound
# before access are returned in `bound_vars` as a side effect.
#
# With works with the surface syntax so it unfortunately has to reproduce some
# of the lowering logic (and consequently likely has bugs!)
function find_var_uses!(varlist, bound_vars, ex, num_esc)
    if isa(ex, Symbol)
        var = Var(ex,num_esc)
        if !(var in bound_vars)
            var ∈ varlist || push!(varlist, var)
        end
        return varlist
    elseif isa(ex, Expr)
        if ex.head == :quote || ex.head == :line || ex.head == :inbounds
            return varlist
        end
        if ex.head == :(=)
            find_var_uses_lhs!(varlist, bound_vars, ex.args[1], num_esc)
            find_var_uses!(varlist, bound_vars, ex.args[2], num_esc)
        elseif ex.head == :kw
            find_var_uses!(varlist, bound_vars, ex.args[2], num_esc)
        elseif ex.head == :for || ex.head == :while || ex.head == :comprehension || ex.head == :let
            # New scopes
            inner_bindings = copy(bound_vars)
            find_var_uses!(varlist, inner_bindings, ex.args, num_esc)
        elseif ex.head == :try
            # New scope + ex.args[2] is a new binding
            find_var_uses!(varlist, copy(bound_vars), ex.args[1], num_esc)
            catch_bindings = copy(bound_vars)
            !isa(ex.args[2], Symbol) || push!(catch_bindings, Var(ex.args[2],num_esc))
            find_var_uses!(varlist,catch_bindings,ex.args[3], num_esc)
            if length(ex.args) > 3
                finally_bindings = copy(bound_vars)
                find_var_uses!(varlist,finally_bindings,ex.args[4], num_esc)
            end
        elseif ex.head == :call
            find_var_uses!(varlist, bound_vars, ex.args[2:end], num_esc)
        elseif ex.head == :local
            foreach(ex.args) do e
                if !isa(e, Symbol)
                    find_var_uses!(varlist, bound_vars, e, num_esc)
                end
            end
        elseif ex.head == :(::)
            find_var_uses_lhs!(varlist, bound_vars, ex, num_esc)
        elseif ex.head == :escape
            # In the 0.7-DEV churn, escapes persist during recursive macro
            # expansion until all macros are expanded.  Therefore, we
            # need to need to keep track of the number of escapes we're
            # currently inside, and to replay these when we construct the let
            # expression. See https://github.com/JuliaLang/julia/issues/23221
            # for additional pain and gnashing of teeth ;-)
            find_var_uses!(varlist, bound_vars, ex.args[1], num_esc+1)
        else
            find_var_uses!(varlist, bound_vars, ex.args, num_esc)
        end
    end
    varlist
end

find_var_uses!(varlist, bound_vars, exs::Vector, num_esc) =
    foreach(e->find_var_uses!(varlist, bound_vars, e, num_esc), exs)

# Find variable uses on the left hand side of an assignment.  Some of what may
# be variable uses turn into bindings in this context (cf. tuple unpacking).
function find_var_uses_lhs!(varlist, bound_vars, ex, num_esc)
    if isa(ex, Symbol)
        var = Var(ex,num_esc)
        var ∈ bound_vars || push!(bound_vars, var)
    elseif isa(ex, Expr)
        if ex.head == :tuple
            find_var_uses_lhs!(varlist, bound_vars, ex.args, num_esc)
        elseif ex.head == :(::)
            find_var_uses!(varlist, bound_vars, ex.args[2], num_esc)
            find_var_uses_lhs!(varlist, bound_vars, ex.args[1], num_esc)
        else
            find_var_uses!(varlist, bound_vars, ex.args, num_esc)
        end
    end
end

find_var_uses_lhs!(varlist, bound_vars, exs::Vector, num_esc) = foreach(e->find_var_uses_lhs!(varlist, bound_vars, e, num_esc), exs)


end # module

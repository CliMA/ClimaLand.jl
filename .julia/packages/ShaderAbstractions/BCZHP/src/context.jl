#=
We need to track the current OpenGL context.
Since we can't do this via pointer identity  (OpenGL may reuse the same pointers)
We go for this slightly ugly version.
=#
const ACTIVE_OPENGL_CONTEXT = Base.RefValue{Any}(nothing)

function current_context()
    isassigned(ACTIVE_OPENGL_CONTEXT) || error("No active context")
    ctx = ACTIVE_OPENGL_CONTEXT[]
    ctx === nothing && error("No active context")
    return ctx
end

function is_current_context(x)
    return x == ACTIVE_OPENGL_CONTEXT[]
end

function native_context_alive(x)
    error("Not implemented for $(typeof(x))")
end

"""
Is current context & is alive
"""
is_context_active(x) = is_current_context(x) && context_alive(x)

"""
Has context been destroyed or is it still living?
"""
context_alive(x) = native_context_alive(x)

function native_switch_context!(x)
    error("Not implemented for $(typeof(x))")
end

"""
Invalidates the current context
"""
function switch_context!()
    # for reverting to no context
    ACTIVE_OPENGL_CONTEXT[] = nothing
end

"""
Switches to a new context `x`. Is a noop if `x` is already current
"""
function switch_context!(x)
    if !is_current_context(x)
        ACTIVE_OPENGL_CONTEXT[] = x
        native_switch_context!(x)
    end
end

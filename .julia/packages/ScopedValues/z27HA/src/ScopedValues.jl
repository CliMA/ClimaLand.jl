module ScopedValues

export ScopedValue, with, @with, ScopedFunctor

if isdefined(Base, :ScopedValues)

import Base.ScopedValues: ScopedValue, with, @with, Scope, get
import Core: current_scope

macro enter_scope(scope, expr)
    Expr(:tryfinally, esc(expr), nothing, :(Scope($(esc(scope))::Union{Nothing, Scope})))
end

else

using HashArrayMappedTries
import HashArrayMappedTries: insert

"""
    ScopedValue(x)

Create a container that propagates values across dynamic scopes.
Use [`with`](@ref) to create and enter a new dynamic scope.

Values can only be set when entering a new dynamic scope,
and the value referred to will be constant during the
execution of a dynamic scope.

Dynamic scopes are propagated across tasks.

# Examples

```jldoctest
julia> const sval = ScopedValue(1);

julia> sval[]
1

julia> with(sval => 2) do
           sval[]
       end
2

julia> sval[]
1
```

!!! compat "Julia 1.11"
    Scoped values were introduced in Julia 1.11. In Julia 1.8+ a compatible
    implementation is available from the package ScopedValues.jl.
"""
mutable struct ScopedValue{T}
    const has_default::Bool
    const default::T
    ScopedValue{T}() where T = new(false)
    ScopedValue{T}(val) where T = new{T}(true, val)
    ScopedValue(val::T) where T = new{T}(true, val)
end

Base.eltype(::ScopedValue{T}) where {T} = T

"""
    isassigned(val::ScopedValue)

Test if the ScopedValue has a default value.
"""
Base.isassigned(val::ScopedValue) = val.has_default

const ScopeStorage = HAMT{ScopedValue, Any}

mutable struct Scope
    values::ScopeStorage
end

function Scope(parent::Union{Nothing, Scope}, key::ScopedValue{T}, value) where T
    val = convert(T, value)
    if parent === nothing
        values = ScopeStorage()
    else
        values = parent.values
    end
    values = insert(values, key, val)
    return Scope(values)
end

function Scope(scope, pairs::Pair{<:ScopedValue}...)
    for pair in pairs
        scope = Scope(scope, pair...)
    end
    return scope::Scope
end
Scope(::Nothing) = nothing

function Base.show(io::IO, scope::Scope)
    print(io, Scope, "(")
    first = true
    for (key, value) in scope.values
        if first
            first = false
        else
            print(io, ", ")
        end
        print(io, typeof(key), "@")
        show(io, Base.objectid(key))
        print(io, " => ")
        show(IOContext(io, :typeinfo => eltype(key)), value)
    end
    print(io, ")")
end

struct NoValue end
const novalue = NoValue()

"""
    get(val::ScopedValue{T})::Union{Nothing, Some{T}}

If the scoped value isn't set and doesn't have a default value,
return `nothing`. Otherwise returns `Some{T}` with the current
value.
"""
function get(val::ScopedValue{T}) where {T}
    # Inline current_scope to avoid doing the type assertion twice.
    scope = current_scope()
    if scope === nothing
        isassigned(val) && return Some(val.default)
        return nothing
    end
    scope = scope::Scope
    if isassigned(val)
        return Some(Base.get(scope.values, val, val.default)::T)
    else
        v = Base.get(scope.values, val, novalue)
        v === novalue || return Some(v::T)
    end
    return nothing
end

function Base.getindex(val::ScopedValue{T})::T where T
    maybe = get(val)
    maybe === nothing && throw(KeyError(val))
    return something(maybe)::T
end

function Base.show(io::IO, val::ScopedValue)
    print(io, ScopedValue)
    print(io, '{', eltype(val), '}')
    print(io, '(')
    v = get(val)
    if v === nothing
        print(io, "undefined")
    else
        show(IOContext(io, :typeinfo => eltype(val)), something(v))
    end
    print(io, ')')
end

"""
    with(f, var::ScopedValue{T} => val::T)

Execute `f` in a new scope with `var` set to `val`.
"""
function with(f, pair::Pair{<:ScopedValue}, rest::Pair{<:ScopedValue}...)
    @nospecialize
    scope = current_scope()
    scope = Scope(scope, pair...)
    for pair in rest
        scope = Scope(scope, pair...)
    end
    enter_scope(scope) do
        f()
    end
end

with(@nospecialize(f)) = f()

"""
    @with vars... expr

Macro version of `with(f, vars...)` but with `expr` instead of `f` function.
This is similar to using [`with`](@ref) with a `do` block, but avoids creating
a closure.
"""
macro with(exprs...)
    if length(exprs) > 1
        ex = last(exprs)
        exprs = exprs[1:end-1]
    elseif length(exprs) == 1
        ex = only(exprs)
        exprs = ()
    else
        error("@with expects at least one argument")
    end
    for expr in exprs
        if expr.head !== :call || first(expr.args) !== :(=>)
            error("@with expects arguments of the form `A => 2` got $expr")
        end
    end
    exprs = map(esc, exprs)

    return quote
        current_scope = $current_scope()
        scope = $(Scope)(current_scope, $(exprs...))
        @enter_scope scope begin
            $(esc(ex))
        end
    end
end

macro __tryfinally(ex, fin)
    Expr(:tryfinally,
        :($(esc(ex))),
        :($(esc(fin))),
    )
end

# Macro version of `Base.CoreLogging.with_logstate`
macro with_logstate(logstate, expr)
    quote
        t = $current_task()
        old = t.logstate
        @__tryfinally (t.logstate = logstate; $(esc(expr))) (t.logstate = old;)
    end
end

# Macro version of `enter_scope`
macro enter_scope(scope, expr)
    return quote
        logstate = $(Base.CoreLogging.LogState)($ScopePayloadLogger($current_logger(), $(esc(scope))))
        @with_logstate logstate begin
            $(esc(expr))
        end
    end
end

include("payloadlogger.jl")

end # isdefined

"""
    ScopedFunctor(f)

Create a functor that records the current dynamic scope, i.e. all current
`ScopedValue`s, along with `f`. When the functor is invoked, it runs `f`
in the recorded dynamic scope.
"""
struct ScopedFunctor{F}
    f::F
    scope::Union{Nothing, Scope}

    ScopedFunctor{F}(f) where {F} = new{F}(f, current_scope())
end
ScopedFunctor(f) = ScopedFunctor{typeof(f)}(f)
(sf::ScopedFunctor)() = @enter_scope sf.scope sf.f()

@deprecate scoped with

end # module ScopedValues

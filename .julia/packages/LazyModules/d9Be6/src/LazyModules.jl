module LazyModules

using Base: invokelatest
using Base: PkgId, UUID

export @lazy

const _LAZYMODE = Ref(true)

const _LOAD_LOCKER = Threads.ReentrantLock()

mutable struct LazyModule
    _lazy_pkgid::PkgId
    _lazy_loaded::Bool
end
LazyModule(uuid::UUID, name::String) = LazyModule(PkgId(uuid, name), false)
Base.Docs.Binding(m::LazyModule, v::Symbol) = Base.Docs.Binding(checked_import(m._lazy_pkgid), v)
function Base.show(io::IO, m::LazyModule)
    print(io, "LazyModule(", m._lazy_pkgid.name, ")")
end

struct LazyFunction
    pkgid::PkgId
    s::Symbol
end

function (f::LazyFunction)(args...; kwargs...)
    m = checked_import(f.pkgid)
    return invokelatest(getfield(m, f.s), args...; kwargs...)
end
function Base.show(io::IO, f::LazyFunction)
    print(io, "LazyFunction(", f.pkgid.name, ".", f.s, ")")
end
Base.Docs.aliasof(f::LazyFunction,   b) = Base.Docs.Binding(checked_import(f.pkgid), f.s)

function Base.getproperty(m::LazyModule, s::Symbol)
    if s in (:_lazy_pkgid, :_lazy_loaded)
        return getfield(m, s)
    end
    checked_import(m)
    lm = Base.root_module(getfield(m, :_lazy_pkgid))
    obj = getfield(lm, s)
    if obj isa Function
        return LazyFunction(getfield(m, :_lazy_pkgid), s)
    else
        return obj
    end
end

function checked_import(pkgid::PkgId)
    mod = if Base.root_module_exists(pkgid)
        Base.root_module(pkgid)
    else
        lock(_LOAD_LOCKER) do
            @debug "loading package: $(pkgid.name)"
            Base.require(pkgid)
        end
    end

    return mod
end

function checked_import(m::LazyModule)
    if !getfield(m, :_lazy_loaded)
        checked_import(getfield(m, :_lazy_pkgid))
        setfield!(m, :_lazy_loaded, true)
    end
    return m
end


"""
    @lazy import PkgName=UUID

Lazily import package `PkgName` with the actual loading delayed to the first usage.

```julia
module MyLazyPkg
    @lazy import Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80"
    draw_figure(data) = Plots.plot(data, title="MyPkg Plot")
end
```
"""
macro lazy(ex)
    if ex.head in (:import, :using)
        error("require `@lazy $(ex)=UUID` format.")
    end
    if ex.head != :(=)
        error("unrecognized expression: $(ex)")
    end
    uuid = UUID(ex.args[2])
    ex = ex.args[1]

    if !_LAZYMODE[]
        @info "disable lazy package loading: environment variable `LazyModules_lazyload=false` is detected" maxlog=1
        return ex
    end
    if ex.head != :import
        @warn "only `import` command is supported, fallback to eager mode"
        return ex
    end
    args = ex.args
    if length(args) != 1
        @warn "only single package import is supported, fallback to eager mode"
        return ex
    end
    x = args[1]
    if x.head == :.
        # usage: @lazy import Foo
        m = _lazy_load(__module__, x.args[1], uuid, x.args[1])
        # TODO(johnnychen94): the background eager loading seems to work only for Main scope
        isa(m, Module) && return m
        isnothing(m) && return ex
        return m
    elseif x.head == :(:)
        # usage: @lazy import Foo: foo, bar
        @warn "lazily importing symbols are not supported, fallback to eager mode"
        return ex
    elseif x.head == :as # compat: Julia at least v1.6
        # usage: @lazy import Foo as LazyFoo
        m = _lazy_load(__module__, x.args[2], uuid, x.args[1].args[1])
        isa(m, Module) && return m
        isnothing(m) && return ex
        return m
    else
        @warn "unrecognized syntax $ex"
        return ex
    end
end

function _lazy_load(mod, name::Symbol, uuid::UUID, sym::Symbol)
    if isdefined(mod, name)
        # otherwise, Revise will constantly trigger the constant redefinition warning
        m = getfield(mod, name)
        if m isa LazyModule || m isa Module
            return m
        else
            @warn "Failed to import module, the name `$name` already exists, do nothing"
            return nothing
        end
    end
    try
        m = LazyModule(uuid, String(sym))
        Core.eval(mod, :(const $(name) = $m))
        return m
    catch err
        @warn err
        return nothing
    end
end

"""
    require(m)

To avoid world-age issues, some package should be eagerly loaded first. This
function checks if the lazy package is loaded already.

The following example unlock the `fancy_color` feature unless users explicitly
load the Colors package.

```julia
using LazyModules

@lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

function fancy_color()
    LazyModules.require(Colors)
    return zero(Colors.RGB)
end
```
"""
function require(m::LazyModule)
    pkgid = m._lazy_pkgid
    if !Base.root_module_exists(pkgid)
        mname = pkgid.name
        error("$mname is required to be loaded first, maybe `using $mname` or `import $mname` and try again.")
    end
end

if VERSION < v"1.1"
    isnothing(x) = x === nothing
end

function __init__()
    if get(ENV, "LazyModules_lazyload", "") == "false"
        _LAZYMODE[] = false
    end
end
end # module

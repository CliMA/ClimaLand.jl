"""
    Domain(::Module)

Get the default domain for a module. If no domain has been defined for the
module, a new one is created.
"""
function Domain(__module__::Module)
    if !isdefined(__module__, :__nvtx_domain__)
        @eval __module__ begin
            const __nvtx_domain__ = $(Domain(string(__module__)))
        end
    end
    return Base.invokelatest(getproperty, __module__, :__nvtx_domain__)
end


file_lineno(__source__) = "$(__source__.file):$(__source__.line)"

# determine the domain and attributes for the macro call
function domain_attrs(__module__, default_message, args)
    domain = nothing
    message = nothing
    color = nothing
    category = nothing
    payload = nothing
    # if not a keyword, first arg is the message
    if length(args) >= 1
        arg = args[1]
        if !(arg isa Expr && arg.head == :(=))
            message = args[1]
            args = args[2:end]
        end
    end
    for arg in args
        if !(arg isa Expr && arg.head == :(=))
            error("$arg is not a keyword")
        end
        kw, val = arg.args
        if kw == :domain && isnothing(domain)
            domain = val
        elseif kw == :message && isnothing(message)
            message = val
        elseif kw == :color && isnothing(color)
            color = val
        elseif kw == :category && isnothing(category)
            category = val
        elseif kw == :payload && isnothing(payload)
            payload = val
        else
            if kw in [:domain, :message, :color, :category, :payload]
                error("$kw already defined")
            else
                error("invalid keyword $kw")
            end
        end
    end
    if isnothing(domain)
        domain = Domain(__module__)
    end
    if isnothing(message)
        message = default_message
    end
    if isnothing(color)
        # generate a unique color from the message
        color = hash(message) % UInt32
    end
    if domain isa Domain && message isa String
        # if domain and message are constant, using a StringHandle
        message = :(init!($(StringHandle(domain, message))))
    else
        message = esc(message)
    end
    # lazily initialize the domain
    domain = :(init!($(esc(domain))))
    color = esc(color)
    category = esc(category)
    payload = esc(payload)
    return domain, message, color, category, payload
end

"""
    NVTX.@mark [message] [domain=...] [color=...] [category=...] [payload=...]

Instruments an instantaneous event.

 - `message` is a string. Default is to use `"file:lineno"``. String
   interpolation is supported, but may incur additional overhead.
 - `domain` is a [`Domain`](@ref). Default is to use the default domain of the
   current module.
 - `color` is either a `Colorant` from the Colors.jl package, or an `UInt32`
   containing an ARGB32 value. Default is to generate one based on the hash of
   the message.
 - `category`: an integer describing the category of the event. Default is 0.
 - `payload`: an optional integer (`Int32`, `UInt32`, `Int64`, `UInt64`) or
   floating point (`Float32`, `Float64`) value to attach to the event.

# Example

```julia
for i = 1:10
    NVTX.@mark "iterate" payload=i
    do_work()
end
```
"""
macro mark(args...)
    domain, message, color, category, payload = domain_attrs(__module__, file_lineno(__source__), args)
    quote
        _isactive = isactive()
        if _isactive
            mark($domain; message=$message, color=$color, category=$category, payload=$payload)
        end
    end
end

"""
    NVTX.@range [message] [domain=...] [color=...] [category=...] [payload=...] expr

Instruments a range over the `expr`.

The default message is the expression, with file and line number. See
[`@mark`](@ref) for the other arguments.

# Example
```julia
for i = 1:10
    NVTX.@range "iterate" payload=i begin
        do_work()
    end
end
```
"""
macro range(args...)
    @assert length(args) >= 1
    expr = args[end]
    args = args[1:end-1]
    domain, message, color, category, payload = domain_attrs(__module__, "$expr $(file_lineno(__source__))", args)
    quote
        _isactive = isactive()
        if _isactive
            rangeid = range_start($domain; message=$message, color=$color, category=$category, payload=$payload)
        end
        # Use Expr(:tryfinally, ...) so we don't introduce a new soft scope (https://github.com/JuliaGPU/NVTX.jl/issues/28)
        # TODO: switch to solution once https://github.com/JuliaLang/julia/pull/39217 is resolved
        $(Expr(:tryfinally, esc(expr), :(_isactive && range_end(rangeid))))
    end
end

"""
    NVTX.@annotate [message] [domain=...] [color=...] [category=...] [payload=...] function ... end

Instruments a range a function definition, so that each invocation of the method
will be marked with a range. Equivalent to using [`@range`](@ref) within the
body of the function.

The default message is the function signature with file and line number.. See
[`@mark`](@ref) for the other arguments. Function arguments can be used as range
arguments.

# Example
```julia
NVTX.@annotate payload=x function foo(x)
    # function body
end
# is equivalent to
function foo(x)
    NVTX.@range payload=x begin
        # function body
    end
end

foo(1)
foo(2)
```
"""
macro annotate(args...)
    @assert length(args) >= 1
    expr = args[end]
    args = args[1:end-1]
    if expr.head in (:function, :(=)) && expr.args[1].head in (:call, :where)
        fsig = expr.args[1]
        body = expr.args[2]
    else
        error("NVTX.@annotate can only be applied to a function definition")
    end
    domain, message, color, category, payload = domain_attrs(__module__, "$fsig $(file_lineno(__source__))", args)
    quote
        $(esc(fsig)) = begin
            _isactive = isactive()
            if _isactive
                rangeid = range_start($domain; message=$message, color=$color, category=$category, payload=$payload)
            end
            try
                $(esc(body))
            finally
                if _isactive
                    range_end(rangeid)
                end
            end
        end
    end
end

"""
    NVTX.@category value name

Annotate an NVTX category `value` with `name` in the modules default domain.
For convenience, this returns the evaluated `value`

```julia
const category_3 = @category 3 "category 3"
```

A category should only be defined once for a particular value within a given domain.

See also [`name_category`](@ref)
"""
macro category(val, name)
    quote
        v = $(esc(val))
        if isactive()
            name_category($(Domain(__module__)), v, $(esc(name)))
        end
        v
    end
end

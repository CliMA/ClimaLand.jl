struct None end
const none = None()

struct FieldExpr
    isconst::Bool
    isatomic::Bool
    name::Symbol
    type::Any # none or Symbol or Expr
    default::Any # literal or Expr or none
    tags::Union{None, Expr}
end

function _expr(f::FieldExpr)
    nm_typ = f.type === none ? f.name : :($(f.name)::$(f.type))
    if f.isconst
        return Expr(:const, nm_typ)
    elseif f.isatomic
        return :(@atomic($nm_typ))
    else
        return nm_typ
    end
end

function _kw(f::FieldExpr)
    # don't include type to allow for more flexible inputs
    # and we can rely on implicit convert call in setfield!
    nm_typ = f.name # f.type === none ? f.name : :($(f.name)::$(f.type))
    return f.default === none ? nm_typ : Expr(:kw, nm_typ, f.default)
end

# parse a single-line field expression
# we progressively peel off layers of the expression to extract field information
# supported field expression types include:
# mutable struct Foo
#     no_type
#     with_type::Int
#     with_default = 1
#     with_type_default::Int = 1
#     with_tag &(xml=(key="with-tag",),)
#     with_tag_type::Int &(xml=(key="with-tag-type",),)
#     with_tag_default = 1 &(xml=(key="with-tag-default",),)
#     with_tag_type_default::Int = 1 &(xml=(key="with-tag-default",),)
#     @atomic no_type_atomic
#     @atomic with_type_atomic::Int
#     @atomic(with_default_atomic) = 1
#     @atomic(with_type_default_atomic::Int) = 1
#     @atomic(with_tag_atomic) &(xml=(key="with-tag-atomic",),)
#     @atomic(with_tag_type_atomic::Int) &(xml=(key="with-tag-type-atomic",),)
#     @atomic(with_tag_default_atomic) = 1 &(xml=(key="with-tag-default-atomic",),)
#     @atomic(with_tag_type_default_atomic::Int) = 1 &(xml=(key="with-tag-default-atomic",),)
#     const no_type_const
#     const with_type_const::Int
#     const with_default_const = 1
#     const with_type_default_const::Int = 1
#     const with_tag_const &(xml=(key="with-tag-const",),)
#     const with_tag_type_const::Int &(xml=(key="with-tag-type-const",),)
#     const with_tag_default_const = 1 &(xml=(key="with-tag-default-const",),)
#     const with_tag_type_default_const::Int = 1 &(xml=(key="with-tag-default-const",),)
# end
function FieldExpr(ex)
    isconst = isatomic = false
    name = Symbol()
    type = none
    default = none
    tags = none
    if Meta.isexpr(ex, :const)
        isconst = true
        ex = ex.args[1]
    end
    if Meta.isexpr(ex, :(=))
        def_and_tags = ex.args[2]
        if Meta.isexpr(def_and_tags, :call) && def_and_tags.args[1] == :&
            default = def_and_tags.args[2]
            tags = def_and_tags.args[3]
        else
            default = def_and_tags
        end
        ex = ex.args[1]
    end
    if Meta.isexpr(ex, :call) && ex.args[1] == :&
        tags = ex.args[3]
        ex = ex.args[2]
    end
    if Meta.isexpr(ex, :atomic)
        isatomic = true
        ex = ex.args[1]
    end
    if ex isa Symbol
        name = ex
    elseif Meta.isexpr(ex, :(::))
        name, type = ex.args
    else
        return nothing
    end
    name = Meta.isexpr(name, :escape) ? name.args[1] : name
    return FieldExpr(isconst, isatomic, name, type, default, tags)
end

function parsefields!(field_exprs::Vector{Any})
    fields = FieldExpr[]
    for (i, fex) in enumerate(field_exprs)
        fex isa LineNumberNode && continue
        f = FieldExpr(fex)
        if f !== nothing
            push!(fields, f)
            # replace field expression that may include defaults/tags w/ just the field name/type
            # + const/atomic annotations
            field_exprs[i] = _expr(f)
        else
            # ignore lines that aren't fields
        end
    end
    return fields
end

function parse_struct_def(kind, src, mod, expr)
    expr = macroexpand(mod, expr)
    Meta.isexpr(expr, :struct) || throw(ArgumentError("Invalid usage of @$kind macro"))
    # kind is: :noarg, :kwarg, :defaults, :tags
    ismutable, T, fieldsblock = expr.args
    if Meta.isexpr(T, :<:)
        T = T.args[1]
    end
    if Meta.isexpr(T, :curly)
        T_with_typeparams = copy(T)
        # sanitize T_with_typeparams to remove any type param bounds
        for i = 2:length(T_with_typeparams.args)
            if T_with_typeparams.args[i] isa Expr
                T_with_typeparams.args[i] = T_with_typeparams.args[i].args[1]
            end
        end
        typeparams = T.args[2:end]
        T = T.args[1]
    end
    ret = Expr(:block)
    # we always want to return original struct definition expression
    push!(ret.args, :($Base.@__doc__ $expr))
    # parse field exprs and sanitize field definitions
    fields = parsefields!(fieldsblock.args)
    if kind == :noarg
        ismutable || throw(ArgumentError("@noarg structs must be mutable"))
        if any(f.isconst for f in fields)
            #TODO: we could allow non-trailing const fields if they have default values
            # by setting the default value in the initial new() call
            # or if all fields have default values, we could allow const fields
            # because we'd have a value for each field to pass to new()
            throw(ArgumentError("const fields are not allowed in @noarg structs"))
        end
        # generate noarg constructor
        if @isdefined(T_with_typeparams)
            sig = Expr(:where, Expr(:call, T_with_typeparams), typeparams...)
            new_expr = Expr(:(=), :x, Expr(:call, Expr(:curly, :new, T_with_typeparams.args[2:end]...)))
        else
            sig = Expr(:call, T)
            new_expr = :(x = new())
        end
        cexpr = Expr(:function, sig, Expr(:block, src, new_expr))
        defs = [:(setfield!(x, $(Meta.quot(f.name)), $(f.default), $(f.isatomic ? Meta.quot(:sequentially_consistent) : Meta.quot(:not_atomic)))) for f in fields if f.default !== none]
        append!(cexpr.args[2].args, defs)
        push!(cexpr.args[2].args, Expr(:return, :x))
        # add inner constructor right after field definitions
        push!(expr.args[3].args, cexpr)
        #TODO: should we also generate an all-arg constructor like default struct constructors
        # that call convert to the field type for each field?
        # override StructUtils.noarg(::Type{nm}) = true and add outside struct definition
        push!(ret.args, :(StructUtils.noarg(::StructUtils.StructStyle, ::Type{<:$T}) = true))
        generate_field_defaults_and_tags!(ret, T, fields)
    elseif kind == :kwarg
        if !isempty(fields)
            # generate outer kwarg constructor, like: Foo(; a=1, b=2, ...) = Foo(a, b, ...)
            params = Expr(:parameters, (_kw(fex) for fex in fields)...)
            sig = Expr(:call, T, params)
            fexpr = Expr(:function, sig, Expr(:block, src, :(return $T($((f.name for f in fields)...)))))
            push!(ret.args, fexpr)
            if @isdefined(T_with_typeparams)
                # generate another kwarg constructor with type parameters
                sig = Expr(:where, Expr(:call, T_with_typeparams, params), typeparams...)
                fexpr = Expr(:function, sig, Expr(:block, src, :(return $T_with_typeparams($((f.name for f in fields)...)))))
                push!(ret.args, fexpr)
            end
        end
        # override StructUtils.kwarg(::Type{T}) = true and add outside struct definition
        push!(ret.args, :(StructUtils.kwarg(::StructUtils.StructStyle, ::Type{<:$T}) = true))
        generate_field_defaults_and_tags!(ret, T, fields)
    elseif kind == :nonstruct
        # Override StructUtils.structlike to return false and add outside struct definition
        push!(ret.args, :(StructUtils.structlike(::StructUtils.StructStyle, ::Type{<:$T}) = false))
        generate_field_defaults_and_tags!(ret, T, fields)
    else
        # if any default are specified, ensure all trailing fields have defaults
        # then generate multiple outer constructors allowing partial specification
        anydefaults = false
        for (i, f) in enumerate(fields)
            anydefaults |= f.default !== none
            if anydefaults && f.default === none
                throw(ArgumentError("All trailing fields must have default values in @$kind structs"))
            end
        end
        if anydefaults
            # find the first field with a default
            first_default_idx = findfirst(f -> f.default !== none, fields)
            if first_default_idx !== nothing
                # generate constructors for each possible number of arguments
                # from minimum required (non-default fields) up to just before all fields
                non_default_count = first_default_idx - 1
                total_fields = length(fields)
                
                # generate constructors that take k arguments where k goes from
                # non_default_count to (total_fields - 1)
                for num_args in non_default_count:(total_fields - 1)
                    if num_args >= 0
                        # constructor takes first num_args fields as arguments
                        sig = Expr(:call, T, (fields[i].name for i in 1:num_args)...)
                        # build the call to the full constructor
                        call_args = []
                        for (i, f) in enumerate(fields)
                            if i <= num_args
                                push!(call_args, f.name)
                            else
                                push!(call_args, f.default)
                            end
                        end
                        fexpr = Expr(:function, sig, Expr(:block, src, :(return $T($(call_args...)))))
                        push!(ret.args, fexpr)
                        
                        # also generate version with type parameters if applicable
                        if @isdefined(T_with_typeparams)
                            sig_tp = Expr(:where, Expr(:call, T_with_typeparams, (fields[i].name for i in 1:num_args)...), typeparams...)
                            fexpr_tp = Expr(:function, sig_tp, Expr(:block, src, :(return $T_with_typeparams($(call_args...)))))
                            push!(ret.args, fexpr_tp)
                        end
                    end
                end
            end
        end
        generate_field_defaults_and_tags!(ret, T, fields)
    end
    return esc(ret)
end

function generate_field_defaults_and_tags!(ret, T, fields)
    # generate fielddefaults override if applicable
    if any(f.default !== none for f in fields)
        defs_nt = Expr(:tuple, Expr(:parameters, [:(($(f.name)=$(f.default))) for f in fields if f.default !== none]...))
        push!(ret.args, :(StructUtils.fielddefaults(::StructUtils.StructStyle, ::Type{<:$T}) = $defs_nt))
    end
    # generate fieldtags override if applicable
    if any(f.tags !== none for f in fields)
        tags_nt = Expr(:tuple, Expr(:parameters, [:($(f.name)=$(f.tags)) for f in fields if f.tags !== none]...))
        push!(ret.args, :(StructUtils.fieldtags(::StructUtils.StructStyle, ::Type{<:$T}) = $tags_nt))
    end
end

const SHARED_MACRO_DOCS = """
The `@noarg`, `@kwarg`, `@defaults`, and `@tags` macros all support
specifying "field tags" for each field in a struct. Field tags are
a NamedTuple prefixed by `&` and are a way to attach metadata to a field. The
field tags are accessible via the `StructUtils.fieldtags` function, and certain
field tags are used by the `StructUtils.make` function to control how fields are
constructed, including:
  * `dateformat`: a `DateFormat` object to use when parsing or formatting a `Dates.TimeType` field
  * `lower`: a function to apply to a field when `applyeach` is called on a struct
  * `lift`: a function to apply to a field when `StructUtils.make` is called on a struct
  * `ignore`: a `Bool` to indicate if a field should be skipped/ignored when `applyeach` or `make` is called
  * `name`: a `Symbol` to be used instead of a defined field name in `applyeach` or used to match a field in `make`
  * `choosetype`: a function to apply to a field when `StructUtils.make` is called to determine the concrete type of an abstract or Union typed field

For example, the following struct definition includes a field with a `dateformat` tag:
```julia
@tags struct MyStruct
    date::Date &(dateformat=dateformat"yyyy-mm-dd",)
end
```
"""

"""
    @noarg mutable struct T
        ...
    end

Macro to enhance a `mutable struct` definition by automatically
generating an empty or "no-argument" constructor. Similar to the
`@kwarg` macro, default values can be specified for fields, which will
be set in the generated constructor. `StructUtils.noarg` trait is also
overridden to return `true` for the struct type. This allows
structs to easily participate in programmatic construction via
`StructUtils.make`.

Note that `const` fields are currently not allowed in `@noarg` structs.

$SHARED_MACRO_DOCS

Example
```julia
@noarg mutable struct Foo
    a::Int
    b::String
    c::Float64 = 1.0
    d::Vector{Int} = [1, 2, 3]
end
```

In the above example, the `@noarg` macro generates the following inner constructor:
```julia
function Foo()
    x = new()
    x.c = 1.0
    x.d = [1, 2, 3]
    return x
end
```
"""
macro noarg(expr)
    parse_struct_def(:noarg, __source__, __module__, expr)
end

"""
    @kwarg struct T
        ...
    end

Macro to enhance a `struct` definition by automatically generating a
keyword argument constructor. Default values can be specified for fields,
which will be set in the generated constructor. `StructUtils.kwarg` trait is
also overridden to return `true` for the struct type. This allows structs
to easily participate in programmatic construction via `StructUtils.make`.

$SHARED_MACRO_DOCS

Example
```julia
@kwarg struct Foo
    a::Int
    b::String = "foo"
    c::Float64 = 1.0
    d::Vector{Int} = [1, 2, 3]
end
```

In the above example, the `@kwarg` macro generates the following inner constructor:
```julia
function Foo(; a, b="foo", c=1.0, d=[1, 2, 3])
    return Foo(a, b, c, d)
end
```
"""
macro kwarg(expr)
    parse_struct_def(:kwarg, __source__, __module__, expr)
end

"""
    @defaults struct T
        ...
    end

Macro to enhance a `struct` definition by automatically generating an
outer constructor with default values for trailing fields. The generated
constructor will accept arguments for non-default fields and pass default
values to the inner constructor. `StructUtils.fielddefaults` trait is also
overridden to return a `NamedTuple` of default values for the struct type.

$SHARED_MACRO_DOCS

Example
```julia
@defaults struct Foo
    a::Int
    b::String = "foo"
    c::Float64 = 1.0
    d::Vector{Int} = [1, 2, 3]
end
```

In the above example, the `@defaults` macro generates the following outer constructor:
```julia
function Foo(a)
    return Foo(a, "foo", 1.0, [1, 2, 3])
end
```
"""
macro defaults(expr)
    parse_struct_def(:defaults, __source__, __module__, expr)
end

"""
    @tags struct T
        ...
    end

Macro to enhance a `struct` definition by allowing field tags to be
specified for each field.

$SHARED_MACRO_DOCS
"""
macro tags(expr)
    parse_struct_def(:tags, __source__, __module__, expr)
end

"""
    @nonstruct struct T
        ...
    end

Macro to mark a struct as not struct-like for StructUtils purposes. This macro
overrides `StructUtils.structlike` to return `false` for the struct type, which
means that `StructUtils.make` will not attempt to construct the struct using its
fields, but will instead use the `lift` function to convert the source directly
to the struct type.

This is useful for "unit" or "atom" types where the fields should be considered
private to the `make` process and the struct should be constructed by lifting
the source object directly.

**Note**: The `@nonstruct` macro does not support field defaults, field tags, or
other StructUtils macros (`@defaults`, `@tags`, `@noarg`, `@kwarg`) because by
using `@nonstruct`, you are explicitly opting out of StructUtils' struct-like
functionality. The struct's fields are considered private implementation details
for the `make` process.

Example
```julia
@nonstruct struct MyUnit
    value::String
end

# This will use StructUtils.lift to convert the source directly to MyUnit
# rather than trying to construct it from field values
x = StructUtils.make(MyUnit, "hello")
```
"""
macro nonstruct(expr)
    parse_struct_def(:nonstruct, __source__, __module__, expr)
end

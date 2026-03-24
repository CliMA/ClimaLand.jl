#=
**NOTE**: It is possible that the required forms of a bound function could evolve over time as ClimaLand
develops, though all utilities functionality should be preserved solely with adequate
changes to the compliance checks carried out within `construct_method_info!()`, should this occur.

The following functions are primarily used for automated handling, checking,
and documenting of different bounds and types that are compliant with the ConstrainedNeuralModel
architecture, and not their operation. This automates some bound compliance nuances away from the user, and enables
straightforward and reproducible model types and documentation in accordance with
our new conventions regarding neural parameterizations (see: https://github.com/CliMA/ClimaParams.jl/issues/277)

bound functions and functor types are specified by the user with the `@bound` or `@bound_type`
macros, respectively, which handle automated assessment of bound capabilities with the possible
configurations of ConstrainedNeuralModel types, and can inform the user on how to edit
bounds if specified bounds are non-compliant. The resulting info is stored in the module's
_BOUND_INFO_ or _BOUND_TYPES_ dictionaries, which are queried to build model metadata
and model documentation for model saving/storage via calls to
build_model_bound_documentation() or build_model_API().

Automated documentation assesses what modules are necessary (if any) beyond ClimaLand
to rebuild the same model in another codespace, creates an API of any affiliated
functions, and/or provides the actual written code for model bounds or any custom bound-types,
allowing for straightforward reproduction and reproducibility when sending saved models
to new users, or implementing models in different languages/systems.

Methods regarding the @bound_type macro and subsequent processing/validation:
- @bound_type
- construct_type_info!

Methods regarding the @bound macro and subsequent processing/validation:
- @bound
- construct_method_info!

Methods involving user-supplied bound types (functors):
- is_bound_type
- get_bound_type_info

Methods involving user-supplied bounds (functions/methods):
- is_valid_bound
- get_bound_info
- get_bound_evaluation_modes
- has_evaluation_mode
- check_evaluation_mode

Methods regarding automated documentation and metadata construction:
- build_parameter_metadata
- build_model_metadata
- build_API
- build_bound_docs
- assess_model_transferability
- get_fixed_layer_info
- populate_functor_methods!

Internal utilities for module mechanics and interaction within the julia ecosystem:
- bound_symbol
- get_val
- get_bounds
- get_modules
- get_methods
- instancemethods
- _check_children
- _is_under_main

Internal utilities for processing user code and associated documentation:
- _strip_location_comments
- _get_func_info
- _get_argtypes
- _get_declarative_docs
- _get_doc_from_method

=#
import InteractiveUtils: methodswith

"""
    _BOUND_INFO_::Dict{Symbol, Dict{Symbol, IdDict{Method, Dict}}}

Internal encapsulated dictionary for storing all valid bounds (functions or functor methods),
as well as information about them relevant for documentation, such as argtype, their actual
code syntax (string), their docstrings, whether they are configured for single or batched inputs,
etc. Is referenced to assess if a bound is valid for usage by ConstrainedNeuralModels, as well
as to build the metadata for model reproducibility when saving a ConstrainedNeuralModel.
"""
const _BOUND_INFO_ = Dict{Symbol, Dict{Symbol, IdDict{Method, Dict}}}(
    :FUNCTION => Dict{Symbol, IdDict{Method, Dict}}(),
    :FUNCTOR => Dict{Symbol, IdDict{Method, Dict}}(),
)


"""
    _BOUND_TYPES_::Dict{Symbol, Dict}()

Internal encapsulated dictionary for storing all valid functor types that can be used to
make bounds for ConstrainedNeuralModels. Stores relevant reproducibility information such
as supertype, mutability, their actual code syntax (string), their field types and names,
and type documentation. Is referenced to assess if a functor is valid for usage by 
ConstrainedNeuralModels, as well as to build the metadata for model reproducibility
when saving a ConstrainedNeuralModel.
"""
const _BOUND_TYPES_ = Dict{Symbol, Dict}()

"""
    is_bound_type(T::Union{Symbol, Type})

Checks if provided type is a valid bound functor type (i.e., if it is in _BOUND_TYPES_).
Defined for passing of a type or the symbolic binding of the type.
"""
is_bound_type(T::Union{Symbol, Type}) = haskey(_BOUND_TYPES_, bound_symbol(T))

"""
    get_bound_type_info(T::Union{Symbol, Type})

Retrieves information on the bound type, if it is a valid bound type.
"""
function get_bound_type_info(T::Union{Symbol, Type})
    if is_bound_type(T)
        return _BOUND_TYPES_[bound_symbol(T)]
    else
        error("""
        Type $(T) is not a valid bound type. Make sure its definition
        includes the @bound_type macro.
        """)
    end
end

"""
    bound_symbol(bound)

Simple utility to return what would be the key in _BOUND_INFO_ or _BOUND_TYPES_ for a given bound.
Handles various argtypes on account of varying usage.
"""
function bound_symbol(bound)
    if bound isa Function || bound isa Module
        return Symbol(bound)
    elseif bound isa Type
        return nameof(bound)
    elseif bound isa Method
        return bound.name
    elseif bound isa Symbol
        return bound
    else
        return nameof(typeof(bound))
    end
end

"""
    is_valid_bound(bound)

Checks if a passed bound is a valid bound for usage by a ConstrainedNeuralModel.
Handles various argtypes (a specific submethod, a name, the object itself, etc.)
Note: Will not work if trying to pass an instance of a callable type as a symbol.
"""
function is_valid_bound(bound)
    return haskey(_BOUND_INFO_[:FUNCTION], bound_symbol(bound)) ||
           haskey(_BOUND_INFO_[:FUNCTOR], bound_symbol(bound))
end

"""
    get_bound_info(bound)

Retrieves all methods and associated information on the given bound, if it is a valid bound.
Dispatched for arbitrary input (a function for a bound, its name, etc.)
"""
function get_bound_info(bound)
    b = bound_symbol(bound)
    if is_valid_bound(bound)
        #if valid, must be in one of them:
        if haskey(_BOUND_INFO_[:FUNCTOR], b)
            return _BOUND_INFO_[:FUNCTOR][b]
        else
            return _BOUND_INFO_[:FUNCTION][b]
        end
    else
        error(
            """
      No bound data exists for $(b) - make sure all methods have been defined with 
      the @bound macro. If $(b) is a functor, also make sure its declaration includes
      the @bound_type macro.
      """,
        )
    end
end


"""
    get_bound_info(bound::Symbol)

Retrieves the associated information for a specific bound method, if it is a valid bound.
Dispatched for specific method inputs (a sub-item of the _BOUND_INFO_ entry on a bound).
"""
function get_bound_info(method::Method)
    b_info = get_bound_info(method.name)
    if haskey(b_info, method)
        return b_info[method]
    else
        error(
            "Supplied sub-method $(method) is not a valid bound of $(method.name) - make sure the @bound macro has been used.",
        )
    end
end

"""
    check_evaluation_mode(input_types, bound)

Utility to warn the user if the evaluation modes of a given bound
only enable single-input or batched-input types exclusively.
"""
function check_evaluation_mode(input_types, bound)
    single_capable = :single in input_types
    batch_capable = :batched in input_types
    generic_capable = :generic in input_types
    if single_capable && !batch_capable && !generic_capable
        @warn """
        Supplied bound $(Symbol(bound)) is only configured for single-input usage.
        Make sure to only pass the model single inputs, or define an additional method
        if batched inputs are desired.
        """
    elseif !single_capable && batch_capable && !generic_capable
        @warn """
        Supplied bound $(Symbol(bound)) is only configured for batched-input usage.
        Make sure to only pass the model batched inputs, or define an additional method
        if single inputs are desired.
        """
    end
end

"""
   get_bound_evaluation_modes(bound_info::IdDict{Method, Dict})

Returns all possible evaluation modes of a given bound
"""
function get_bound_evaluation_modes(bound_info::IdDict{Method, Dict})
    modes = [x[:type] for x in values(bound_info)]
    return map(i -> getindex.(modes, i), eachindex(first(modes)))
end

"""
   get_bounds(c::Union{UpperOnly, LowerOnly})

Utility to return bounds from a given ConstraintType in
an iterable way. Dispatched for one-sided constraints.
"""
get_bounds(c::Union{UpperOnly, LowerOnly}) = (c.bound,)

"""
   get_bounds(c::TwoSided)

Utility to return bounds from a given ConstraintType in
an iterable way. Dispatched for two-sided constraints.
"""
get_bounds(c::TwoSided) = (c.upper_bound, c.lower_bound)

"""
   has_evaluation_mode(constraints, type::Symbol)

Checks if a all bounds in a given constraint set satisfy
a type of evaluation mode (only options are :static or :dynamic),
and assesses if that evaluation mode is sufficient for both
single- and batched-inputs.
"""
function has_evaluation_mode(constraints, type::Symbol)
    @assert type in [:static, :dynamic] """
    Can only check if bound types are :dynamic or :static.
    """
    results = []
    for bound in get_bounds(constraints)
        bound_info = get_bound_info(bound)
        inps, classes = get_bound_evaluation_modes(bound_info)
        check_evaluation_mode(inps, bound)
        push!(results, type in classes)
    end
    return all(results)
end

"""
   @bound_type (expr)

Macro to apply to the declaration of functor types used by ConstrainedNeuralModels.
Stores the code syntax of the declaration for reproducibility, as well as extracting
the type name from the declaration's abstract syntax tree (AST). This information
is then processed after declaring the type to create the _BOUND_TYPES_ entry for
the functor type; see `construct_type_info!()`
"""
macro bound_type(expr::Expr)
    #get the string of the actual type's code:
    code_lines = string(expr) #get the string
    code_str = "@bound_type " * _strip_location_comments(code_lines) #strip out julia-added tags to the string

    @assert expr.head == :struct "@bound_type can only be applied to type declaration."
    type_name = expr.args[2] #skip mutability arg
    #the name is always the base of the nested expression as a Symbol, so just iterate into expression until its not an expression anymore:
    while type_name isa Expr
        type_name = type_name.args[1]
    end
    type_data = :(Dict(:name => $(QuoteNode(type_name)), :code => $code_str))
    quote
        Base.@__doc__ ($(esc(expr))) #declares the type, associating the provided docstring with it, if any
        ConstrainedNeuralModels.construct_type_info!( #processes the type into _BOUND_TYPES_
            $(esc(type_name)),
            $(esc(type_data)),
        )
    end
end

"""
   construct_type_info!(type::Type, data::Dict)

Used by the `@bound_type` macro to process the information about
a bound type and build its entry in _BOUND_TYPES_ after the type
is declared and its docstring bound to the type. Collects and
stores information like the types and names of the fields of
the type, its documentation, mutability and supertype, as well
as the code syntax itself as parsed by the macro.
"""
function construct_type_info!(type::Type, data::Dict)
    type_params = collect(zip(fieldnames(type), fieldtypes(type)))
    type_doc = _get_declarative_docs(type)

    type_data = Dict(
        :name => data[:name],
        :code => data[:code],
        :params => type_params,
        :is_mutable => ismutable(type),
        :supertype => supertype(type),
        :typevars => _get_argtypes(type),
        :doc => type_doc,
    )
    #methods relating to this type populated later
    _BOUND_TYPES_[data[:name]] = type_data
    return nothing
end

"""
   populate_functor_methods!(type::Type)

Populates the :methods entry in _BOUND_TYPES_ for a valid
functor type with all methods related to that bound type,
including its constructors, any methods that take the type
as an argument, or the instance methods (the methods of instances
of the type, the functor methods) of the type. This cannot be
done at the time of its entry to _BOUND_TYPES_ as such methods
are usually defined after the declaration of the type.
"""
function populate_functor_methods!(type::Type)
    if is_bound_type(type)
        methodslist =
            vcat(methods(type), methodswith(type), instancemethods(type))
        _BOUND_TYPES_[nameof(type)][:methods] = Dict(
            (method = m, doc = _get_doc_from_method(m)) for m in methodslist
        )
    else
        error("$(type) is not a valid bound type.")
    end
    return nothing
end

"""
   @bound (expr)

Macro to apply to the declaration of function methods used by ConstrainedNeuralModels.
The macro serves to:
    - parse the function's abstract syntax tree (AST), extracting information
    like the function's name, its declared return type, whether it is a standalone
    function or a functor of a functor type, and the actual code syntax of
    the function (as a string) for reproducibility; see `_get_func_info()`
    - declare the function, and subsequently identify the correspondingly created
    Method (identified by Set differences, as new methods are not always appended
    to the end of the function's MethodsTable); see `get_methods()` and `instancemethods()`
This information is then processed after declaring the method to create the
_BOUND_INFO_ entry for the method, which tests to see if the method is compliant with
the needs of ConstrainedNeuralModels; see `construct_method_info!()`
"""
macro bound(expr::Expr)
    func_expr = expr
    #get the function's code as a string:
    code_lines = string(expr)
    code_str = "@bound " * _strip_location_comments(code_lines) #strip julia-added tags off the code string

    @assert func_expr.head in [:function, :(=)] "@bound can only be used on function definitions."
    #get the name of the thing being defined, whether its a functor or method, and the declared return type:
    func_name, bound_type, ret_type_expr = _get_func_info(func_expr)

    func_data = :(Dict(
        :code => $code_str,
        :name => $(QuoteNode(func_name)),
        :class => $(QuoteNode(bound_type)),
        :ret_type => $ret_type_expr,
    ))

    # Actually declare the function, and get the method associated with it and pass it to the check function:
    # Methods are not necessarily added to Julia's function tables in order (i.e. not appended to the end), so
    # the way to do this in a stable way (even if julia Base changes) is to compare what methods exist for the
    # name before and after the function declaration (the odd one out is the new one):
    quote
        local _before = try
            Set(get_methods($(esc(func_name)))) #get methods of the name before declaration (empty if no associated methods)
        catch
            Set{Method}()
        end
        @inline Base.@__doc__ ($(esc(func_expr))) #declare the function and make sure associated docstring is attached to it, if it exists
        local _after = Set(get_methods($(esc(func_name)))) #find the methods after declaration
        local _new_method = setdiff(_after, _before) #get the actual method
        @assert length(_new_method) == 1 "Expected exactly one new method"
        ConstrainedNeuralModels.construct_method_info!( #run the check function on the new method
            $(QuoteNode(func_name)),
            first(_new_method),
            $(esc(func_data)),
        )
    end
end

"""
    instancemethods(t::Type, mod = nothing)

Returns the list of methods of instances of a provided type (i.e. the
functor methods of a given functor type), either in the entire codespace
or within a specified module `mod`. This allows determination of the
methods of a callable type without needing to first instantiate an example,
which is otherwise necessary to see such methods with the `methods()`
function and not visible with the methodswith() function.

This is analagous to `methods()` for function names with multiple methods,
except with functor types. This has been a much-requested addition to Base julia
(https://github.com/JuliaLang/julia/pull/52304) and might exist in future
releases, at which point this function can be removed. For now, this implements
a simplified version specifically for the needs of the module.
"""
function instancemethods(t::Type, mod = nothing)
    searchtype = Tuple{t, Vararg{Any}} #search for all methods involving this type and any other args
    world = Base.get_world_counter()
    m_list = Base._methods_by_ftype(searchtype, -1, world) #same underlying syntax used by `methods()` but using a different search-type, so it should be a relatively stable function
    ms = Method[]
    #filter results by specified module:
    for m in m_list::Vector
        (isnothing(mod) || m.method.module ∈ mod) && push!(ms, m.method)
    end
    return ms
end

"""
    get_methods(x::Union{Function, Type})

A simple wrapper utility for usage by the @bound macro
in order to gather the right methods list depending on
whether the users applies the macro to a function declaration
or to the declaration of a functor method. If called on a
function, the methods(x) list is returned, and if called
on a callable type, the instancemethods(x) list is returned.
"""
function get_methods(x::Union{Function, Type})
    if x isa Function
        return methods(x)
    elseif x isa Type
        return instancemethods(x)
    else
        error("Invalid input for method search")
    end
end

"""
    _strip_location_comments(str::AbstractString)

A simple utility to filter out the the line and location
numbers in the string form of a function's code
via regex (i.e. just a regex filter on a string), leaving
behind the same code as typed by the user. Used by
`@bound` and `@bound_type`.
"""
function _strip_location_comments(str::AbstractString)
    replace(str, r"[ \t]*#=\s*[^=]*?:\d+\s*=#[ \t]*\n?" => "")
end

"""
    _get_func_info(expr::Expr)

Reads a declared function to determine if a function declaration is
a function method (:FUNCTION) or a method of a callable
type (:FUNCTOR), determine the declared return type of
the function, and the name of the function (or callable
type) the method is declared for.
"""
function _get_func_info(expr::Expr)
    sig = expr.args[1]
    #defines return-type defaults if user hasn't declared it
    where_params = [:(nothing)]
    ret_type = :(nothing)

    #All functions in julia with a parametric "where" will begin with the "where" Expr:
    if sig isa Expr && (sig.head == :where)
        #get the "where" declaraton so the return-type is fully known
        #(sig.args[1] is the actual function definition and return type, the rest pertains to the where)
        where_params = sig.args[2:end]
        #go to the next part of the function:
        sig = sig.args[1]
    end

    #If return type is delcared, all functions in julia will have that return type next as a head (func()::return_type, we look for the "::"):
    if sig isa Expr && (sig.head == :(::))
        #get the return type and combine it with the "where" if it exists (all but first arg pertain to this type)
        ret_type = Expr(:where, sig.args[2:end]..., where_params...)
        # go to the next part of the function (the actual name/args)
        sig = sig.args[1]
    end

    #All declared functions in julia have all their names/args expressed as a :call expression:
    @assert sig isa Expr && sig.head == :call
    if sig.args[1] isa Symbol #its a function method (the symbol is the function name)
        return sig.args[1], :FUNCTION, ret_type
    elseif sig.args[1] isa Expr && sig.args[1].head == :(::) #its a functor, the head is not a name but instead a certain type (i.e. name is like "(::MyType)")
        sig = sig.args[1].args[end] #skips the local name of the type's instance, if its defined, i.e. skips the "m" of `function (m::MyType)(pred, input)` since it is the first arg in a "::" expression
        if sig isa Symbol #its a basic type (no parametric types at the end)
            return sig, :FUNCTOR, ret_type
        else #its a parametric type (like "MyType{FT}"), in which case we just need the first part to get the type name
            return sig.args[1], :FUNCTOR, ret_type
        end
    else
        error(
            "Unexpected function type: @bound can only be applied to function or functor methods.",
        )
    end
end

"""
    _get_argtypes(x::UnionAll, transform = nothing)

Utility function to extract the individual argument types from a grouped
UnionAll type from a method signature, combining them
with the parametric types declared for the whole group. Once the
list is determined from unwrapping the UnionAll and thus identifying
the tuple of arg types (instead of the parametric types alone, which precede
the tuple of argument types), an optional functional transform can
be applied to each Type. Necessary in comparing method argument
types or the signatures for particular docstrings.
"""
function _get_argtypes(x::UnionAll, transform = nothing)
    u = x
    #move the external parametric types into each arg:
    while u isa UnionAll
        u = u.body
    end
    #get to the arg tuple (at the end)
    while u isa Union
        u = u.b #unions are stored as nested {{a, b}, b}... so you get to the end
    end
    #now get the fieldtypes from the resulting data type
    return _get_argtypes(u, transform) #not recursive, just dispatched for readability (see below)
end

"""
    _get_argtypes(x::Core.TypeofBottom, transform = nothing)

Dispatched form of _get_argtypes() for the case of an empty
signature (an instance Union{}, which is type Core.TypeofBottom).
"""
_get_argtypes(x::Core.TypeofBottom, transform = nothing) = ()

"""
    _get_argtypes(x::DataType, transform = nothing)

Dispatched form of _get_argtypes() for the case of a DataType
or a non-parametric signature instead of a grouped, parametric UnionAll.
Extracts the types of the fields in the function signature, with
the option to broadcast a final function over the tuple of Types.
"""
_get_argtypes(x::DataType, transform = nothing) =
    isnothing(transform) ? fieldtypes(x) : transform.(fieldtypes(x))

"""
    _get_declarative_docs(x; in_module::Union{Module, Nothing} = nothing)

Returns the documentation string pertaining specifically to the
declaration of a Type, as opposed to the docstrings associated with
callable methods of a Type or those of functions/methods, if it exists (and `nothing`
otherwise). The function will search within Main whenever `parentmodule(x)` 
is not defined, though the user can override these with the `in_module` arg. 
Does not apply for Types whose parent module is Core that are not extended 
in Base. For functions and Methods, see `_get_doc_from_method()` instead. 

    Examples:
    - _get_declarative_docs(ConstrainedNeuralModel)
    - _get_declarative_docs(Flux.Optimisers.RMSProp)
    - _get_declarative_docs(AbstractArray, in_module = Base)

This can also be used to grab the specific documentation of other declarations,
like the instance of a type (e.g., docstrings attached to the declaration/assignment
of a variable), a blank function interface (e.g., `:asdf` for `function asdf end`),
an anonymous function (but not a nullary function, e.g.`g = x -> x^3`, 
not `f() = (something)`), or anything defined with `const`, though for these
the Symbol of the variable's name must be passed, not the variable itself
or Symbol(variable) (i.e., if documentation was attached to the line `x = 3`
or `x = Flux.Bilinear(3, 4, 5)`, you must use `_get_declarative_docs(:x)`, 
not `_get_declarative_docs(Symbol(x))` nor `_get_declarative_docs(x)`). For such instances
or working with Symbol inputs, when working in a subsidiary module of Main
(e.g., you are running Main on a GPU), you will likely need to pass in_module = @__MODULE__
to achieve the desired result.

This is included since the @docs macro or Docs module always returns
ALL documentation associated with a particular binding, and does
not build in the capability to grab the specific documentation of a particular
signature. This is important in ConstrainedNeuralModels for attaching the
right documentation to the right methods, or grabbing the specific documentation
associated with the definition of a type, when building metadata for bounds and
bound types.

This works because if the declaration statement is succesful (i.e., it runs), the
same name (binding) can't also be bound to another declarative statement
(doing so would just be overwriting the previous declaration), so its associated
docstring will always be the only doc with a signature of Union{} (no signature,
while no-arg nullary functions have signature Tuple{}).
"""
function _get_declarative_docs(x; in_module::Union{Module, Nothing} = nothing)
    return try
        cond_1 = x isa Type || x isa Module || x isa Function #parentmodule() only exists for these types
        cond_2 = isnothing(in_module)
        mod = cond_2 ? (cond_1 ? parentmodule(x) : Main) : in_module
        b = Docs.Binding(mod, bound_symbol(x))
        Docs.meta(mod)[b].docs[Union{}].text[1]
    catch
        nothing
    end
end

"""
    _get_doc_from_method(m::Method)

Returns the documentation string pertaining specifically to the passed
method, if it exists (and `nothing` otherwise).

This method is included since the @docs macro or Docs module always returns
ALL documentation relevant to a specified binding, even when specifying
a particular signature type. `Docs.doc(binding, signature)` attempts to give
ALL docs that might be relevant to the signature type, which can return
multiple methods when writing dispatched methods for subtypes. Docs does
not build in the ability to get the doc of a particular method, as bindings
are attached to function names. This is important in ConstrainedNeuralModels
for attaching the right documentation to the right methods, when building
metadata for bounds and bound types.

Method signatures are compared against those in the documentation as strings,
due to the nuances of comparing TypeVars in julia (an example of this can be found
at https://discourse.julialang.org/t/equality-of-unwrapped-unionalls/128716, or
https://discourse.julialang.org/t/docstring-for-method/92675/6 regarding docs
signature matching). As the method signature upon declaration is the same form
of that which is added into the Docs multidoc dictionary, comparing types as
strings is consistent for finding the right docstring match.
"""
function _get_doc_from_method(m::Method)
    mod = parentmodule(m)
    binding = Docs.Binding(mod, m.name) #get how Docs stores this item
    if haskey(Docs.meta(mod), binding) #find it in that module's Docs
        docs_dict = Docs.meta(mod)[binding].docs #return all docs associated with that item
    else
        return nothing
    end
    #strings to avoid comparison nuances bewteen identical UnionAll types
    m_args = _get_argtypes(m.sig, string)[2:end] #self-arg is first for a method, which is not in the Docs signature
    argtypes = _get_argtypes.(keys(docs_dict), string) #no self-args here

    idx = findfirst(==(m_args), argtypes)
    if isnothing(idx)
        return nothing
    else
        docstrings = collect(values(docs_dict))
        return docstrings[idx].text[1]
    end
end

"""
    construct_method_info!(fname::Symbol, m::Method, data::Dict)

Used by the @bound macro to check and process the user's declared bound
method and check compliance for usage in ClimaLand simulations,
before creating a corresponding entry in _BOUND_INFO_. It is possible
that the required forms of a bound function could evolve over time as ClimaLand
develops, though all utilities functionality should be preserved with adequate
changes to the compliance checks carried out within this function.

Function args, as well as their input, output, and
return types are checked, as well as their eltypes. Functor methods are first
checked to see if their functor type is a valid bound type. If valid, the
method's associated docstring (if it exists), inferred evaluation mode, and
other information are bundled and added to the valid methods in _BOUND_INFO_.
"""
function construct_method_info!(fname::Symbol, m::Method, data::Dict)
    #Tests if bound method is compliant for usage in ClimaLand:
    #1) check if the type is a valid bound, if its a FUNCTOR:
    if data[:class] == :FUNCTOR && !is_bound_type(fname)
        error("""
        Type `$(fname)` is not a valid bound functor type.
        Please make to define it with the @bound_type macro.
        """)
    end
    m_argnames = split(m.slot_syms, "\0", keepempty = false) #instead of using the (undocumented) Base.method_argnames()
    #2) make sure the only args to the function are "pred" and "input":
    if m_argnames[2:end] != ["pred", "input"]
        error(
            """
      Insufficient bound function. 
      Please make sure your constraint method for `$(m.name)`
      takes only two args, in order: (pred = , input =).
      For extra args, custom extension of the module is required, though
      most capabilities involving extra args can be involved through the creation
      of a functor method and creating a bound type with @bound_type.
      """,
        )
    end
    pred_type, inp_type = fieldtypes(m.sig)[2:3]
    #3) Ensure user specified a return-type of their function:
    if isnothing(data[:ret_type])
        error(
            """
      Insufficient bound function.
      Please make sure you define a return type for $(m.name) taking
      argtypes $((pred_type, inp_type)) in accordance to the style of desired bound:
      - For generic bounds, return type should be <: typeof(pred)
      - For GPU-performant bounds with batched inputs, return type should match that of `pred`.
      - For GPU-performant bounds with single inputs, return type should be <:AbstractFloat.
      """,
        )
    end
    #4) Check if user-specifed argtypes that will ruin the compiler's ability
    # to specialize their float-type, making training immensely slow:
    generic_pred_check =
        occursin("<:AbstractFloat", string(pred_type)) &&
        !occursin("where", string(pred_type))
    generic_return_check =
        occursin("<:AbstractFloat", string(data[:ret_type])) &&
        !occursin("where", string(data[:ret_type]))
    if generic_pred_check || generic_return_check
        @warn """
        Detected non-parametric generic float-type arg in $(m.name) args/return. Using 
        generic types for the args and returns, such as 

            (pred::AbstractArray{<:AbstractArray}, input::AbstractArray{<:AbstractFloat})::AbstractArray{<:AbstractFloat}

        hints to the compiler that `pred` and `input` and the return type might all have
        different float types. This can cause immense slowdowns when training a model, as
        the compiler has no way to specialize calculated gradients. Unless this is 
        specifically the desired behvaior, ensure consistent float-type among all three
        via one parametric arg, e.g.:

            (pred::AbstractArray{T}, input::AbstractArray{T})::AbstractArray{T} where {T<:AbstractFloat}
            (pred::SVector{N, T}, input::SArray{S, T})::T where {N, S, T<:AbstractFloat}

        """
    end
    #start with generic bound-type until conditions indicate otherwise:
    bound_type = (:generic, :dynamic)
    if pred_type <: SArray
        #5) Check if compliant with single-input :static evaluation:
        if pred_type <: SVector
            inp_check = SVector{1, <:AbstractFloat} <: pred_type
            out_check = data[:ret_type] <: AbstractFloat
            @assert inp_check && out_check """
            Insufficient bound function. Ensure your arg `pred` of `$(m.name)`
            satisfies SVector{1, <:AbstractFloat} <: typeof(pred), and that
            the return type is a float.
            """
            bound_type = (:single, :static)
        else
            #6) Check if compliant with batched-input :static evaluation:
            inp_check =
                SMatrix{1, N, <:AbstractFloat, N} where {N <: Int} <: pred_type
            @assert inp_check """
            Insufficient bound function. Ensure your arg `pred` of `$(m.name)`
            satisfies (SMatrix{1, N, <:AbstractFloat, N} where {N<:Int}) <: typeof(pred),
            and that the return type is the same type as `pred`.
            """
            bound_type = (:batched, :static)
        end
    else
        #7) If not-static and single, check compliance:
        if pred_type <: Vector{<:AbstractFloat}
            @assert data[:ret_type] <: AbstractFloat "Ensure dynamic single-input function $(m.name) returns a float."
            bound_type = (:single, :dynamic)
            #8) If not-static and batched, check compliance:
        elseif pred_type <: Matrix{<:AbstractFloat}
            bound_type = (:batched, :dynamic)
        end
    end
    #9) Ensure return-type is compliant with the discovered bound-type (handled :single above already)
    if bound_type[1] != :single && !(data[:ret_type] <: pred_type)
        error(
            """
     Insufficient bound function. Make sure return type of $(m.name) taking argtypes
     $((pred_type, inp_type)) satifies {return_type} <: $(pred_type).
     """,
        )
    end
    #10) Ensure eltype of all args is the same
    t1 = eltype(pred_type) <: eltype(inp_type)
    t2 = eltype(data[:ret_type]) <: eltype(pred_type)
    @assert all([t1, t2]) """
    Insufficient bound function. Make sure return type of $(m.name) taking argtypes
    $((pred_type, inp_type)) satifies all args and return type having the same float type.
    """

    #Without an instance of the specified argtypes, we cannot run any more
    #assessments for :generic types, or if :dynamic outputs are row-like.

    #With a fully compliant bound, get the necessary info from the bound
    #and add it to _BOUND_INFO_.
    docstring = _get_doc_from_method(m)
    data[:method] = m
    data[:docs] = docstring
    data[:argtypes] = (input = inp_type, pred = pred_type)
    data[:type] = bound_type

    if haskey(_BOUND_INFO_[data[:class]], fname)
        _BOUND_INFO_[data[:class]][fname][m] = data
    else
        _BOUND_INFO_[data[:class]][fname] = IdDict(m => data)
    end
    return nothing
end

"""
    build_parameter_metadata(
        build_func::Flux.Optimisers.Restructure,
        creator_metadata::String,
        model_struct_filepath::String
    )::String

Creates the metadata description to include with saved parameter files for a
ConstrainedNeuralModel. Custom metadata strings can also be supplied and
added to the metadata saved with the parameters.

The description includes the name of the associated file for the structure
it was saved with, as well as a brief description of the model structure
and the index markers for breaking up the flattened parameter vector
into the different component parameter sets. Also includes a brief description
of how to load the ConstrainedNeuralModel from these parameters and a corresponding
build-function found in the associated structural save file.
"""
function build_parameter_metadata(
    build_func::Flux.Optimisers.Restructure,
    creator_metadata::String,
    model_struct_filepath::String,
)::String
    custom_metadata = ifelse(
        creator_metadata == "",
        "The creator has not specified any custom metadata for these parameters.",
        "The creator has specified the additional metadata for these parameters:\n" *
        creator_metadata,
    )
    ps_meta = """
    Flattened vector of trainable parameters of the predictive model for a ConstrainedNeuralModel type.
    The model structure looks like:
        $(string(build_func.model.predictive_model))
    Index markers for the different model pieces (using Flux layer notation) are as follows:
        $(string(build_func.offsets.predictive_model.layers))

    $(custom_metadata)

    MODEL BUILDING: model structure info was saved in \"$(basename(model_struct_filepath))\" when creating this file.

    To build the associated ConstrainedNeuralModel with these parameters, call `load_function()` with these
    parameters and the data of the model structure as follows:
        model = load_model(THIS_DATA_OBJECT["trainable_params"], "filepath/to/structure/data")
    """
    return ps_meta
end

"""
    get_modules(m::Module)

Internal utility function to grab the modules utilized by a given
module. Used to help assess if any of the internals of a saved
ConstrainedNeuralModel require additional modules to be loaded
into a new codespace.
"""
get_modules(m::Module) = ccall(:jl_module_usings, Any, (Any,), m)

"""
    assess_model_transferability(model::ConstrainedNeuralModel)

Recursively checks all subcomponents of a passed ConstrainedNeuralModel
to determine what modules (and methods, if `get_api_data = true`) are
needed for the model that are not capable within ClimaLand or this module.
Makes use of the recursive `_check_children()` utility.
"""
function assess_model_transferability(
    model::ConstrainedNeuralModel;
    get_api_data = false,
)
    base_list = vcat(
        Base,
        Core,
        get_modules(ConstrainedNeuralModels),
        get_modules(ClimaLand),
    )
    compare_list = Set(
        vcat(
            ConstrainedNeuralModels,
            ClimaLand,
            base_list,
            [get_modules(x) for x in base_list]...,
        ),
    )
    func_dict = Dict{Module, Dict{Type, Vector{Method}}}()
    for child in Flux.Functors.children(model)
        _check_children(child, func_dict, compare_list, get_api_data)
    end
    #without knowing a concrete example of the input type, we can't
    #check @code_typed to see if other modules' methods get called
    #inside the identified functions - while this might be unfortunate for
    #some cases, we wouldn't want to recurse any deeper than that anyway,
    #and the depth of this at present should adequately handle most cases.
    return func_dict
end

"""
    _check_children(child, data, comp_list, get_api_data)

Recursive function which assesses the modules associated with the
types of all present fields and compares them to a provided
comparison list. Any fieldtype not in the list makes a note for 
the type and/or the module, and, if `get_api_data = true`,
all methods of that type, methods using that type, and instancemethods (see `instancemethods()`)
of the type. All subcomponents of the component are then additionally
checked. These entries are aggregated downstream into an API should it be
added to the ConstrainedNeuralModel structural save file for
reproducibility (see `build_API()`, `build_model_API()`).
"""
function _check_children(child, data, comp_list, get_api_data)
    child_type_source = parentmodule(typeof(child))
    if !(child_type_source in comp_list)
        child_data = typeof(child).name.wrapper
        child_info =
            get_api_data ?
            collect(
                Set(
                    vcat(
                        collect(methods(child_data)),
                        collect(methodswith(child_data)),
                        collect(instancemethods(child_data)),
                    ),
                ),
            ) : []
        if haskey(data, child_type_source)
            if !haskey(data[child_type_source], child_data)
                data[child_type_source][child_data] = child_info
            end
        else
            data[child_type_source] = Dict(child_data => child_info)
        end
    end
    for sub_child in Flux.Functors.children(child)
        _check_children(sub_child, data, comp_list, get_api_data)
    end
    return nothing
end

"""
    get_fixed_layer_info(m::ConstrainedNeuralModel)

Simple utility function to build the string in the metadata
regarding the fixed layers of a ConstrainedNeuralModel. If
default layers were used, the string specifies as such, or
otherwise lists the layers and parameter values of the
custom fixed layers.
"""
function get_fixed_layer_info(m::ConstrainedNeuralModel)
    fixed_meta = "FIXED LAYERS: "
    if m.using_default_fixed_layers
        fixed_meta *= "The creator used the default fixed layers for this model.\n"
    else
        fixed_meta *= "The creator used custom fixed layers, defined as follows:\n"
        for l in m.fixed_layers.layers
            fixed_meta *= "LAYER: $(l)\n$(join(Flux.trainables(l), ", "))\n"
        end
    end
    return fixed_meta
end

"""
    build_API(data::Dict{Module, Dict{Type, Vector{Method}}})

Using the transferability data from `assess_model_transferability()`,
this function constructs a textual API listing the modules and types
needed for a given ConstrainedNeuralModel to work in a new codespace
that are not already provided by ClimaLand and/or this module. By module
and type, the declarative docs (see `_get_declarative_docs()`) are listed,
along with a description of all methods associated with that type and any
documentation associated with those methods (see `_get_doc_from_method()`).

Note: This will not capture/list any methods from other modules
nested within the methods found via `assess_model_transferability()`, such
as utility functions for readability. As such, we recommend thorough
documentation or including custom metadata with a ConstrainedNeuralModel and
its associated methods that lists/describes any such functions, if they are necessary.

Note: This only builds a reference API, and does not capture the code of
any of the listed methods. Only those marked with `@bound` or `@bound_type` will
have their code syntax made available, which is another reason why including
functional descriptions (or even the code syntax itself) in the documentation strings
for affiliated functions when building custom fixed layers, bound types, methods,
and predictive models is crucial.
"""
function build_API(data::Dict{Module, Dict{Type, Vector{Method}}})
    api = "API:\n Note that not every method listed will have documentation \
    or available code.\nIf this model was saved from the GPU, modules and types might be shown \
    as subchildren of their home modules.\n"
    for (mod, mod_dict) in data
        print_mod = _is_under_main(mod) ? Main : mod
        api *= "*IN MODULE - $(print_mod):\n"
        for (type, method_list) in mod_dict
            api *= "TYPE: $(nameof(type))\n"
            type_doc = _get_declarative_docs(type)
            api *= isnothing(type_doc) ? "" : type_doc
            api *= "\n"
            for m in sort(collect(method_list), by = x -> string(x.sig))
                d = _get_doc_from_method(m)
                if !isnothing(d)
                    api *= "\"\"\"\n" * d * "\n\"\"\"\n"
                end
                api *= split(string(m), "@")[1] * "\n\n"
            end
        end
        api *= "\n"
    end
    return api
end

"""
    get_val(::Val{B})

Extracts the value from a Val{} type.
"""
get_val(::Val{B}) where {B} = B

"""
    build_bound_docs(m::ConstrainedNeuralModel)

Builds the string describing the constraints of the passed
ConstrainedNeuralModel, grabbing the information on their
affiliated bounds from _BOUND_TYPES_ and _BOUND_INFO_.
The string includes a description of the bounds and/or bound
types, as well as their code syntax and other properties, to
create full reproducilbility of the code immediately surrounding
the saved model.
"""
function build_bound_docs(m::ConstrainedNeuralModel)
    doc = ""
    were_trainable = get_val(m.trainable_constraints)
    for bound in get_bounds(m.constraints)
        doc *= "BOUND: $(bound)\n"
        bound_info = get_bound_info(bound)
        if is_bound_type(nameof(typeof(bound)))
            bt = typeof(bound).name.wrapper
            populate_functor_methods!(bt)
            bt_info = get_bound_type_info(bt)
            params =
                join(["\t$(t[1])::$(t[2])" for t in bt_info[:params]], "\n")
            doc *= """
             Bound is a functor of type $(nameof(typeof(bound))):
                 SUPERTYPE: $(bt_info[:supertype])
                 MUTABLE: $(bt_info[:is_mutable])
                 TRAINABLE: $(were_trainable)
                 PARAMS:
                 $(params)
                 DOCSTRING:
                 $(bt_info[:doc])
                 CODE:
             $(bt_info[:code])

             All methods and documentation for $(nameof(typeof(bound))) can be found in the API.
             bound methods for this type are provided below:\n
             """
        else
            doc *= "Bound is a function, with the following bound methods:\n"
        end
        for (method, method_data) in
            sort(collect(bound_info), by = x -> string(x.first.sig))
            if !isnothing(method_data[:docs])
                doc *= "\"\"\"\n"
                doc *= method_data[:docs] * "\"\"\"\n"
            end
            doc *= method_data[:code] * "\n\n"
        end
    end
    return doc
end

"""
    _is_under_main(m::Module)

utility function to establish if a given module or submodule
is under Main. Used to distinguish between modules when running
on GPU, as the GPU can create nested submodules.
"""
function _is_under_main(m::Module)
    mod = m
    # traverse up the submodule tree until reaching
    # the root - if the root is Main, return true - if the root
    # of the module is a non-Main module, its parent will be
    # itself (e.g., parentmodule(Flux) = Flux), which is how
    # we know when to stop traversing and return false,
    # instead of iterating forever.
    while mod !== nothing
        mod === Main && return true
        mod === parentmodule(mod) && return false
        mod = parentmodule(mod)
    end
    return false
end

"""
    build_model_metadata(m::ConstrainedNeuralModel, creator_metadata::String, param_filepath::String)

Builds the metadata about a ConstrainedNeuralModel structure to accompany the saved model,
as well as its saved parameter file. The metadata describes required modules to load the model,
the code and descriptions of its associated bounds, any user-supplied metadata, and descriptions
of all model subcomponents. If an API of possible other methods is detected, an API is included
in the return as well.

Version numbers for Flux, ClimaLand, and julia when making the model are specified, as well as
a description on how to rebulid the model with known parameters (providing the name of its
associated parameter file)
"""
function build_model_metadata(
    m::ConstrainedNeuralModel,
    creator_metadata::String,
    param_filepath::String,
)
    transfer_data = assess_model_transferability(m)
    need_api = length(keys(transfer_data)) != 0
    needed_modules = filter(m -> !_is_under_main(m), keys(transfer_data))
    module_metadata = ifelse(
        length(needed_modules) > 0,
        "MODULE INFO: The following modules in addition to ClimaLand are needed to load this model: $(join(needed_modules, ", "))",
        "MODULE INFO: This model can be loaded directly with ClimaLand.ConstrainedNeuralModels without the addition of extra modules.",
    )
    api_metadata = ifelse(
        need_api,
        "CUSTOM TYPES: Additional methods or types are necessary for the loading of this ConstrainedNeuralModel. If creator did not specify any custom metadata, contact them about necessary steps for loading the model.\n",
        "CUSTOM TYPES: No additional methods or types are necessary for loading this model.\n",
    )
    custom_metadata = ifelse(
        creator_metadata == "",
        "CREATOR METADATA: The creator has not specified any custom metadata for this model type.",
        "CREATOR METADATA: The creator has specified the additional metadata for this model type:\n" *
        creator_metadata,
    )
    fixed_meta = get_fixed_layer_info(m)
    scaling_tag =
        m.scaling isa ConstScaling ?
        " ConstScaling - [$(join(m.scaling.in_scales, ","))]" : "NoScaling"
    m_meta = """
        $(module_metadata)
        $(api_metadata)
        FLOAT TYPE: $(eltype(m.out_scale.sc))
        CONSTRAINT TYPE: $(nameof(typeof(m.constraints)))

        $(fixed_meta)

        PREDICTIVE MODEL STRUCTURE:
        $(m.predictive_model)
        model parameters were saved in \"$(basename(param_filepath))\" when creating this file.

        INPUT SCALING: $(scaling_tag)
        OUTPUT SCALING: $(m.out_scale.sc[1])

        $(custom_metadata)

        This model was constructed using Julia v$(VERSION), Flux v$(pkgversion(Flux)), and ClimaLand v$(pkgversion(ClimaLand)).
        Trying to load with older versions may exhibit unexpected behavior.

        To build this ConstrainedNeuralModel with this setup, call `load_model()` on
        the path to this file or the data dictionary itself, and a Vector of floats 
        (same type as FLOAT TYPE) `params` of size $(sum(length, Flux.trainables(m.predictive_model))), as follows:
            
            model = load_model(params, (THIS DATA, OR FILEPATH))

        See `load_model()` documentation for more details.
        """
    return m_meta
end

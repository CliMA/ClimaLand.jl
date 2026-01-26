# Docs plugin (from DocStringExtensions, but it relies on Git, so...)
function codeblock(code)
    return "```julia\n$code\n```\n"
end

struct DefSignature end

"""
$(DefSignature())

Similar to `SIGNATURES` but has more precise method/struct
information obtained directly from the [`@pub`](@ref)
macro.
"""
const DEF = DefSignature()

function Docs.formatdoc(buf::IOBuffer, doc::Docs.DocStr, ::DefSignature)
    binding = doc.data[:binding]
    typesig = doc.data[:typesig]
    mod = doc.data[:module]
    object = Docs.resolve(binding)
    if !isdefined(mod, JIEKO_STUB)
        error("expected @pub in $mod")
    end
    stub = getfield(mod, JIEKO_STUB)

    name = binding.var
    if haskey(stub.consts, name)
        return print(buf, codeblock(stub.consts[name].doc))
    elseif haskey(stub.macros, name)
        return print(buf, codeblock(stub.macros[name].doc))
    elseif haskey(stub.structs, name)
        return print(buf, codeblock(stub.structs[name].doc))
    end

    typesig <: Tuple || error("expected typesig, got $typesig")
    type = Tuple{typeof(object),typesig.parameters...}
    sig, def = Any, nothing
    for (t, captured) in stub.interface
        if type <: t && t <: sig
            sig = t
            def = captured
        end
    end
    isnothing(def) && error("expected @pub on method $type in $mod")
    print(buf, codeblock(def.doc))
end

struct DefList end

"""
$DEF

List all the `@pub` definitions of a module. It shows
nothing if the binded object is not a module.
"""
const DEFLIST = DefList()

function Docs.formatdoc(buf::IOBuffer, doc::Docs.DocStr, ::DefList)
    binding = doc.data[:binding]
    mod = Docs.resolve(binding)
    mod isa Module || error("expected module, got $mod")
    # NOTE: not error here cuz we may not yet define
    # those @pub in the module
    isdefined(mod, JIEKO_STUB) || return nothing
    if isdefined(mod, :Prelude)
        println(
            buf,
            """
            ### Prelude

            Contains $mod.Prelude, all public definitions can be
            imported by `using $mod.Prelude`.
            """
        )
    end

    stub = getfield(mod, JIEKO_STUB)
    println(buf, "### Definitions\n\n")

    if !isempty(stub.consts)
        println(buf, "#### Constants\n\n")
        for (name, captured) in stub.consts
            println(buf, codeblock(captured.doc))
        end
    end

    if !isempty(stub.macros)
        println(buf, "#### Macros\n\n")
        for (name, captured) in stub.macros
            println(buf, codeblock(captured.doc))
        end
    end

    if !isempty(stub.structs)
        println(buf, "#### Structs\n\n")
        for (name, captured) in stub.structs
            println(buf, codeblock(captured.doc))
        end
    end

    if !isempty(stub.interface)
        println(buf, "#### Interfaces\n\n")
        for (type, captured) in stub.interface
            println(buf, codeblock(captured.doc))
        end
    end
end # format

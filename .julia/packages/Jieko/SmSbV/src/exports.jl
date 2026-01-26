"""
$DEF

Wrap all interfaces from the module in which this macro is called in a `Prelude` module
and export them. This is useful for exporting all interfaces from a module to be used
in other modules.

Optionally, you can pass a list of extra symbols to be exported along with the interfaces.

The concept of `Prelude` is borrowed from rust, where it is used to give users an explicit
way to import all the symbols from a package. This is useful for the package authors to
provide a default set of symbols to be used by the users via `using MyPackage.Prelude`
without giving them the easy way of importing everything from the package via `using MyPackage`.

Because `Prelude` module only contains the API symbols, it makes it easier for the toolchain
to check the APIs without mixing them with the implementation details.

# Example

```julia
module MyModule
using Jieko.Prelude # load everything you need from Jieko

@interface foo(x::Float64)::Int = 2

# export interface `foo` and some extra symbols
@prelude_module

end
```
"""
@pub macro prelude_module()
    return esc(prelude_module_m(__module__))
end

function relative_using(parent::Symbol, names::Symbol...)
    body = Expr(:(:), Expr(:., :., :., parent))
    for name in names
        push!(body.args, Expr(:., name))
    end
    return Expr(:using, body)
end

function prelude_module_m(mod::Module)
    isdefined(mod, JIEKO_STUB) || return xbaremodule(:Prelude, quote end)
    stub = getfield(mod, JIEKO_STUB)::JiekoStub
    stmts = expr_map(allcaptured(stub)) do captured
        quote
            $(relative_using(nameof(mod), captured.name))
            export $(captured.name)
        end
    end

    return Expr(:block,
        xbaremodule(:Prelude, quote
            $(relative_using(nameof(mod), nameof(mod)))
            export $(nameof(mod))
            $stmts
        end),
        emit_public(:Prelude),
    )

    return Expr(:block,
        Expr(:toplevel, Expr(:module, true, :Prelude, quote
            $(relative_using(nameof(mod), nameof(mod)))
            export $(nameof(mod))
            $stmts
            $extra_stmts
        end)),
        emit_public(:Prelude),
    )
end

function xbaremodule(name::Symbol, stmts::Expr)
    Expr(:toplevel, Expr(:module, true, :Prelude, stmts))
end

# interface, overloadable methods, docstring is the
# expression being collected by the macro

"""
$DEF

Mark a definition as public. This will help
the toolchain generating docstring and other things.

If definition is a function. It will be treated as an interface
definition as overloaded function cannot be denoted by public.
"""
macro pub(expr)
    return esc(pub_m(__module__, expr))
end

function pub_m(mod::Module, expr)
    expr = macroexpand(mod, expr)
    return quote
        $(emit_stub(mod))
        $Core.@__doc__ $expr
        $(emit_capture(mod, expr))
    end
end

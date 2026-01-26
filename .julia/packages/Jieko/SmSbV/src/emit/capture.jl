function emit_capture(mod::Module, expr::Expr)
    ret = capture_nothrow(mod, expr)
    isnothing(ret) && return xcall(Base, :error, "Unsupported expression: $expr")
    return ret
end

# NOTE: struct and macro does not need to support overloading
# thus just use the name as the key
function capture_nothrow(mod::Module, expr::Expr)
    if Meta.isexpr(expr, :struct)
        return capture_struct(mod, expr)
    elseif is_function(expr)
        return capture_function(mod, expr)
    elseif Meta.isexpr(expr, :macro)
        return capture_macro(mod, expr)
    elseif Meta.isexpr(expr, :const)
        return capture_const(mod, expr)
    elseif Meta.isexpr(expr, :block)
        for each in expr.args
            ret = capture_nothrow(mod, each)
            isnothing(ret) || return ret
        end
    else
        return nothing
    end
end

function capture_const(mod::Module, expr::Expr)
    body = expr.args[1].args[1]
    doc = sprint(print, body)
    name = name_only(body)

    return quote
        $JIEKO_STUB[$(QuoteNode(name))] = $CapturedConst($mod, $(QuoteNode(name)), $doc)
        $(emit_public(name))
    end
end

function capture_struct(mod::Module, expr::Expr)
    ismutable = expr.args[1]
    name = name_only(expr.args[2])
    doc = sprint(print, expr.args[2])
    stub = if ismutable
        CapturedStruct(mod, name, "mutable struct " * doc)
    else
        CapturedStruct(mod, name, "struct " * doc)
    end
    return quote
        $JIEKO_STUB[$(QuoteNode(name))] = $stub
        $(emit_public(name))
    end
end

function capture_function(mod::Module, expr::Expr)
    _, call, _ = split_function(expr)
    name = name_only(call)
    if Meta.isexpr(call, :(::))
        rettype = call.args[2]
        call = call.args[1]
        doc = sprint(print, call)
        doc *= " -> $rettype"
    else
        rettype = Any
        doc = sprint(print, call)
    end
    sig = split_signature(call)

    return quote
        $JIEKO_STUB[$sig] = $CapturedFunction($mod, $(QuoteNode(name)), $sig, $rettype, $doc)
        $(emit_method_check(mod, name))
        $(emit_public(name))
    end
end

function capture_macro(mod::Module, expr::Expr)
    head = expr.args[1]
    name = name_only(head)
    name = Symbol("@", name)
    doc = string("@", head.args[1])
    for idx in 2:length(head.args)
        doc *= " " * marg2hint(head.args[idx])
    end
    return quote
        $JIEKO_STUB[$(QuoteNode(name))] = $(CapturedMacro(mod, name, doc))
        $(emit_public(name))
    end
end

function marg2hint(expr)::String
    if expr isa Symbol
        return "<$expr>"
    elseif Meta.isexpr(expr, :(::))
        return "<$(expr.args[1])::$(expr.args[2])>"
    elseif Meta.isexpr(expr, :...)
        return marg2hint(expr.args[1]) * "..."
    elseif Meta.isexpr(expr, :kw)
        return "[$(marg2hint(expr.args[1])) = $(repr(expr.args[2]))]"
    else
        error("Unsupported macro argument: $(expr)")
    end
end

# General case: do nothing (identity)
substitute(x, code) = x
substitute(x::Core.SSAValue, code) = substitute(code[x.id], code)
substitute(x::Core.ReturnNode, code) = substitute(code[x.val.id], code)
substitute(s::Symbol, code) = s
# Expression: recursively substitute for Expr
substitute(e::Expr, code) =
    Expr(substitute(e.head, code), substitute.(e.args, Ref(code))...)

function code_info(expr::Expr)
    lc = Base.Meta.lower(Main, expr)
    if lc isa Symbol
        return code_info(lc)
    else
        return code_info(lc.args[1])
    end
end
code_info(s::Symbol) = s
code_info(ci::Core.CodeInfo) = ci.code
function code_lowered_single_expression(expr)
    code = code_info(expr) # vector
    if code isa Symbol
        return code
    else
        s = string(substitute(code[end], code))
        return Base.Meta.parse(s)
    end
end

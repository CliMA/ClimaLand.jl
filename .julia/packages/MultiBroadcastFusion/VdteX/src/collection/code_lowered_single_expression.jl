# General case: do nothing (identity)
substitute(x, code) = x
substitute(x::Core.SSAValue, code) = substitute(code[x.id], code)
substitute(x::Core.ReturnNode, code) = substitute(code[x.val.id], code)
substitute(s::Symbol, code) = s
# Expression: recursively substitute for Expr
substitute(e::Expr, code) =
    Expr(substitute(e.head, code), substitute.(e.args, Ref(code))...)

code_info(expr) = Base.Meta.lower(Main, expr).args[1]
function code_lowered_single_expression(expr)
    code = code_info(expr).code # vector
    s = string(substitute(code[end], code))
    return Base.Meta.parse(s)
end

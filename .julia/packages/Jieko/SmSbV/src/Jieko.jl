"""
$DEFLIST
"""
module Jieko

if VERSION < v"1.9-"
    using Base: @kwdef
end

using ExproniconLite: xcall, name_only, expr_map, split_function, is_function, split_signature

const JIEKO_STUB = Symbol("##Jieko#STUB#")

include("types.jl")
include("emit/stub.jl")
include("emit/check.jl")
include("emit/public.jl")
include("emit/capture.jl")
include("doc.jl")
include("pub.jl")
include("exports.jl")
include("reflect.jl")
include("err.jl")

var"##Jieko#STUB#"[Symbol("@pub")] = CapturedMacro(Jieko, Symbol("@pub"), "@pub <definition>")
var"##Jieko#STUB#"[:DEF] = CapturedConst(Jieko, :DEF, "DEF")
var"##Jieko#STUB#"[:DEFLIST] = CapturedConst(Jieko, :DEFLIST, "DEFLIST")

@static if VERSION > v"1.11-"
    @eval $(Expr(:public, Symbol("@pub")))
    @eval $(Expr(:public, :DEF))
    @eval $(Expr(:public, :DEFLIST))
end

@prelude_module

end # Jieko

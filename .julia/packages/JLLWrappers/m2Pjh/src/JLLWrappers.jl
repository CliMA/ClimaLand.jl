module JLLWrappers

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@compiler_options"))
    @eval Base.Experimental.@compiler_options compile=min optimize=0 infer=false
end

if VERSION >= v"1.6.0-DEV"
    using Preferences
end

const global_typeassert_available = VERSION >= v"1.9.0-"

# We need to glue expressions together a lot
function excat(exs::Union{Expr,Nothing}...)
    ex = Expr(:block)
    for exn in exs
        exn === nothing && continue
        if Meta.isexpr(exn, :block)
            append!(ex.args, exn.args)
        else
            push!(ex.args, exn)
        end
    end
    return esc(ex)
end

include("toplevel_generators.jl")
include("wrapper_generators.jl")
include("runtime.jl")

end # module

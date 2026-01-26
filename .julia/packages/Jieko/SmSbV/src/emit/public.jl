function emit_public(name::Symbol)
    # see JuliaLang/julia/issues/51450
    @static if VERSION > v"1.11-"
        return Expr(:public, name)
    else
        return nothing
    end
end

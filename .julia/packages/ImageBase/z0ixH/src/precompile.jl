macro warnpcfail(ex::Expr)
    modl = __module__
    file = __source__.file === nothing ? "?" : String(__source__.file)
    line = __source__.line
    quote
        $(esc(ex)) || @warn """precompile directive
     $($(Expr(:quote, ex)))
 failed. Please report an issue in $($modl) (after checking for duplicates) or remove this directive.""" _file=$file _line=$line
    end
end

function _precompile_()
    for ST in (N0f8, Float32, Float64)
        for T in (Gray{ST}, RGB{ST})
            @warnpcfail precompile(restrict, (Vector{T},))
            @warnpcfail precompile(restrict, (Matrix{T},))
        end
    end
end # _precompile_

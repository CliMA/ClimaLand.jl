function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    @assert precompile(Tuple{typeof(findfont),String})   # time: 0.12886831
    @assert precompile(Tuple{typeof(try_load),String})   # time: 0.033520337
    @assert precompile(Tuple{typeof(renderface),FTFont,Char,Int64})   # time: 0.019107351
end

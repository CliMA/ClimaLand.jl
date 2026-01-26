function _precompile_()
  ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
  for T in (Bool, Int, Float32, Float64)
    for A in (Vector, Matrix)
      precompile(stridedpointer, (A{T},))
      precompile(stridedpointer, (LinearAlgebra.Adjoint{T,A{T}},))
    end
  end
end

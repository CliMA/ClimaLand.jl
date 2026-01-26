module ArrayLayoutsSparseArraysExt

using ArrayLayouts
using ArrayLayouts: copyto!_layout
using SparseArrays
# Specifying the full namespace is necessary because of https://github.com/JuliaLang/julia/issues/48533
# See https://github.com/JuliaStats/LogExpFunctions.jl/pull/63
import ArrayLayouts.LinearAlgebra

import Base: copyto!

# ambiguity from sparsematrix.jl
copyto!(dest::LayoutMatrix, src::SparseArrays.AbstractSparseMatrixCSC) =
	copyto!_layout(dest, src)

copyto!(dest::SubArray{<:Any,2,<:LayoutMatrix}, src::SparseArrays.AbstractSparseMatrixCSC) =
	copyto!_layout(dest, src)

@inline LinearAlgebra.dot(a::LayoutArray{<:Number}, b::SparseArrays.SparseVectorUnion{<:Number}) =
	ArrayLayouts.dot(a,b)

@inline LinearAlgebra.dot(a::SparseArrays.SparseVectorUnion{<:Number}, b::LayoutArray{<:Number}) =
	ArrayLayouts.dot(a,b)

end

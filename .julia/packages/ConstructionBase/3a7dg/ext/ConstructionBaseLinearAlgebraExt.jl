module ConstructionBaseLinearAlgebraExt

import ConstructionBase
import LinearAlgebra

### Tridiagonal

function tridiagonal_constructor(dl::V, d::V, du::V) where {V<:AbstractVector{T}} where T
    LinearAlgebra.Tridiagonal{T,V}(dl, d, du)
end
function tridiagonal_constructor(dl::V, d::V, du::V, du2::V) where {V<:AbstractVector{T}} where T
    LinearAlgebra.Tridiagonal{T,V}(dl, d, du, du2)
end

# `du2` may be undefined, so we need a custom `getfields` that checks `isdefined`
function ConstructionBase.getfields(o::LinearAlgebra.Tridiagonal)
    if isdefined(o, :du2)
        (dl=o.dl, d=o.d, du=o.du, du2=o.du2)
    else
        (dl=o.dl, d=o.d, du=o.du)
    end
end

ConstructionBase.constructorof(::Type{<:LinearAlgebra.Tridiagonal}) = tridiagonal_constructor

### Cholesky

ConstructionBase.setproperties(C::LinearAlgebra.Cholesky, patch::NamedTuple{()}) = C

function ConstructionBase.setproperties(C::LinearAlgebra.Cholesky, patch::NamedTuple{(:L,),<:Tuple{<:LinearAlgebra.LowerTriangular}})
    return LinearAlgebra.Cholesky(C.uplo === 'U' ? copy(patch.L.data') : patch.L.data, C.uplo, C.info)
end
function ConstructionBase.setproperties(C::LinearAlgebra.Cholesky, patch::NamedTuple{(:U,),<:Tuple{<:LinearAlgebra.UpperTriangular}})
    return LinearAlgebra.Cholesky(C.uplo === 'L' ? copy(patch.U.data') : patch.U.data, C.uplo, C.info)
end
function ConstructionBase.setproperties(
    C::LinearAlgebra.Cholesky,
    patch::NamedTuple{(:UL,),<:Tuple{<:Union{LinearAlgebra.LowerTriangular,LinearAlgebra.UpperTriangular}}}
)
    return LinearAlgebra.Cholesky(patch.UL.data, C.uplo, C.info)
end
function ConstructionBase.setproperties(C::LinearAlgebra.Cholesky, patch::NamedTuple)
    throw(ArgumentError("Invalid patch for `Cholesky`: $(patch)"))
end

end #module
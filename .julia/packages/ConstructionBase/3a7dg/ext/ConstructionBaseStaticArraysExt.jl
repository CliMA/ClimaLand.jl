module ConstructionBaseStaticArraysExt

using ConstructionBase
using StaticArrays

# general static arrays need to keep the size parameter
ConstructionBase.constructorof(sa::Type{<:SArray{S}}) where {S} = SArray{S}
ConstructionBase.constructorof(sa::Type{<:MArray{S}}) where {S} = MArray{S}
ConstructionBase.constructorof(sa::Type{<:SizedArray{S}}) where {S} = SizedArray{S}

# static vectors don't even need the explicit size specification
ConstructionBase.constructorof(::Type{<:SVector}) = SVector
ConstructionBase.constructorof(::Type{<:MVector}) = MVector

# set properties by name: x, y, z, w
@generated function ConstructionBase.setproperties(obj::Union{SVector{N}, MVector{N}}, patch::NamedTuple{KS}) where {N, KS}
    if KS == (:data,)
        :( constructorof(typeof(obj))(only(patch)) )
    else
        N <= 4 || error("type $obj does not have properties $(join(KS, ", "))")
        propnames = (:x, :y, :z, :w)[1:N]
        KS ⊆ propnames || error("type $obj does not have properties $(join(KS, ", "))")
        field_exprs = map(enumerate(propnames)) do (i, p)
            from = p ∈ KS ? :patch : :obj
            :( $from.$p )
        end
        :( constructorof(typeof(obj))($(field_exprs...)) )
    end
end

end

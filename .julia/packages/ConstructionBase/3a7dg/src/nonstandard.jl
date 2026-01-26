
### SubArray
# `offset1` and `stride1` fields are calculated from parent indices.
# Setting them has no effect.
subarray_constructor(parent, indices, args...) = SubArray(parent, indices)

constructorof(::Type{<:SubArray}) = subarray_constructor

### ReinterpretArray
struct ReinterpretArrayConstructor{T} end
# `readable` and `writeable` are calculated from `T` and `parent`
# Setting them has no effect.
function (::ReinterpretArrayConstructor{T})(parent, args...) where T
    reinterpret(T, parent)
end

constructorof(::Type{<:Base.ReinterpretArray{T}}) where T = ReinterpretArrayConstructor{T}()

### PermutedDimsArray
struct PermutedDimsArrayConstructor{N,perm,iperm} end
# Parent must have the same N - it has to match length(perm)
function (::PermutedDimsArrayConstructor{N,perm,iperm})(parent::AA
) where {N,perm,iperm,AA<:AbstractArray{T,N}} where T
    PermutedDimsArray{T,N,perm,iperm,AA}(parent)
end

constructorof(::Type{<:PermutedDimsArray{<:Any,N,perm,iperm,<:Any}}) where {N,perm,iperm} =
    PermutedDimsArrayConstructor{N,perm,iperm}()



### LinRange
# `lendiv` is a calculated field
linrange_constructor(start, stop, len, lendiv) = LinRange(start, stop, len)

constructorof(::Type{<:LinRange}) = linrange_constructor

### Expr: args get splatted
# ::Expr annotation is to make it type-stable on Julia 1.3-, probably not necessary anymore
constructorof(::Type{<:Expr}) = (head, args) -> Expr(head, args...)::Expr

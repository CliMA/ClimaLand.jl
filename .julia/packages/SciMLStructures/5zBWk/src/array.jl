hasportion(::Tunable, ::AbstractArray) = true
hasportion(::Constants, ::AbstractArray) = false
hasportion(::Caches, ::AbstractArray) = false
hasportion(::Discrete, ::AbstractArray) = false
hasportion(::Initials, ::AbstractArray) = false

struct ArrayRepack{T}
    x::T
end
function (f::ArrayRepack)(A)
    @assert length(A) == prod(size(f.x))
    if has_trivial_array_constructor(typeof(f.x), A)
        restructure(f.x, A)
    else
        error("The original type $(typeof(f.x)) does not support the SciMLStructures interface via the AbstractArray `repack` rules. No method exists to take in a regular array and construct the parent type back. Please define the SciMLStructures interface for this type.")
    end
end

canonicalize(::Tunable, p::AbstractArray) = vec(p), ArrayRepack(p), true
canonicalize(::Constants, p::AbstractArray) = nothing, nothing, nothing
canonicalize(::Caches, p::AbstractArray) = nothing, nothing, nothing
canonicalize(::Discrete, p::AbstractArray) = nothing, nothing, nothing

isscimlstructure(::AbstractArray) = false
isscimlstructure(::AbstractArray{<:Number}) = true

function SciMLStructures.replace(
        ::SciMLStructures.Tunable, arr::AbstractArray, new_arr::AbstractArray)
    restructure(arr, new_arr)
end

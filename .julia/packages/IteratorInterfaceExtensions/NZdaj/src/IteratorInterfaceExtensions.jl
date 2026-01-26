module IteratorInterfaceExtensions

export getiterator, isiterable, IteratorSize2

isiterable(x::T) where {T} = Base.isiterable(T)

function getiterator(x)
    if !isiterable(x)
        error("Can't get iterator for non iterable source.")
    end
    return x
end

struct HasLengthAfterStart <: Base.IteratorSize end

IteratorSize2(x) = IteratorSize2(typeof(x))
IteratorSize2(::Type{T}) where {T} = Base.IteratorSize(T)

end # module

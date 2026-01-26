# rawview
const AAFixed{T<:FixedPoint,N} = AbstractArray{T,N}
function Base.showarg(io::IO, A::MappedArray{T,N,AA,typeof(reinterpret)}, toplevel=false) where {T<:Integer,N,AA<:AAFixed}
    print(io, "rawview(")
    Base.showarg(io, parent(A), false)
    print(io, ')')
    toplevel && print(io, " with eltype ", T)
end

# normedview
const AAInteger{T<:Integer,N} = AbstractArray{T,N}
function Base.showarg(io::IO, A::MappedArray{T,N,AA,F,typeof(reinterpret)}, toplevel=false) where {T<:FixedPoint,N,AA<:AAInteger,F}
    print(io, "normedview(")
    ColorTypes.showcoloranttype(io, T)
    print(io, ", ")
    Base.showarg(io, parent(A), false)
    print(io, ')')
    toplevel && print(io, " with eltype ", T)
end

function Base.showarg(io::IO, r::Base.ReinterpretArray{T}, toplevel) where {T<:Colorant}
    print(io, "reinterpret(")
    ColorTypes.colorant_string_with_eltype(io, T)
    print(io, ", ")
    Base.showarg(io, parent(r), false)
    print(io, ')')
end

function Base.showarg(io::IO, r::Base.ReinterpretArray{T}, toplevel) where {T<:FixedPoint}
    print(io, "reinterpret(")
    ColorTypes.showcoloranttype(io, T)
    print(io, ", ")
    Base.showarg(io, parent(r), false)
    print(io, ')')
end

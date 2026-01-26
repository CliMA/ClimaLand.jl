module StructUtilsMeasurementsExt

using StructUtils, Measurements

StructUtils.structlike(::Type{<:Measurements.Measurement}) = false
StructUtils.lower(m::Measurements.Measurement) = string(m)
StructUtils.lift(::Type{T}, m::String) where {T<:Measurements.Measurement} = parse(T, m)

end # module

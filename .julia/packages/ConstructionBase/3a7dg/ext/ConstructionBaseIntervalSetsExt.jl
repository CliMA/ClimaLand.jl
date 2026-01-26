module ConstructionBaseIntervalSetsExt

using ConstructionBase
using IntervalSets

ConstructionBase.constructorof(::Type{<:Interval{L, R}}) where {L, R} = Interval{L, R}

end

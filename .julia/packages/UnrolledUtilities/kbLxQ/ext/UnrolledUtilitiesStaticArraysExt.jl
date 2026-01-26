module UnrolledUtilitiesStaticArraysExt

import UnrolledUtilities
import StaticArrays: SVector, MVector

@inline UnrolledUtilities.output_type_for_promotion(::SVector) = SVector
@inline UnrolledUtilities.constructor_from_tuple(::Type{SVector}) = SVector

@inline UnrolledUtilities.output_type_for_promotion(::MVector) = MVector
@inline UnrolledUtilities.constructor_from_tuple(::Type{MVector}) = MVector

end

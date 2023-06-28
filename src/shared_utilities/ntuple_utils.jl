"""
### Copied from the Clima EDMF codebase ###

    Cent{I <: Integer}

A helper struct which is used in setting indices of ClimaCore
Fields of NTuples.
"""
struct Cent{I <: Integer}
    i::I
end

"""
### Copied from the Clima EDMF codebase ###

    Base.setindex!(field::ClimaCore.Fields.Field,
                   v,
                   i::Cent,
                  ) 

Sets a component of a ClimaCore Field of Ntuples by indexing.
"""
Base.@propagate_inbounds Base.setindex!(
    field::ClimaCore.Fields.Field,
    v,
    i::Cent,
) = Base.setindex!(ClimaCore.Fields.field_values(field), v, i.i)

"""
### Copied from the Clima EDMF codebase ###

    Base.getindex(field::ClimaCore.Fields.Field,
                   i::Integer,
                  ) 

Returns the element of the ClimaCore.Fields.Field by index.
"""
Base.@propagate_inbounds Base.getindex(
    field::ClimaCore.Fields.Field,
    i::Integer,
) = Base.getproperty(field, i)

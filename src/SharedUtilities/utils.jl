"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function.
"""
function heaviside(x::FT)::FT where {FT}
    if x >= FT(0.0)
        return FT(1.0)
    else
        return FT(0.0)
    end
end

"""
    to_scalar_coefs(vector_coefs)

Helper function used in computing tendencies of vertical diffusion terms.
"""
to_scalar_coefs(vector_coefs) = map(vector_coef -> vector_coef.u₃, vector_coefs)

"""
   dss!(dY::ClimaCore.Fields.FieldVector, p::ClimaCore.Fields.FieldVector, T::FT)

Computes the weighted direct stiffness summation and updates `dY` in place.
In the case of a column domain, no dss operations are performed.
"""
function dss!(
    dY::ClimaCore.Fields.FieldVector,
    _::ClimaCore.Fields.FieldVector,
    _::FT,
) where {FT}
    for key in propertynames(dY)
        property = getproperty(dY, key)
        dss_helper!(property, axes(property))
    end
end

"""
    dss_helper!(field_vec::ClimaCore.Fields.FieldVector, _)

Method of `dss_helper!` which unpacks properties of Y when on a
domain that is 2-dimensional in the horizontal.

The assumption is that Y contains FieldVectors which themselves contain either
FieldVectors or Fields, and that the final unpacked variable is a Field.
This method is invoked when the current property itself contains additional
property(ies).
"""
function dss_helper!(field_vec::ClimaCore.Fields.FieldVector, _)
    for key in propertynames(field_vec)
        property = getproperty(field_vec, key)
        dss_helper!(property, axes(property))
    end
end

"""
    dss_helper!(field::ClimaCore.Fields.Field,
            domain::Union{ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace, ClimaCore.Spaces.AbstractSpectralElementSpace})

Method of `dss_helper!` which performs dss on fields of Y when on a
domain that is 2-dimensional in the horizontal.

The assumption is that Y contains FieldVectors which themselves contain either
FieldVectors or Fields, and that the final unpacked variable is a Field.
This method is invoked when the element cannot be unpacked further.
"""
function dss_helper!(
    field::ClimaCore.Fields.Field,
    domain::Union{
        ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.AbstractSpectralElementSpace,
    },
)
    Spaces.weighted_dss!(field)
end

"""
    dss_helper!(_,
        domain::Union{ClimaCore.Spaces.FiniteDifferenceSpace, ClimaCore.Spaces.PointSpace})

Computes the appropriate weighted direct stiffness summation based on
the domain type, updates `Y` in place.

For column (FiniteDifferenceSpace) domains, no dss is needed.
"""
function dss_helper!(
    _,
    domain::Union{
        ClimaCore.Spaces.FiniteDifferenceSpace,
        ClimaCore.Spaces.PointSpace,
    },
) end

# TODO this function is used for canopy conductance/photosynthesis/rad transfer in ozark - is this the best way?
function dss_helper!(
    vector::Vector,
    domain::Tuple,
) end

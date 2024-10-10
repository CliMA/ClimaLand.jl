export boundary_flux,
    AbstractBC,
    AbstractBoundary,
    TopBoundary,
    BottomBoundary,
    diffusive_flux,
    boundary_var_domain_names,
    boundary_var_types,
    boundary_vars


### These methods and types are not required for all models. They
### are only useful for models using ClimaCore's differential operators
### for solving PDEs. Your model may not need them.
## Currently, only the soil and soil CO2 models use these.
"""
    AbstractBC

An abstract type for types of boundary conditions, which will include
prescribed functions of space and time as Dirichlet conditions or
Neumann conditions, in addition to other 
convenient conditions.

"""
abstract type AbstractBC end

"""
    AbstractBoundary

An abstract type to indicate which boundary we are doing calculations for.
Currently, we support the top boundary (TopBoundary)
and bottom boundary (BottomBoundary).
"""
abstract type AbstractBoundary end


"""
    TopBoundary{} <: AbstractBoundary{}

A simple object which should be passed into a function to
indicate that we are considering the top boundary.
"""
struct TopBoundary <: AbstractBoundary end

"""
    BottomBoundary{} <: AbstractBoundary{}

A simple object which should be passed into a function to
indicate that we are considering the bottom boundary.
"""
struct BottomBoundary <: AbstractBoundary end


bc_name(::BottomBoundary) = :bottom_bc
bc_name(::TopBoundary) = :top_bc

"""
    diffusive_flux(K, x_2, x_1, Δz)

Calculates the diffusive flux of a quantity x (water content, temp, etc).
Here, x_2 = x(z + Δz) and x_1 = x(z), so x_2 is at a larger z by convention.
"""
function diffusive_flux(K, x_2, x_1, Δz)
    return @. -K * (x_2 - x_1) / Δz
end

"""
    boundary_flux(bc::AbstractBC, bound_type::AbstractBoundary, Δz, _...)::ClimaCore.Fields.Field
A function which returns the correct boundary flux  given
    any boundary condition (BC). 
"""
function boundary_flux(
    bc::AbstractBC,
    bound_type::AbstractBoundary,
    Δz,
    _...,
)::ClimaCore.Fields.Field end


"""
    boundary_vars(::AbstractBC , ::ClimaLand.TopBoundary)

The list of symbols for additional variables to add to the
model auxiliary state, for models solving PDEs, which defaults 
to adding storage for the top boundary flux fields,
but which can be extended depending on the type of 
boundary condition used.

For the Soil and SoilCO2 models - which solve PDEs - the 
tendency functions and 
update_boundary_fluxes functions are coded to access the field 
`:top_bc`  to be present in the model cache, which is why 
this is the default.  If this is not your (PDE) model's 
desired behavior, you can extend this function with a new method.

The field `top_bc_wvec` is created to prevent allocations only; it is used
in the tendency function only.

Use this function in the exact same way you would use `auxiliary_vars`.
"""
boundary_vars(::AbstractBC, ::ClimaLand.TopBoundary) = (:top_bc, :top_bc_wvec)

"""
    boundary_vars(::AbstractBC, ::ClimaLand.BottomBoundary)

The list of symbols for additional variables to add to the
model auxiliary state, for models solving PDEs, which defaults 
to adding storage for the bottom boundary flux fields,
but which can be extended depending on the type of 
boundary condition used.

For the Soil and SoilCO2 models - which solve PDEs - the 
tendency functions and 
update_boundary_fluxes functions are coded to access the field 
`:bottom_bc`  to be present in the model cache, which is why 
this is the default.  If this is not your (PDE) model's 
desired behavior, you can extend this function with a new method.

The field `bottom_bc_wvec` is created to prevent allocations only; it is used
in the tendency function only.

Use this function in the exact same way you would use `auxiliary_vars`.
"""
boundary_vars(::AbstractBC, ::ClimaLand.BottomBoundary) =
    (:bottom_bc, :bottom_bc_wvec)

"""
    boundary_var_domain_names(::AbstractBC, ::ClimaLand.AbstractBoundary)

The list of domain names for additional variables to add to the
model auxiliary state, for models solving PDEs, which defaults 
to adding storage on the surface domain 
for the top or bottom boundary flux fields,
but which can be extended depending on the type of 
boundary condition used.

Use in conjunction with `boundary_vars`, in the same way you would use
`auxiliary_var_domain_names`. 
"""
boundary_var_domain_names(::AbstractBC, ::ClimaLand.AbstractBoundary) =
    (:surface, :surface)

"""
    boundary_var_types(model::AbstractModel{FT}, ::AbstractBC, ::ClimaLand.AbstractBoundary) where {FT}

The list of types for additional variables to add to the
model auxiliary state, for models solving PDEs, which defaults 
to adding a scalar variable on the surface domain 
for the top or bottom boundary flux fields,
but which can be extended depending on the type of 
boundary condition used.

Use in conjunction with `boundary_vars`, in the same way you would use
`auxiliary_var_types`. The use of a scalar is appropriate for
models with a single PDE; models with multiple PDEs will need to supply
multiple scalar fields.
"""
function boundary_var_types(
    model::AbstractModel{FT},
    ::AbstractBC,
    ::ClimaLand.AbstractBoundary,
) where {FT}
    (FT, ClimaCore.Geometry.WVector{FT})
end

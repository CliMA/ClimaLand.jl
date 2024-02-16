import ClimaLand: AbstractBC, AbstractBoundary, boundary_flux
export TemperatureStateBC,
    MoistureStateBC,
    FreeDrainage,
    FluxBC,
    AtmosDrivenFluxBC,
    RichardsAtmosDrivenFluxBC,
    boundary_vars,
    boundary_var_domain_names,
    boundary_var_types

"""
    get_top_surface_field(
        center_field::ClimaCore.Fields.Field,
        surface_space,
    )

A helper function to extract the top level of a center field and
cast it onto the surface face space.
"""
function get_top_surface_field(
    center_field::ClimaCore.Fields.Field,
    surface_space,
)
    # TODO: Find cleaner way to do this
    nz = Spaces.nlevels(axes(center_field))
    return Fields.Field(
        Fields.field_values(Fields.level(center_field, nz)),
        surface_space,
    )
end

"""
    get_top_surface_field(
        center_val,
        _,
    )

A helper function for the case where we use `get_top_surface_field` on
a parameter that is a scalar rather than a field. Returns the scalar value.
"""
function get_top_surface_field(center_val, _)
    return center_val
end

"""
    get_bottom_surface_field(
        center_field::ClimaCore.Fields.Field,
        bottom_space,
    )

A helper function to extract the bottom level of a center field and
cast it onto the bottom face space.
"""
function get_bottom_surface_field(
    center_field::ClimaCore.Fields.Field,
    bottom_space,
)
    # TODO: Find cleaner way to do this
    return Fields.Field(
        Fields.field_values(Fields.level(center_field, 1)),
        bottom_space,
    )
end

"""
    get_bottom_surface_field(
        center_val,
        _,
    )

A helper function for the case where we use `get_bottom_surface_field` on
a parameter that is a scalar rather than a field. Returns the scalar value.
"""
function get_bottom_surface_field(center_val, _)
    return center_val
end

"""
    AbstractSoilBC <: ClimaLand. AbstractBC

An abstract type for soil-specific types of boundary conditions, like free drainage.
"""
abstract type AbstractSoilBC <: ClimaLand.AbstractBC end


"""
    boundary_vars(::NamedTuple, ::ClimaLand.TopBoundary)

The list of symbols for additional variables to add to the soil
model auxiliary state, which defaults to adding storage for the
 top boundary flux fields,
but which can be extended depending on the type of boundary condition used.

Note that `:top_bc` must be present, with the default type and domain name,
for both the RichardsModel and the 
EnergyHydrology soil models.

 Use this function in the exact same way you would use `auxiliary_vars`.
"""
boundary_vars(::NamedTuple, ::ClimaLand.TopBoundary) = (:top_bc,)

"""
    boundary_vars(::NamedTuple, ::ClimaLand.BottomBoundary)

The list of symbols for additional variables to add to the soil
model auxiliary state, which defaults to adding storage for the
 bottom boundary flux fields,
but which can be extended depending on the type of boundary condition used.

Use this function in the exact same way you would use `auxiliary_vars`.

Note that `:bottom_bc` must be present, with the default type and domain name,
for both the RichardsModel and the 
EnergyHydrology soil models.
"""
boundary_vars(::NamedTuple, ::ClimaLand.BottomBoundary) = (:bottom_bc,)

"""
    boundary_var_domain_names(::AbstractSoilBC, ::ClimaLand.AbstractBoundary)

The list of domain names for additional variables to add to the soil
model auxiliary state,  which defaults to adding storage for the
 boundary flux fields,
but which can be extended depending on the type of boundary condition used.
Note that extensions to this function must still include a flux bc defined on
the surface in addition to other variables.

Use this function in the exact same way you would use `auxiliary_var_domain_names`.
"""
boundary_var_domain_names(::NamedTuple, ::ClimaLand.AbstractBoundary) =
    (:surface,)

"""
    boundary_var_types(::Soil.EnergyHydrology{FT}, ::NamedTuple, ::ClimaLand.AbstractBoundary) where {FT}

The list of variable types for additional variables to add to the EnergyHydrology
model auxiliary state, which defaults to adding storage for the
 boundary flux fields,
but which can be extended depending on the type of boundary condition used.

Use this function in the exact same way you would use `auxiliary_types`.
Note that extensions to this function must still include a flux bc defined on
the surface in addition to other variables.
"""
function boundary_var_types(
    ::Soil.EnergyHydrology{FT},
    ::NamedTuple,
    ::ClimaLand.AbstractBoundary,
) where {FT}
    (NamedTuple{(:water, :heat), Tuple{FT, FT}},)
end

"""
    boundary_var_types(::Soil.RichardsModel{FT}, ::NamedTuple, ::ClimaLand.AbstractBoundary) where {FT}

The list of variable types for additional variables to add to the Richards
model auxiliary state,  which defaults to adding storage for the
 boundary flux fields,
but which can be extended depending on the type of boundary condition used.

Since the only prognostic variable for the Richards-Richardson equation is the
volumetric water content, only a water flux boundary condition is required per boundary.

Use this function in the exact same way you would use `auxiliary_types`.
Note that extensions to this function must still include a flux bc defined on
the surface in addition to other variables.
"""
boundary_var_types(
    ::Soil.RichardsModel{FT},
    ::NamedTuple,
    ::ClimaLand.AbstractBoundary,
) where {FT} = (NamedTuple{(:water,), Tuple{FT}},)


"""
   MoistureStateBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
state boundary condition ϑ_l = f(p,t) at either the top or bottom of the domain.
"""
struct MoistureStateBC{F <: Function} <: AbstractSoilBC
    bc::F
end

"""
   TemperatureStateBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
state boundary condition T = f(p,t) at either the top or bottom of the domain.
"""
struct TemperatureStateBC{F <: Function} <: AbstractSoilBC
    bc::F
end

"""
   FluxBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
normal flux boundary condition f(p,t) at either the top or bottom of the domain.
"""
struct FluxBC{F <: Function} <: AbstractSoilBC
    bc::F
end

"""
    FreeDrainage <: AbstractSoilBC
A concrete type of soil boundary condition, for use at
the BottomBoundary only, where the flux is set to be
`F = -K∇h = -K`.
"""
struct FreeDrainage <: AbstractSoilBC end


"""
   RichardsAtmosDrivenFluxBC{F <: Function, R <: AbstractRunoffModel} <: AbstractSoilBC

A concrete type of boundary condition intended only for use with the RichardsModel,
which uses a prescribed precipitation rate (m/s) to compute the infiltration
into the soil.

A runoff model is used
to simulate surface and subsurface runoff and this is accounted
for when setting boundary conditions. In order to run the simulation
*without* runoff, choose runoff = NoRunoff() - this is also the default.

If you wish to simulate precipitation and runoff in the full `EnergyHydrology` model,
you must use the `AtmosDrivenFluxBC` type.
$(DocStringExtensions.FIELDS)
"""
struct RichardsAtmosDrivenFluxBC{F <: Function, R <: AbstractRunoffModel} <:
       AbstractSoilBC
    "The prescribed liquid water precipitation rate f(t) (m/s); Negative by convention."
    precip::F
    "The runoff model. The default is no runoff."
    runoff::R
end

function RichardsAtmosDrivenFluxBC(precip::Function; runoff = NoRunoff())
    return RichardsAtmosDrivenFluxBC{typeof(precip), typeof(runoff)}(
        precip,
        runoff,
    )
end



"""
    AtmosDrivenFluxBC{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
        R <: AbstractRunoffModel
    } <: AbstractSoilBC

A concrete type of soil boundary condition for use at the top
of the domain. This holds the conditions for the atmosphere
`AbstractAtmosphericDrivers`, for the radiation state
`AbstractRadiativeDrivers`. This is only supported for the
`EnergyHydrology` model.

This choice indicates the Monin-Obukhov Surface Theory will
be used to compute the sensible and latent heat fluxes, as
well as evaporation,
 and that the net radiation and precipitation
will also be computed. The net energy and water fluxes
are used as boundary conditions.

A runoff model is used
to simulate surface and subsurface runoff and this is accounted
for when setting boundary conditions. The default is to have no runoff
accounted for.

$(DocStringExtensions.FIELDS)
"""
struct AtmosDrivenFluxBC{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
    R <: AbstractRunoffModel,
} <: AbstractSoilBC
    "The atmospheric conditions driving the model"
    atmos::A
    "The radiative fluxes driving the model"
    radiation::B
    "The runoff model. The default is no runoff."
    runoff::R
end


function AtmosDrivenFluxBC(atmos, radiation; runoff = NoRunoff())
    args = (atmos, radiation, runoff)
    return AtmosDrivenFluxBC{typeof.(args)...}(args...)
end

"""
    create_soil_bc_named_tuple(water_flux, heat_flux)

A helper function which takes in two scalar values of `water_flux`
and `heat_flux`, and creates a named tuple out of them.

When broadcasted over two ClimaCore.Fields.Field objects, 
this returns a Field of NamedTuples which we can access like
x.water, x.heat, to obtain the boundary condition fields.
"""
function create_soil_bc_named_tuple(water_flux, heat_flux)
    return (; :water => water_flux, :heat => heat_flux)
end

"""
    boundary_vars(::AtmosDrivenFluxBC, ::ClimaLand.TopBoundary)

An extension of the `boundary_vars` method for AtmosDrivenFluxBC. This
adds the surface conditions (SHF, LHF, evaporation, and resistance) and the
net radiation to the auxiliary variables.

These variables are updated in place in `soil_boundary_fluxes`.
"""
boundary_vars(bc::AtmosDrivenFluxBC, ::ClimaLand.TopBoundary) =
    (:turbulent_fluxes, :R_n, :top_bc)

"""
    boundary_var_domain_names(::AtmosDrivenFluxBC, ::ClimaLand.TopBoundary))

An extension of the `boundary_var_domain_names` method for AtmosDrivenFluxBC. This
specifies the part of the domain on which the additional variables should be
defined.
"""
boundary_var_domain_names(bc::AtmosDrivenFluxBC, ::ClimaLand.TopBoundary) =
    (:surface, :surface, :surface)
"""
    boundary_var_types(
        ::AtmosDrivenFluxBC{
            <:PrescribedAtmosphere{FT},
            <:AbstractRadiativeDrivers{FT},
            <:AbstractRunoffModel,
        }, ::ClimaLand.TopBoundary,
    ) where {FT}

An extension of the `boundary_var_types` method for AtmosDrivenFluxBC. This
specifies the type of the additional variables.
"""
boundary_var_types(
    model::EnergyHydrology{FT},
    bc::AtmosDrivenFluxBC{
        <:PrescribedAtmosphere{FT},
        <:AbstractRadiativeDrivers{FT},
        <:AbstractRunoffModel,
    },
    ::ClimaLand.TopBoundary,
) where {FT} = (
    NamedTuple{(:lhf, :shf, :vapor_flux, :r_ae), Tuple{FT, FT, FT, FT}},
    FT,
    NamedTuple{(:water, :heat), Tuple{FT, FT}},
)


"""
    soil_boundary_fluxes(
        bc::AtmosDrivenFluxBC{
            <:PrescribedAtmosphere,
            <:PrescribedRadiativeFluxes,
        },
        boundary::ClimaLand.TopBoundary,
        model::EnergyHydrology,
        Δz,
        Y,
        p,
        t,
    )

Returns the net volumetric water flux (m/s) and net energy
flux (W/m^2) for the soil `EnergyHydrology` model at the top
of the soil domain.

If you wish to compute surface fluxes taking into account the
presence of a canopy, snow, etc, as in a land surface model,
this is not the correct method to be using.

This function calls the `turbulent_fluxes` and `net_radiation`
functions, which use the soil surface conditions as well as
the atmos and radiation conditions in order to
compute the surface fluxes using Monin Obukhov Surface Theory.
"""
function soil_boundary_fluxes(
    bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:PrescribedRadiativeFluxes},
    boundary::ClimaLand.TopBoundary,
    model::EnergyHydrology,
    Δz,
    Y,
    p,
    t,
)

    p.soil.turbulent_fluxes .= turbulent_fluxes(bc.atmos, model, Y, p, t)
    p.soil.R_n .= net_radiation(bc.radiation, model, Y, p, t)
    # We are ignoring sublimation for now
    precip = p.drivers.P_liq
    infiltration = soil_surface_infiltration(
        bc.runoff,
        precip .+ p.soil.turbulent_fluxes.vapor_flux,
        Y,
        p,
        model.parameters,
    )
    # We do not model the energy flux from infiltration
    return @. create_soil_bc_named_tuple(
        infiltration,
        p.soil.R_n + p.soil.turbulent_fluxes.lhf + p.soil.turbulent_fluxes.shf,
    )
end

"""
    ClimaLand.boundary_flux(bc::FluxBC,  _...)::ClimaCore.Fields.Field

A method of boundary fluxes which returns the desired flux.

We add a field of zeros in order to convert the bc (float) into
a field.
"""
function ClimaLand.boundary_flux(
    bc::FluxBC,
    boundary::ClimaLand.AbstractBoundary,
    model,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    FT = eltype(Δz)
    return FT.(bc.bc(p, t)) .+ FT.(ClimaCore.Fields.zeros(axes(Δz)))
end


"""
    ClimaLand.boundary_flux(bc::RichardsAtmosDrivenFluxBC,
                           boundary::ClimaLand.AbstractBoundary,
                           model::RichardsModel{FT},
                           Δz::ClimaCore.Fields.Field,
                           Y::ClimaCore.Fields.FieldVector,
                           p::NamedTuple,
                           t,
                           )::ClimaCore.Fields.Field where {FT}

A method of boundary fluxes which returns the desired water volume flux for
the RichardsModel, at the top of the domain, in the case of a prescribed
precipitation flux.

If `model.runoff` is not of type `NoRunoff`, surface runoff is accounted for
when computing the infiltration.
"""
function ClimaLand.boundary_flux(
    bc::RichardsAtmosDrivenFluxBC,
    boundary::ClimaLand.AbstractBoundary,
    model::RichardsModel,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    FT = eltype(Δz)
    precip = FT.(bc.precip(t)) .+ FT.(ClimaCore.Fields.zeros(axes(Δz)))
    return soil_surface_infiltration(bc.runoff, precip, Y, p, model.parameters)
end

"""
    ClimaLand.boundary_flux(rre_bc::MoistureStateBC,
                           ::ClimaLand.TopBoundary,
                           model::AbstractSoilModel,
                           Δz::ClimaCore.Fields.Field,
                           Y::ClimaCore.Fields.FieldVector,
                           p::NamedTuple,
                           t,
                           )::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on θ_l at the top of the
domain into a flux of liquid water.
"""
function ClimaLand.boundary_flux(
    rre_bc::MoistureStateBC,
    ::ClimaLand.TopBoundary,
    model::AbstractSoilModel,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    FT = eltype(Δz)
    surface_space = axes(Δz)
    # First extract the value of the top layer of the pressure head
    # and cast onto the face space
    ψ_c = get_top_surface_field(p.soil.ψ, surface_space)

    # Calculate pressure head using boundary condition on ϑ_l = θ_bc
    # We first need to extract the parameters of the soil in the top layer
    # Again, we need to cast them onto the face space
    (; hydrology_cm, θ_r, ν, S_s) = model.parameters
    hcm_bc = get_top_surface_field(hydrology_cm, surface_space)
    θ_r_bc = get_top_surface_field(θ_r, surface_space)
    ν_bc = get_top_surface_field(ν, surface_space)
    S_s_bc = get_top_surface_field(S_s, surface_space)

    θ_bc = FT.(rre_bc.bc(p, t))
    ψ_bc = @. pressure_head(hcm_bc, θ_r_bc, θ_bc, ν_bc, S_s_bc)

    # Lastly, we need an effective conductivity to compute the flux that results from
    # the gradient in pressure.
    # currently we approximate this as equal to the center value at the top layer (K_c)
    # More accurate would be to compute the mean between K_c and K evaluated at the boundary
    # condition.
    K_eff = get_top_surface_field(p.soil.K, surface_space)

    # Pass in (ψ_bc .+ Δz) as to account for contribution of gravity (∂(ψ+z)/∂z
    return ClimaLand.diffusive_flux(K_eff, ψ_bc .+ Δz, ψ_c, Δz)
end

"""
    ClimaLand.boundary_flux(rre_bc::MoistureStateBC,
                           ::ClimaLand.BottomBoundary,
                           model::AbstractSoilModel,
                           Δz::ClimaCore.Fields.Field,
                           Y::ClimaCore.Fields.FieldVector,
                           p::NamedTuple,
                           t,
                           )::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on θ_l at the bottom of the
domain into a flux of liquid water.
"""
function ClimaLand.boundary_flux(
    rre_bc::MoistureStateBC,
    ::ClimaLand.BottomBoundary,
    model::AbstractSoilModel,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    FT = eltype(Δz)
    surface_space = axes(Δz)
    # First extract the value of the bottom layer of the pressure head
    # and cast onto the face space
    ψ_c = get_bottom_surface_field(p.soil.ψ, surface_space)


    # Calculate pressure head using boundary condition on ϑ_l = θ_bc
    # We first need to extract the parameters of the soil in the bottom layer
    # Again, we need to cast them onto the face space
    (; hydrology_cm, θ_r, ν, S_s) = model.parameters
    hcm_bc = get_bottom_surface_field(hydrology_cm, surface_space)
    θ_r_bc = get_bottom_surface_field(θ_r, surface_space)
    ν_bc = get_bottom_surface_field(ν, surface_space)
    S_s_bc = get_bottom_surface_field(S_s, surface_space)

    θ_bc = FT.(rre_bc.bc(p, t))
    ψ_bc = @. pressure_head(hcm_bc, θ_r_bc, θ_bc, ν_bc, S_s_bc)

    # Lastly, we need an effective conductivity to compute the flux that results from
    # the gradient in pressure.
    # currently we approximate this as equal to the center value at the top layer (K_c)
    # More accurate would be to compute the mean between K_c and K evaluated at the boundary
    # condition.
    K_eff = get_bottom_surface_field(p.soil.K, surface_space)

    # At the bottom boundary, ψ_c is at larger z than ψ_bc
    #  so we swap their order in the derivative calc
    # Pass in (ψ_c .+ Δz)  to account for contribution of gravity  (∂(ψ+z)/∂z
    return ClimaLand.diffusive_flux(K_eff, ψ_c .+ Δz, ψ_bc, Δz)
end


"""
    ClimaLand.boundary_flux(heat_bc::TemperatureStateBC,
                           ::ClimaLand.TopBoundary,
                           model::EnergyHydrology,
                           Δz::ClimaCore.Fields.Field,
                           Y::ClimaCore.Fields.FieldVector,
                           p::NamedTuple,
                           t,
                           )::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on temperature at the top of the
domain into a flux of energy.
"""
function ClimaLand.boundary_flux(
    heat_bc::TemperatureStateBC,
    ::ClimaLand.TopBoundary,
    model::EnergyHydrology,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    FT = eltype(Δz)
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    surface_space = axes(Δz)
    # We need to project the center values onto the face space.
    T_c = get_top_surface_space(p.soil.T, surface_space)
    κ_c = get_top_surface_space(p.soil.κ, surface_space)

    T_bc = FT.(heat_bc.bc(p, t))
    return ClimaLand.diffusive_flux(κ_c, T_bc, T_c, Δz)
end

"""
    ClimaLand.boundary_flux(heat_bc::TemperatureStateBC,
                           ::ClimaLand.BottomBoundary,
                           model::EnergyHydrology,
                           Δz::ClimaCore.Fields.Field,
                           Y::ClimaCore.Fields.FieldVector,
                           p::NamedTuple,
                           t,
                           )::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on temperature at the bottom of the
domain into a flux of energy.
"""
function ClimaLand.boundary_flux(
    heat_bc::TemperatureStateBC,
    ::ClimaLand.BottomBoundary,
    model::EnergyHydrology,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    FT = eltype(Δz)
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    surface_space = axes(Δz)
    # We need to project the center values onto the face space.
    T_c = get_bottom_surface_space(p.soil.T, surface_space)
    κ_c = get_bottom_surface_space(p.soil.κ, surface_space)
    T_bc = FT.(heat_bc.bc(p, t))
    return ClimaLand.diffusive_flux(κ_c, T_c, T_bc, Δz)
end


"""
    ClimaLand.boundary_flux(bc::FreeDrainage,
                           boundary::ClimaLand.BottomBoundary,
                           model::AbstractSoilModel,
                           Δz::ClimaCore.Fields.Field,
                           Y::ClimaCore.Fields.FieldVector,
                           p::NamedTuple,
                           t,
                           )::ClimaCore.Fields.Field

A method of boundary fluxes which enforces free drainage at the bottom
of the domain.
"""
function ClimaLand.boundary_flux(
    bc::FreeDrainage,
    boundary::ClimaLand.BottomBoundary,
    model::AbstractSoilModel,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    FT = eltype(Δz)
    K_c = Fields.level(p.soil.K, 1)
    return FT.(-1 .* K_c)
end

"""
    soil_boundary_fluxes(bc::NamedTuple, boundary, model, Δz, Y, p, t)

Returns the boundary fluxes for ϑ_l and ρe_int, in that order.
"""
function soil_boundary_fluxes(bc::NamedTuple, boundary, model, Δz, Y, p, t)
    params = model.parameters
    return create_soil_bc_named_tuple.(
        ClimaLand.boundary_flux(bc.water, boundary, model, Δz, Y, p, t),
        ClimaLand.boundary_flux(bc.heat, boundary, model, Δz, Y, p, t),
    )
end

"""
    ClimaLand.∂tendencyBC∂Y(
        model::RichardsModel,
        ::MoistureStateBC,
        boundary::ClimaLand.TopBoundary,
        Δz,
        Y,
        p,
        t,
)

Computes and returns the derivative of the part of the
implicit tendency in the top layer, due to the boundary
condition, with respect to the state variable in the top layer.

For a diffusion equation like Richards equation with a single state
variable, this is given by
`∂T_N∂Y_N = [-∂/∂z(∂F_bc/∂Y_N)]_N`, where `N` indicates the top
layer cell index.
"""
function ClimaLand.∂tendencyBC∂Y(
    model::RichardsModel,
    ::MoistureStateBC,
    boundary::ClimaLand.TopBoundary,
    Δz,
    Y,
    p,
    t,
)
    (; ν, hydrology_cm, S_s, θ_r) = model.parameters
    fs = ClimaLand.Domains.obtain_face_space(model.domain.space.subsurface)
    face_len = ClimaCore.Utilities.PlusHalf(ClimaCore.Spaces.nlevels(fs) - 1)
    interpc2f_op = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    K = Fields.level(interpc2f_op.(p.soil.K), face_len)
    dψ = Fields.level(
        interpc2f_op.(dψdϑ.(hydrology_cm, Y.soil.ϑ_l, ν, θ_r, S_s)),
        face_len,
    )
    return ClimaCore.Fields.FieldVector(;
        :soil => (; :ϑ_l => @. -K / Δz * dψ / (2 * Δz)),
    )
end

"""
    ClimaLand.∂tendencyBC∂Y(
        ::AbstractSoilModel,
        ::AbstractSoilBC,
        boundary::ClimaLand.TopBoundary,
        Δz,
        Y,
        p,
        t,
)

A default method which computes and returns the zero for the
derivative of the part of the
implicit tendency in the top layer, due to the boundary
condition, with respect to the state variable in the top layer.

For a diffusion equation like Richards equation with a single state
variable, this is given by
`∂T_N∂Y_N = [-∂/∂z(∂F_bc/∂Y_N)]_N`, where `N` indicates the top
layer cell index.

If `F_bc` can be approximated as independent of `Y_N`, the derivative
is zero.
"""
function ClimaLand.∂tendencyBC∂Y(
    ::AbstractSoilModel,
    ::AbstractSoilBC,
    boundary::ClimaLand.TopBoundary,
    Δz,
    Y,
    p,
    t,
)
    return ClimaCore.Fields.FieldVector(; :soil => (; :ϑ_l => zeros(axes(Δz))))
end

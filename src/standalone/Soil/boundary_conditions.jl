import ClimaLand:
    AbstractBoundary,
    AbstractBC,
    boundary_flux,
    boundary_vars,
    boundary_var_domain_names,
    boundary_var_types
export TemperatureStateBC,
    MoistureStateBC,
    FreeDrainage,
    HeatFluxBC,
    WaterFluxBC,
    AtmosDrivenFluxBC,
    RichardsAtmosDrivenFluxBC,
    WaterHeatBC

# Helper functions
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

# New BC type for Richards Equation (AbstractWaterBC)
"""
    AbstractWaterBC <: ClimaLand.AbstractBC

An abstract type for boundary conditions for Richards equation.
"""
abstract type AbstractWaterBC <: ClimaLand.AbstractBC end

"""
   MoistureStateBC <: AbstractWaterBC

A simple concrete type of boundary condition, which enforces a
state boundary condition ϑ_l = f(p,t) at either the top or bottom of the domain.
"""
struct MoistureStateBC{F <: Function} <: AbstractWaterBC
    bc::F
end

"""
   WaterFluxBC <: AbstractWaterBC

A simple concrete type of boundary condition, which enforces a
normal flux boundary condition f(p,t) at either the top or bottom of the domain.
"""
struct WaterFluxBC{F <: Function} <: AbstractWaterBC
    bc::F
end

"""
    FreeDrainage <: AbstractWaterBC
A concrete type of soil boundary condition, for use at
the BottomBoundary only, where the flux is set to be
`F = -K∇h = -K`.
"""
struct FreeDrainage <: AbstractWaterBC end


"""
   RichardsAtmosDrivenFluxBC{F <: PrescribedPrecipitation, R <: AbstractRunoffModel} <: AbstractWaterBC

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
struct RichardsAtmosDrivenFluxBC{
    F <: PrescribedPrecipitation,
    R <: AbstractRunoffModel,
} <: AbstractWaterBC
    "The prescribed liquid water precipitation rate f(t) (m/s); Negative by convention."
    precip::F
    "The runoff model. The default is no runoff."
    runoff::R
end

function RichardsAtmosDrivenFluxBC(
    precip::PrescribedPrecipitation;
    runoff = NoRunoff(),
)
    return RichardsAtmosDrivenFluxBC{typeof(precip), typeof(runoff)}(
        precip,
        runoff,
    )
end

# Methods
"""
    boundary_vars(::RichardsAtmosDrivenFluxBC{<:PrescribedPrecipitation,
                                              <:Runoff.TOPMODELRunoff,
                                              }, ::ClimaLand.TopBoundary)

An extension of the `boundary_vars` method for RichardsAtmosDrivenFluxBC with 
TOPMODELRunoff.

These variables are updated in place in `boundary_flux`.
"""
boundary_vars(
    bc::RichardsAtmosDrivenFluxBC{
        <:PrescribedPrecipitation,
        <:Runoff.TOPMODELRunoff,
    },
    ::ClimaLand.TopBoundary,
) = (:top_bc, :h∇, :R_s, :R_ss, :infiltration)

"""
    boundary_var_domain_names(::RichardsAtmosDrivenFluxBC{<:PrescribedPrecipitation,
                                              <:Runoff.TOPMODELRunoff,
                                              },
                              ::ClimaLand.TopBoundary)

An extension of the `boundary_var_domain_names` method for RichardsAtmosDrivenFluxBC
with TOPMODELRunoff. 
"""
boundary_var_domain_names(
    bc::RichardsAtmosDrivenFluxBC{
        <:PrescribedPrecipitation,
        <:Runoff.TOPMODELRunoff,
    },
    ::ClimaLand.TopBoundary,
) = (:surface, :surface, :surface, :surface, :surface)
"""
    boundary_var_types(::RichardsAtmosDrivenFluxBC{<:PrescribedPrecipitation,
                                                   <:Runoff.TOPMODELRunoff{FT},
                                                  },
                       ::ClimaLand.TopBoundary,
                       ) where {FT}

An extension of the `boundary_var_types` method for RichardsAtmosDrivenFluxBC
with TOPMODELRunoff.
"""
boundary_var_types(
    model::RichardsModel{FT},
    bc::RichardsAtmosDrivenFluxBC{
        <:PrescribedPrecipitation,
        <:Runoff.TOPMODELRunoff{FT},
    },
    ::ClimaLand.TopBoundary,
) where {FT} = (FT, FT, FT, FT, FT)

"""
    boundary_vars(::RichardsAtmosDrivenFluxBC{<:PrescribedPrecipitation,
                                              <:Runoff.AbstractRunoffModel,
                                              }, ::ClimaLand.TopBoundary)

An extension of the `boundary_vars` method for RichardsAtmosDrivenFluxBC with 
no runoff modeled.

These variables are updated in place in `boundary_flux`.
"""
boundary_vars(
    bc::RichardsAtmosDrivenFluxBC{
        <:PrescribedPrecipitation,
        <:Runoff.AbstractRunoffModel,
    },
    ::ClimaLand.TopBoundary,
) = (:top_bc, :infiltration)

"""
    boundary_var_domain_names(::RichardsAtmosDrivenFluxBC{<:PrescribedPrecipitation,
                                              <:Runoff.AbstractRunoffModel,
                                              },
                              ::ClimaLand.TopBoundary)

An extension of the `boundary_var_domain_names` method for RichardsAtmosDrivenFluxBC
with no runoff modeled. 
"""
boundary_var_domain_names(
    bc::RichardsAtmosDrivenFluxBC{
        <:PrescribedPrecipitation,
        <:Runoff.AbstractRunoffModel,
    },
    ::ClimaLand.TopBoundary,
) = (:surface, :surface)
"""
    boundary_var_types(::RichardsAtmosDrivenFluxBC{<:PrescribedPrecipitation,
                                                   <:Runoff.AbstractRunoffModel,
                                                  },
                       ::ClimaLand.TopBoundary,
                       ) where {FT}

An extension of the `boundary_var_types` method for RichardsAtmosDrivenFluxBC
with no runoff modeled.
"""
boundary_var_types(
    model::RichardsModel{FT},
    bc::RichardsAtmosDrivenFluxBC{
        <:PrescribedPrecipitation,
        <:Runoff.AbstractRunoffModel,
    },
    ::ClimaLand.TopBoundary,
) where {FT} = (FT, FT)

"""
    boundary_flux(bc::WaterFluxBC,  _...)::ClimaCore.Fields.Field

A method of boundary fluxes which returns the desired flux.

We add a field of zeros in order to convert the bc (float) into
a field.
"""
function boundary_flux(
    bc::WaterFluxBC,
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
    boundary_flux(bc::RichardsAtmosDrivenFluxBC,
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
function boundary_flux(
    bc::RichardsAtmosDrivenFluxBC,
    boundary::ClimaLand.AbstractBoundary,
    model::RichardsModel,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)::ClimaCore.Fields.Field
    update_runoff!(p, bc.runoff, Y, t, model)
    return p.soil.infiltration
end

"""
    boundary_flux(rre_bc::MoistureStateBC,
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
function boundary_flux(
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
    boundary_flux(rre_bc::MoistureStateBC,
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
function boundary_flux(
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
    boundary_flux(bc::FreeDrainage,
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
function boundary_flux(
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
        ::RichardsModel,
        ::AbstractWaterBC,
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
    ::RichardsModel,
    ::AbstractWaterBC,
    boundary::ClimaLand.TopBoundary,
    Δz,
    Y,
    p,
    t,
)
    return ClimaCore.Fields.FieldVector(; :soil => (; :ϑ_l => zeros(axes(Δz))))
end

# BC type for the soil heat equation
"""
    AbstractHeatBC <: ClimaLand.AbstractBC

An abstract type for boundary conditions for the soil heat equation.
"""
abstract type AbstractHeatBC <: ClimaLand.AbstractBC end

"""
   TemperatureStateBC <: AbstractHeatBC

A simple concrete type of boundary condition, which enforces a
state boundary condition T = f(p,t) at either the top or bottom of the domain.
"""
struct TemperatureStateBC{F <: Function} <: AbstractHeatBC
    bc::F
end

"""
   HeatFluxBC <: AbstractHeatBC

A simple concrete type of boundary condition, which enforces a
normal flux boundary condition f(p,t) at either the top or bottom of the domain.
"""
struct HeatFluxBC{F <: Function} <: AbstractHeatBC
    bc::F
end


# Methods

"""
    boundary_flux(bc::HeatFluxBC,  _...)::ClimaCore.Fields.Field

A method of boundary fluxes which returns the desired flux.

We add a field of zeros in order to convert the bc (float) into
a field.
"""
function boundary_flux(
    bc::HeatFluxBC,
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
    boundary_flux(heat_bc::TemperatureStateBC,
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
function boundary_flux(
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
    boundary_flux(heat_bc::TemperatureStateBC,
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
function boundary_flux(
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


# BC type for the combined energy-hydrology equations
"""
    AbstractEnergyHydrologyBC <: ClimaLand.AbstractBC

An abstract type for boundary conditions for the combined energy
and hydrology equations.

In general these can consist of independent boundary conditions for
water (<: AbstractWaterBC) and for heat (<: AbstractHeatBC) independently,
, or a boundary condition type which sets both using the same parameterizations.
"""
abstract type AbstractEnergyHydrologyBC <: ClimaLand.AbstractBC end

"""
    WaterHeatBC{W <: AbstractWaterBC, H <: AbstractHeatBC} <:
       AbstractEnergyHydrologyBC

A general struct used to store the boundary conditions for Richards 
and the soil heat equations separately; useful when the boundary
conditions for each component are independent of each other.
"""
struct WaterHeatBC{W <: AbstractWaterBC, H <: AbstractHeatBC} <:
       AbstractEnergyHydrologyBC
    water::W
    heat::H
end
function WaterHeatBC(; water, heat)
    return WaterHeatBC{typeof(water), typeof(heat)}(water, heat)
end


"""
    AtmosDrivenFluxBC{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
        R <: AbstractRunoffModel
    } <: AbstractEnergyHydrologyBC

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
} <: AbstractEnergyHydrologyBC
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

# Special methods for EnergyHydrology BC
"""
    boundary_var_types(::Soil.EnergyHydrology{FT}, ::AbstractEnergyHydrologyBC, ::ClimaLand.AbstractBoundary) where {FT}

The list of domain names for additional variables added to the
EnergyHydrology model auxiliary state, which defaults to adding storage for the
 boundary flux field. 

Because we supply boundary conditions for water and heat, we found it convenient to
have these stored as a NamedTuple under the names `top_bc` and `bottom_bc`.
"""
function boundary_var_types(
    ::Soil.EnergyHydrology{FT},
    ::AbstractEnergyHydrologyBC,
    ::ClimaLand.AbstractBoundary,
) where {FT}
    (NamedTuple{(:water, :heat), Tuple{FT, FT}},)
end

"""
    soil_boundary_fluxes!(bc::WaterHeatBC, boundary::TopBoundary, model, Δz, Y, p, t)

updates the boundary fluxes for ϑ_l and ρe_int.
"""
function soil_boundary_fluxes!(
    bc::WaterHeatBC,
    boundary::AbstractBoundary,
    model,
    Δz,
    Y,
    p,
    t,
)
    params = model.parameters
    name = ClimaLand.bc_name(boundary)
    water_bc = getproperty(p.soil, name).water
    heat_bc = getproperty(p.soil, name).heat
    water_bc .= boundary_flux(bc.water, boundary, model, Δz, Y, p, t)
    heat_bc .= boundary_flux(bc.heat, boundary, model, Δz, Y, p, t)
end

"""
    boundary_vars(::AtmosDrivenFluxBC{<:AbstractAtmosphericDrivers,
                                    <:AbstractRadiativeDrivers,
                                    <:AbstractRunoffModel,
                                    }, ::ClimaLand.TopBoundary)

An extension of the `boundary_vars` method for AtmosDrivenFluxBC. This
adds the surface conditions (SHF, LHF, evaporation, and resistance) and the
net radiation to the auxiliary variables.

These variables are updated in place in `soil_boundary_fluxes!`.
"""
boundary_vars(
    bc::AtmosDrivenFluxBC{
        <:AbstractAtmosphericDrivers,
        <:AbstractRadiativeDrivers,
        <:AbstractRunoffModel,
    },
    ::ClimaLand.TopBoundary,
) = (:turbulent_fluxes, :ice_frac, :R_n, :top_bc, :infiltration, :sfc_scratch)

"""
    boundary_var_domain_names(::AtmosDrivenFluxBC{<:AbstractAtmosphericDrivers,
                                                  <:AbstractRadiativeDrivers,
                                                  <:AbstractRunoffModel,
                                                  },
                              ::ClimaLand.TopBoundary)

An extension of the `boundary_var_domain_names` method for AtmosDrivenFluxBC. This
specifies the part of the domain on which the additional variables should be
defined.
"""
boundary_var_domain_names(
    bc::AtmosDrivenFluxBC{
        <:AbstractAtmosphericDrivers,
        <:AbstractRadiativeDrivers,
        <:AbstractRunoffModel,
    },
    ::ClimaLand.TopBoundary,
) = (:surface, :surface, :surface, :surface, :surface, :surface)
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
        <:AbstractAtmosphericDrivers{FT},
        <:AbstractRadiativeDrivers{FT},
        <:AbstractRunoffModel,
    },
    ::ClimaLand.TopBoundary,
) where {FT} = (
    NamedTuple{(:lhf, :shf, :vapor_flux, :r_ae), Tuple{FT, FT, FT, FT}},
    FT,
    FT,
    NamedTuple{(:water, :heat), Tuple{FT, FT}},
    FT,
    FT,
)

"""
    soil_boundary_fluxes!(
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
function soil_boundary_fluxes!(
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
    update_runoff!(p, bc.runoff, Y, t, model)
    # We do not model the energy flux from infiltration. We multiply
    # the vapor flux by the ice fraction in order to get the contribution
    # from liquid water
    @. p.soil.top_bc.water =
        p.soil.infiltration +
        p.soil.turbulent_fluxes.vapor_flux * (1 - p.soil.ice_frac)
    @. p.soil.top_bc.heat =
        p.soil.R_n + p.soil.turbulent_fluxes.lhf + p.soil.turbulent_fluxes.shf
end

"""
    boundary_vars(::AtmosDrivenFluxBC{<:AbstractAtmosphericDrivers,
                                    <:AbstractRadiativeDrivers,
                                    <:Runoff.TOPMODELRunoff,
                                    }, ::ClimaLand.TopBoundary)

An extension of the `boundary_vars` method for AtmosDrivenFluxBC with 
TOPMODELRunoff. This
adds the surface conditions (SHF, LHF, evaporation, and resistance) and the
net radiation to the auxiliary variables.

These variables are updated in place in `soil_boundary_fluxes!`.
"""
boundary_vars(
    bc::AtmosDrivenFluxBC{
        <:AbstractAtmosphericDrivers,
        <:AbstractRadiativeDrivers,
        <:Runoff.TOPMODELRunoff,
    },
    ::ClimaLand.TopBoundary,
) = (
    :turbulent_fluxes,
    :ice_frac,
    :R_n,
    :top_bc,
    :h∇,
    :R_s,
    :R_ss,
    :infiltration,
    :sfc_scratch,
)

"""
    boundary_var_domain_names(::AtmosDrivenFluxBC{<:AbstractAtmosphericDrivers,
                                                  <:AbstractRadiativeDrivers,
                                                  <:Runoff.TOPMODELRunoff,
                                                  },
                              ::ClimaLand.TopBoundary)

An extension of the `boundary_var_domain_names` method for AtmosDrivenFluxBC
with TOPMODELRunoff. This
specifies the part of the domain on which the additional variables should be
defined.
"""
boundary_var_domain_names(
    bc::AtmosDrivenFluxBC{
        <:AbstractAtmosphericDrivers,
        <:AbstractRadiativeDrivers,
        <:Runoff.TOPMODELRunoff,
    },
    ::ClimaLand.TopBoundary,
) = (
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
)

"""
    boundary_var_types(
        ::AtmosDrivenFluxBC{
            <:PrescribedAtmosphere{FT},
            <:AbstractRadiativeDrivers{FT},
            <:Runoff.TOPMODELRunoff{FT},
        }, ::ClimaLand.TopBoundary,
    ) where {FT}

An extension of the `boundary_var_types` method for AtmosDrivenFluxBC
with TOPMODELRunoff. This
specifies the type of the additional variables.
"""
boundary_var_types(
    model::EnergyHydrology{FT},
    bc::AtmosDrivenFluxBC{
        <:AbstractAtmosphericDrivers{FT},
        <:AbstractRadiativeDrivers{FT},
        <:Runoff.TOPMODELRunoff{FT},
    },
    ::ClimaLand.TopBoundary,
) where {FT} = (
    NamedTuple{(:lhf, :shf, :vapor_flux, :r_ae), Tuple{FT, FT, FT, FT}},
    FT,
    FT,
    NamedTuple{(:water, :heat), Tuple{FT, FT}},
    FT,
    FT,
    FT,
    FT,
    FT,
)

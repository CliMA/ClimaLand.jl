import ClimaLSM: AbstractBC, boundary_flux
export TemperatureStateBC,
    MoistureStateBC, FreeDrainage, FluxBC, AtmosDrivenFluxBC

"""
    AbstractSoilBC <: ClimaLSM. AbstractBC

An abstract type for soil-specific types of boundary conditions, like free drainage.
"""
abstract type AbstractSoilBC <: ClimaLSM.AbstractBC end

"""
   MoistureStateBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
state boundary condition ϑ_l = f(p,t) at either the top or bottom of the domain.
"""
struct MoistureStateBC <: AbstractSoilBC
    bc::Function
end

"""
   TemperatureStateBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
state boundary condition T = f(p,t) at either the top or bottom of the domain.
"""
struct TemperatureStateBC <: AbstractSoilBC
    bc::Function
end

"""
   FluxBC <: AbstractSoilBC

A simple concrete type of boundary condition, which enforces a
normal flux boundary condition f(p,t) at either the top or bottom of the domain.
"""
struct FluxBC <: AbstractSoilBC
    bc::Function
end

"""
    AtmosDrivenFluxBC{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
    } <: AbstractSoilBC

A concrete type of soil boundary condition for use at the top
of the domain. This holds the conditions for the atmosphere
`AbstractAtmosphericDrivers` and for the radiation state 
`AbstractRadiativeDrivers`.

This choice indicates the Monin-Obukhov Surface Theory will
be used to compute the sensible and latent heat fluxes, as 
well as evaporation, and that the net radiation and precipitation
will also be computed. The net energy and water fluxes
are used as boundary conditions.
"""
struct AtmosDrivenFluxBC{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
} <: AbstractSoilBC
    atmos::A
    radiation::B
end

"""
    soil_boundary_fluxes(
        bc::AtmosDrivenFluxBC{
            <:PrescribedAtmosphere,
            <:PrescribedRadiativeFluxes,
        },
        boundary::ClimaLSM.TopBoundary,
        model::EnergyHydrology{FT},
        Y,
        Δz,
        p,
        t,
    ) where {FT}

Returns the net volumetric water flux (m/s) and net energy 
flux (W/m^2) for the soil `EnergyHydrology` model at the top 
of the soil domain.

This  method of `soil_boundary_fluxes` is for use with
a  `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
struct; for example, this is to be used when driving
the soil model with reanalysis data.

This function calls the `surface_fluxes` and `net_radiation`
functions, which use the soil surface conditions as well as 
the prescribed atmos and radiation conditions in order to
compute the surface fluxes.
"""
function soil_boundary_fluxes(
    bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:PrescribedRadiativeFluxes},
    boundary::ClimaLSM.TopBoundary,
    model::EnergyHydrology{FT},
    Y,
    Δz,
    p,
    t,
) where {FT}

    conditions = surface_fluxes(bc.atmos, model, Y, p, t)
    R_n = net_radiation(bc.radiation, model, Y, p, t)
    # We are ignoring sublimation for now, as well as g_soil
    net_water_flux = @. bc.atmos.liquid_precip(t) + conditions.vapor_flux
    net_energy_flux = @. R_n + conditions.lhf + conditions.shf
    return net_water_flux, net_energy_flux

end


"""
    FreeDrainage <: AbstractSoilBC
A concrete type of soil boundary condition, for use at 
the BottomBoundary only, where the flux is set to be
`F = -K∇h = -K`.
"""
struct FreeDrainage <: AbstractSoilBC end

"""
    ClimaLSM.boundary_flux(bc::FluxBC, _, Δz, _...)::ClimaCore.Fields.Field

A method of boundary fluxes which returns the desired flux.

We add a field of zeros in order to convert the bc (float) into
a field.
"""
function ClimaLSM.boundary_flux(
    bc::FluxBC,
    boundary::ClimaLSM.AbstractBoundary,
    Δz::ClimaCore.Fields.Field,
    p::ClimaCore.Fields.FieldVector,
    t,
    params,
)::ClimaCore.Fields.Field
    return bc.bc(p, t) .+ ClimaCore.Fields.zeros(axes(Δz))
end

"""
    ClimaLSM.boundary_flux(rre_bc::MoistureStateBC, ::ClimaLSM.TopBoundary, Δz, p, t, params)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on θ_l at the top of the
domain into a flux of liquid water.
"""
function ClimaLSM.boundary_flux(
    rre_bc::MoistureStateBC,
    ::ClimaLSM.TopBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field
    # Approximate K_bc ≈ K_c, ψ_bc ≈ ψ_c (center closest to the boundary)
    p_len = Spaces.nlevels(axes(p.soil.K))
    K_c = Fields.level(p.soil.K, p_len)
    ψ_c = Fields.level(p.soil.ψ, p_len)

    # Calculate pressure head using boundary condition
    @unpack vg_α, vg_n, vg_m, θ_r, ν, S_s = params
    θ_bc = rre_bc.bc(p, t)
    ψ_bc = @. pressure_head(vg_α, vg_n, vg_m, θ_r, θ_bc, ν, S_s)

    # Pass in (ψ_bc .+ Δz) as x_2 to account for contribution of gravity in RRE
    return ClimaLSM.diffusive_flux(K_c, ψ_bc .+ Δz, ψ_c, Δz)
end

"""
    ClimaLSM.boundary_flux(rre_bc::MoistureStateBC, ::ClimaLSM.BottomBoundary, Δz, p, t, params)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on θ_l at the bottom of the
domain into a flux of liquid water.
"""
function ClimaLSM.boundary_flux(
    rre_bc::MoistureStateBC,
    ::ClimaLSM.BottomBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field
    # Approximate K_bc ≈ K_c, ψ_bc ≈ ψ_c (center closest to the boundary)
    K_c = Fields.level(p.soil.K, 1)
    ψ_c = Fields.level(p.soil.ψ, 1)

    # Calculate pressure head using boundary condition
    @unpack vg_α, vg_n, vg_m, θ_r, ν, S_s = params
    θ_bc = rre_bc.bc(p, t)
    ψ_bc = @. pressure_head(vg_α, vg_n, vg_m, θ_r, θ_bc, ν, S_s)

    # At the bottom boundary, ψ_c is at larger z than ψ_bc
    #  so we swap their order in the derivative calc
    # Pass in (ψ_c .+ Δz) as x_2 to account for contribution of gravity in RRE
    return ClimaLSM.diffusive_flux(K_c, ψ_c .+ Δz, ψ_bc, Δz)
end


"""
    ClimaLSM.boundary_flux(heat_bc::TemperatureStateBC, ::ClimaLSM.TopBoundary, Δz, p, t)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on temperature at the top of the
domain into a flux of energy.
"""
function ClimaLSM.boundary_flux(
    heat_bc::TemperatureStateBC,
    ::ClimaLSM.TopBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    p_len = Spaces.nlevels(axes(p.soil.T))
    T_c = Fields.level(p.soil.T, p_len)
    κ_c = Fields.level(p.soil.κ, p_len)

    T_bc = heat_bc.bc(p, t)
    return ClimaLSM.diffusive_flux(κ_c, T_bc, T_c, Δz)
end

"""
    ClimaLSM.boundary_flux(heat_bc::TemperatureStateBC, ::ClimaLSM.BottomBoundary, Δz, p, t)::ClimaCore.Fields.Field

A method of boundary fluxes which converts a state boundary condition on temperature at the bottom of the
domain into a flux of energy.
"""
function ClimaLSM.boundary_flux(
    heat_bc::TemperatureStateBC,
    ::ClimaLSM.BottomBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field
    # Approximate κ_bc ≈ κ_c (center closest to the boundary)
    T_c = Fields.level(p.soil.T, 1)
    κ_c = Fields.level(p.soil.κ, 1)
    T_bc = heat_bc.bc(p, t)
    return ClimaLSM.diffusive_flux(κ_c, T_c, T_bc, Δz)
end


"""
    ClimaLSM.boundary_flux(bc::FreeDrainage{FT}, ::ClimaLSM.BottomBoundary, Δz, p, t)::ClimaCore.Fields.Field

A method of boundary fluxes which enforces free drainage at the bottom
of the domain.
"""
function ClimaLSM.boundary_flux(
    bc::FreeDrainage,
    boundary::ClimaLSM.BottomBoundary,
    Δz,
    p,
    t,
    params,
)::ClimaCore.Fields.Field
    K_c = Fields.level(p.soil.K, 1)
    return -1 .* K_c
end

"""
    soil_boundary_fluxes(bc::NamedTuple, boundary, model, Y, Δz, p, t)

Returns the boundary fluxes for ϑ_l and ρe_int, in that order.
"""
function soil_boundary_fluxes(bc::NamedTuple, boundary, model, Y, Δz, p, t)
    params = model.parameters
    return ClimaLSM.boundary_flux(bc.water, boundary, Δz, p, t, params),
    ClimaLSM.boundary_flux(bc.heat, boundary, Δz, p, t, params)
end

# TODO: generalize for atmos BC and for heat equation.
function ∂tendency_bc_∂ϑN(::MoistureStateBC, boundary::ClimaLSM.TopBoundary, Δz, Y, p, t, params)
    (; ν, vg_α, vg_n, vg_m, S_s, θ_r) = params
    K = ClimaLSM.Domains.top_center_to_surface(p.soil.K)
    ϑ_l = ClimaLSM.Domains.top_center_to_surface(Y.soil.ϑ_l)
    return @. -K/Δz *dψdθ(ϑ_l, ν, θ_r, vg_α, vg_n, vg_m, S_s)/ (2*Δz)
end

function ∂tendency_bc_∂ϑN(::AbstractSoilBC, boundary::ClimaLSM.TopBoundary, Δz, Y, p, t, params)
    return ClimaLSM.Domains.top_center_to_surface(zeros(axes(Y.soil.ϑ_l)))
end

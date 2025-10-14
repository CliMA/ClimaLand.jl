export LandHydrology, infiltration_capacity, infiltration_at_point

"""
    struct LandHydrology{
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        SW <: Pond.AbstractSurfaceWaterModel{FT},
    } <: AbstractLandModel{FT}

A concrete type of land model used for simulating systems with a
soil and surface water component.
$(DocStringExtensions.FIELDS)"""
struct LandHydrology{
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    SW <: Pond.AbstractSurfaceWaterModel{FT},
} <: AbstractLandModel{FT}
    "The soil model"
    soil::SM
    "The surface water model"
    surface_water::SW
end

"""
    LandHydrology{FT}(;
        land_args::NamedTuple = (;),
        soil_model_type::Type{SM},
        soil_args::NamedTuple = (;),
        surface_water_model_type::Type{SW},
        surface_water_args::NamedTuple = (;),
    ) where {
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        SW <: Pond.AbstractSurfaceWaterModel{FT},
    }
A constructor for the `LandHydrology` model, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `LandHydrology` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.

Additional arguments, like parameters and driving atmospheric data, are passed
in as `land_args`.
"""
function LandHydrology{FT}(;
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    surface_water_model_type::Type{SW},
    surface_water_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    SW <: Pond.AbstractSurfaceWaterModel{FT},
}

    (; precip) = land_args

    sources = ()
    surface_runoff = PrognosticRunoff{FT}(precip)
    boundary_conditions = (; top = RunoffBC(), bottom = Soil.FreeDrainage())

    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )
    domain = soil_args.domain
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_water = surface_water_model_type(;
        runoff = surface_runoff,
        domain = surface_domain,
    )
    args = (soil, surface_water)
    return LandHydrology{FT, typeof.(args)...}(args...)
end

lsm_aux_vars(m::LandHydrology) = (:soil_infiltration,)

lsm_aux_types(m::LandHydrology{FT}) where {FT} = (FT,)

lsm_aux_domain_names(m::LandHydrology) = (:surface,)

#=
If there is a pond present, flux BC should be -i_c

If there is no pond present (we won't model the standing surface water in full LSM)

 _______________________________P-E______________________________________________
                      |   P- E >0 (into soil)      |   P-E <0 (out of soil)
 |________________________________________________________________________________
 | Into Soil, i_c >0  |   min(P-E, i_c)            |   P-E -> min(P-E, i_c) works
I|________________________________________________________________________________
 |Out of soil, i_c <0 |   i_c; min(P-E, i_c) works |   Should be E? rare? unclear?
 |________________________________________________________________________________
=#

"""
    infiltration_at_point(η::FT, i_c::FT, P::FT)

Returns the infiltration given pond height η, infiltration capacity,
and precipitation.


This is defined such that positive means into soil.
"""
infiltration_at_point(η::FT, i_c::FT, P::FT) where {FT <: AbstractFloat} =
    η > eps(FT) ? i_c : min(i_c, P)

"""
    function infiltration_capacity(
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
    )

Function which computes the infiltration capacity of the soil based on
soil characteristics, moisture levels, and pond height.

Defined such that positive means into soil.
"""
function infiltration_capacity(Y::ClimaCore.Fields.FieldVector, p::NamedTuple)
    FT = eltype(Y.soil.ϑ_l)
    face_space = ClimaCore.Spaces.face_space(axes(Y.soil.ϑ_l))
    N = ClimaCore.Spaces.nlevels(face_space)
    space = axes(Y.surface_water.η)
    z = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z

    z_center = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(ClimaCore.Fields.level(z, N - 1)),
        space,
    )
    ψ_center = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(ClimaCore.Fields.level(p.soil.ψ, N - 1)),
        space,
    )
    K_center = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(ClimaCore.Fields.level(p.soil.K, N - 1)),
        space,
    )
    z_face = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(
            ClimaCore.Fields.level(
                ClimaCore.Fields.coordinate_field(face_space).z,
                N - ClimaCore.Utilities.half,
            ),
        ),
        space,
    )
    ψ_face = max.(FT(0), Y.surface_water.η)
    return @. (
        K_center * (ψ_face + z_face - (ψ_center + z_center)) /
        (z_face - z_center)
    )
end

"""
    make_update_explicit_boundary_fluxes(
        land::LandHydrology{FT, SM, SW},
    ) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}

A method which makes a function; the returned function
updates the auxiliary variable `p.soil_infiltration`, which
is needed for both the boundary condition for the soil model and the source
term (runoff) for the surface water model.

This function is called each ode function evaluation.
"""
function make_update_explicit_boundary_fluxes(
    land::LandHydrology{FT, SM, SW},
) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}
    update_soil_bf! = make_update_explicit_boundary_fluxes(land.soil)
    update_pond_bf! = make_update_explicit_boundary_fluxes(land.surface_water)
    function update_boundary_fluxes!(p, Y, t)
        update_soil_bf!(p, Y, t)
        update_pond_bf!(p, Y, t)
    end
    return update_boundary_fluxes!
end


"""
    make_update_implicit_boundary_fluxes(
        land::LandHydrology{FT, SM, SW},
    ) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}

A method which makes a function; the returned function
updates the auxiliary variable `p.soil_infiltration`, which
is needed for both the boundary condition for the soil model and the source
term (runoff) for the surface water model.

This function is called each ode function evaluation.
"""
function make_update_implicit_boundary_fluxes(
    land::LandHydrology{FT, SM, SW},
) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}
    update_soil_bf! = make_update_implicit_boundary_fluxes(land.soil)
    update_pond_bf! = make_update_implicit_boundary_fluxes(land.surface_water)
    function update_boundary_fluxes!(p, Y, t)
        i_c = infiltration_capacity(Y, p)
        @. p.soil_infiltration =
            -infiltration_at_point(
                Y.surface_water.η,
                i_c,
                -FT(land.surface_water.runoff.precip(t)),
            )
        update_soil_bf!(p, Y, t)
        update_pond_bf!(p, Y, t)
    end
    return update_boundary_fluxes!
end

"""
    PrognosticRunoff <: Pond.AbstractSurfaceRunoff

Concrete type of `Pond.AbstractSurfaceRunoff` for use in LSM models,
where precipitation is passed in, but infiltration is computed
prognostically.

This is paired with `Soil.RunoffBC`: both are used at the same time,
ensuring the infiltration used for the boundary condition of soil
is also used to compute the runoff for the surface water.

"""
struct PrognosticRunoff{FT, F <: Function} <: Pond.AbstractSurfaceRunoff
    precip::F
end

function PrognosticRunoff{FT}(precip) where {FT <: AbstractFloat}
    return PrognosticRunoff{FT, typeof(precip)}(precip)
end

"""
    function Pond.surface_runoff(
        runoff::PrognosticRunoff,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

Extension of the `Pond.surface_runoff` function, which computes
 the surface runoff, for use in an LSM when the runoff is determined
prognostically.
"""
function Pond.surface_runoff(
    runoff::PrognosticRunoff,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    FT = FTfromY(Y)
    return @. FT(p.soil_infiltration - runoff.precip(t))
end

"""
    RunoffBC <: Soil.AbstractSoilBC

Concrete type of `Soil.AbstractSoilBC` for use in LSM models,
where precipitation is passed in, but infiltration is computed
prognostically. This infiltration is then used to set an upper
boundary condition for the soil.


This is paired with `Pond.PrognosticRunoff`: both are used at the same
time,
ensuring that the infiltration used for the boundary condition of soil
is also used to compute the runoff for the surface water.
"""
struct RunoffBC <: Soil.AbstractWaterBC end

"""
    function ClimaLand.boundary_flux!(bc_field,
        bc::RunoffBC,
        ::TopBoundary,
        model::Soil.RichardsModel,
        Δz::FT,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
        params,
    )

Extension of the `ClimaLand.boundary_flux` function, which returns the water volume
boundary flux for the soil.
At the top boundary, return the soil infiltration (computed each step and
stored in `p.soil_infiltration`).
"""
function ClimaLand.boundary_flux!(
    bc_field,
    bc::RunoffBC,
    ::TopBoundary,
    model::Soil.RichardsModel,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    bc_field .= p.soil_infiltration
end

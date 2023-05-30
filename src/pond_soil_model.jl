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

    @unpack precip = land_args

    sources = ()
    surface_runoff = PrognosticRunoff{FT}(precip)
    boundary_conditions =
        (; top = (water = RunoffBC(),), bottom = (water = Soil.FreeDrainage(),))

    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )
    surface_water = surface_water_model_type(;
        runoff = surface_runoff,
        surface_water_args...,
    )
    args = (soil, surface_water)
    return LandHydrology{FT, typeof.(args)...}(args...)
end

interaction_vars(m::LandHydrology) = (:soil_infiltration,)

interaction_types(m::LandHydrology{FT}) where {FT} = (FT,)

interaction_domains(m::LandHydrology) = (:surface,)

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
infiltration_at_point(η::FT, i_c::FT, P::FT) where {FT} =
    η > eps(FT) ? i_c : min(i_c, P)

"""
    function infiltration_capacity(
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
    )

Function which computes the infiltration capacity of the soil based on
soil characteristics, moisture levels, and pond height.

Defined such that positive means into soil.
"""
function infiltration_capacity(
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
)
    face_space = ClimaLSM.Domains.obtain_face_space(axes(Y.soil.ϑ_l))
    N = ClimaCore.Spaces.nlevels(face_space)
    space = axes(Y.surface_water.η)
    z = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z

    z_center = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(ClimaCore.Fields.level(z, N - 1)),
        space,
    )
    ψ_center = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(ClimaCore.Fields.level(p.ψ, N - 1)),
        space,
    )
    K_center = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(ClimaCore.Fields.level(p.K, N - 1)),
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
    ψ_face = max.(0.0, Y.surface_water.η)
    return @. (
        K_center * (ψ_face + z_face - (ψ_center + z_center)) /
        (z_face - z_center)
    )
end

"""
    function make_interactions_update_aux(
        land::LandHydrology{FT, SM, SW},
    ) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}

A method which makes a function; the returned function 
updates the auxiliary variable `p.soil_infiltration`, which
is needed for both the boundary condition for the soil model and the source
term (runoff) for the surface water model.

This function is called each ode function evaluation.
"""
function make_interactions_update_aux(
    land::LandHydrology{FT, SM, SW},
) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}
    function update_aux!(p, Y, t)
        i_c = infiltration_capacity(Y, p.soil)
        @. p.soil_infiltration =
            -infiltration_at_point(
                Y.surface_water.η,
                i_c,
                -land.surface_water.runoff.precip(t),
            )
    end
    return update_aux!
end

"""
    PrognosticRunoff{FT} <: Pond.AbstractSurfaceRunoff{FT}

Concrete type of `Pond.AbstractSurfaceRunoff` for use in LSM models,
where precipitation is passed in, but infiltration is computed
prognostically.

This is paired with `Soil.RunoffBC`: both are used at the same time,
ensuring the infiltration used for the boundary condition of soil
is also used to compute the runoff for the surface water.

"""
struct PrognosticRunoff{FT} <: Pond.AbstractSurfaceRunoff{FT}
    precip::Function
end

"""
    function Pond.surface_runoff(
        runoff::PrognosticRunoff{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    ) where {FT}

Extension of the `Pond.surface_runoff` function, which computes the surface runoff, for use in an LSM when the runoff is determined
prognostically.
"""
function Pond.surface_runoff(
    runoff::PrognosticRunoff{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}
    return @. -(runoff.precip(t) - p.soil_infiltration)
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
struct RunoffBC <: Soil.AbstractSoilBC end

"""
    function ClimaLSM.boundary_flux(
        bc::RunoffBC,
        ::TopBoundary,
        Δz::FT,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
        params,
    )::ClimaCore.Fields.Field where {FT}

Extension of the `ClimaLSM.boundary_flux` function, which returns the water volume
boundary flux for the soil.
At the top boundary, return the soil infiltration (computed each step and
stored in `p.soil_infiltration`).
"""
function ClimaLSM.boundary_flux(
    bc::RunoffBC,
    ::TopBoundary,
    Δz::ClimaCore.Fields.Field,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
    params,
)::ClimaCore.Fields.Field where {FT}
    return p.soil_infiltration
end

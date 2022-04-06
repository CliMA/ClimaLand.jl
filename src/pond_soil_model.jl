export LandHydrology, infiltration_capacity, infiltration_at_point


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
    boundary_conditions = RunoffBC{FT}()

    soil = soil_model_type(;
        boundary_conditions = RunoffBC{FT}(),
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

land_components(land::LandHydrology) = (:soil, :surface_water)

function initialize_interactions(
    land::LandHydrology{FT, SM, SW},
) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}

    surface_coords = coordinates(land.surface_water)#In non-column case, would this be level(soil_coords)?
    return (soil_infiltration = similar(surface_coords),)
end

# defined such that positive means into soil
# If there is a pond present, flux BC should be -i_c
# If there is no pond present (we won't model the standing surface water in full LSM)
#=
 _______________________________P-E______________________________________________
                      |   P- E >0 (into soil)      |   P-E <0 (out of soil)
 |________________________________________________________________________________
 | Into Soil, i_c >0  |   min(P-E, i_c)            |   P-E -> min(P-E, i_c) works
I|________________________________________________________________________________
 |Out of soil, i_c <0 |   i_c; min(P-E, i_c) works |   Should be E? rare? unclear?
 |________________________________________________________________________________
=#

infiltration_at_point(η::FT, i_c::FT, P::FT) where {FT} =
    η > eps(FT) ? i_c : min(i_c, P)

# defined such that positive means into soil
function infiltration_capacity(
    model::Soil.AbstractSoilModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
) where {FT}
    N = length(axes(p.ψ).center_local_geometry)
    Δz = last(axes(p.ψ).face_local_geometry.WJ)
    z_surf = Δz + last(parent(model.coordinates))
    ψ, K =
        ClimaCore.Operators.getidx.(
            [p.ψ, p.K],
            Ref(ClimaCore.Operators.Interior()),
            N,
        )
    return @. (K * (max(0.0, Y.surface_water.η) - ψ + Δz) / Δz)
end


function make_interactions_update_aux(
    land::LandHydrology{FT, SM, SW},
) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}
    function update_aux!(p, Y, t)
        i_c = infiltration_capacity(land.soil, Y, p.soil)
        @. p.soil_infiltration =
            -infiltration_at_point(
                Y.surface_water.η,
                i_c,
                -land.surface_water.runoff.precip(t),
            )
    end
    return update_aux!
end

struct PrognosticRunoff{FT} <: Pond.AbstractSurfaceRunoff{FT}
    precip::Function
end


function Pond.surface_runoff(
    runoff::PrognosticRunoff{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}
    return @. -(runoff.precip(t) - p.soil_infiltration)
end


struct RunoffBC{FT} <: Soil.AbstractSoilBoundaryConditions{FT} end

function Soil.boundary_fluxes(
    bc::RunoffBC{FT},
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}
    return p.soil_infiltration[1], FT(0.0)
end

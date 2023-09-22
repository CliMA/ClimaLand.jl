module Runoff
using ClimaLSM
import ClimaLSM: source!
using ..ClimaLSM.Soil: AbstractSoilSource, AbstractSoilModel
export soil_surface_infiltration,
    TOPMODELRunoff,
    AbstractRunoffModel,
    TOPMODELSubsurfaceRunoff,
    subsurface_runoff_source,
    append_source

"""
    AbstractRunoffModel

The abstract type for soil runoff models to be used with
the following boundary condition types:
- `ClimaLSM.Soil.AtmosDrivenFluxBC`
- `ClimaLSM.Soil.RichardsAtmosDrivenFluxBC`,

and for these functions:
-`ClimaLSM.Soil.soil_surface_infiltration`
- `ClimaLSM.Soil.subsurface_runoff_source`
- `ClimaLSM.source!`.

Please see the documentation for these for more details.
The model should specify the subsurface runoff sink term as well
as the surface runoff implementation.
"""
abstract type AbstractRunoffModel end

"""
    NoRunoff <: AbstractRunoffModel

A concrete type of soil runoff model; the 
default choice which does not include the 
effects of runoff.
"""
struct NoRunoff <: AbstractRunoffModel end

"""
    soil_surface_infiltration(::NoRunoff, net_water_flux, _...)

A function which computes the infiltration into the soil
 for the default of `NoRunoff`.

If `net_water_flux = P+E`, where `P` is the precipitation and `E`
is the evaporation (both negative if towards the soil), 
this returns `P+E` as the water boundary flux for the soil.
"""
soil_surface_infiltration(::NoRunoff, net_water_flux, _...) = net_water_flux

"""
    subsurface_runoff_source(runoff::AbstractRunoffModel)::Union{Nothing, AbstractSoilSource} 

A function which returns the soil source for the runoff model 
`runoff`; the default returns nothing in which case no source is added.

"""
subsurface_runoff_source(
    runoff::AbstractRunoffModel,
)::Union{Nothing, AbstractSoilSource} = nothing

"""
    append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple

Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple = (srcs..., src)

"""
    append_source(src::Nothing , srcs::Tuple)::Tuple

Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::Nothing, srcs::Tuple)::Tuple = srcs

struct TOPMODELSubsurfaceRunoff{FT} <: AbstractSoilSource{FT} end

struct TOPMODELRunoff{FT <: AbstractFloat} <: AbstractRunoffModel
    f_over::FT
    regrid_dirpath::String
    raw_datapath::String
    subsurface_source::TOPMODELSubsurfaceRunoff{FT}
end

function set_initial_parameter_field!(
    parameterization::TOPMODELRunoff{FT},
    p,
    surface_space,
) where {FT}
    #    comms =
    p.soil.fmax .= regrid_netcdf_to_field(
        FT,
        parameterization.regrid_dirpath,
        comms,
        parameterization.raw_datapath,
        :fmax,
        surface_space,
    )

    p.soil.landsea_mask .= regrid_netcdf_to_field(
        FT,
        parameterization.regrid_dirpath,
        comms,
        parameterization.raw_datapath,
        :landsea_mask,
        surface_space,
    )
end

"""
TOPMODEL infiltration
"""
function soil_surface_infiltration(
    runoffmodel::TOPMODELRunoff,
    net_water_flux,
    Y,
    p,
    soil_parameters,
)

    (; hydrology_cm, K_sat, ν, θ_r) = soil_parameters
    # we need this at the surface only
    infiltration_capacity = @. -p.soil.K / hydraulic_conductivity(
        hydrology_cm,
        K_sat,
        effective_saturation(ν, Y.soil.ϑ_l, θ_r),
    )
    # net_water_flux is negative if towards the soil; take smaller in magnitude -> max
    # net_water_flux is positive if away from soil -> use as BC.
    return (1 - p.soil.fsat) * max(infiltration_capacity, net_water_flux)
end

"""
    TOPMODEL sink term for baseflow - standin
"""
function ClimaLSM.source!(
    dY,
    src::TOPMODELSubsurfaceRunoff,
    Y,
    p,
    model::AbstractSoilModel,
)
    #=
    h_∇ = sum(heaviside.(Y.soil.ϑ_l .- ν))
    @. dY.soil.ϑ_l -=
        -p.soil.K / hydraulic_conductivity(
            hydrology_cm,
            K_sat,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        ) / model.boundary_conditions.top.runoff.f * exp(-p.soil.ϕ_sat) / h_∇ *
        heaviside.(Y.soil.ϑ_l .- ν)
    =#

end

subsurface_runoff_source(runoff::TOPMODELRunoff) = runoff.subsurface_source
end

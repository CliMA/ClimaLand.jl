module DALEC

using DocStringExtensions
using ClimaCore
using ClimaLSM

export DALEC811Parameters, DALECModel, ACM
export dalec_811_parmin, dalec_811_parmax, dalec_811_parnames,  check_dalec_811_parameter_bounds
export DALEC811AtmosphericDrivers, DALEC811FireDrivers, DALEC811RadiativeDrivers, DALEC811TemporalDrivers
export load_initial_condition!

include("./auxi/ACM.jl")
include("./auxi/seasonality.jl")


"""
    AbstractDALECParameters{FT <: AbstractFloat}

An abstract type of DALEC Parameters.
"""
abstract type AbstractDALECParameters{FT <:AbstractFloat} end


"""
    DALEC811AtmosphericDrivers{T <: AbstractFloat, TMIN, TMAX, C, V, P} <:AbstractAtmosphericDrivers{T}    

Atmospheric drivers for DALEC811 model.
$(DocStringExtensions.FIELDS)
"""
struct DALEC811AtmosphericDrivers{T <: AbstractFloat, TMIN, TMAX, C, V, P} <:AbstractAtmosphericDrivers{T}  
    T_MIN::TMIN
    T_MAX::TMAX
    ATMOSPHERIC_CO2::C
    VPD::V
    PRECIP::P
    LATITUDE::T
    MEAN_TEMP::T
    MEAN_PRECIP::T
    FT::Type{T}
end


"""
    DALEC811RadiativeDrivers{T <: AbstractFloat, SRD} <:AbstractRadiativeDrivers{T}

Radiative drivers for the DALEC811 model.
$(DocStringExtensions.FIELDS)
"""
struct DALEC811RadiativeDrivers{T <: AbstractFloat, SRD} <:AbstractRadiativeDrivers{T}  
    SSRD::SRD
    FT::Type{T}
end

"""
    DALEC811FireDrivers{T <: AbstractFloat, BURN}

Fire drivers for the DALEC811 model.
$(DocStringExtensions.FIELDS)
"""
struct DALEC811FireDrivers{T <: AbstractFloat, BURN} 
    BURNED_AREA::BURN
    FT::Type{T}
end

"""
    DALEC811TemporalDrivers{T <: AbstractFloat, TI, D, I}

Temporal drivers for the DALEC811 model.
$(DocStringExtensions.FIELDS)
"""
struct DALEC811TemporalDrivers{T <: AbstractFloat, TI, D, I}
    TIME::TI
    DOY::D
    DELTA_T::T
    NODAYS::I
    FT::Type{T}
end

"""
    DALEC811Parameters{FT <: AbstractFloat} <: AbstractDALECParameters{FT} 

Parameters for the DALEC811 model.
$(DocStringExtensions.FIELDS)
"""
struct DALEC811Parameters{FT <: AbstractFloat} <: AbstractDALECParameters{FT} 
     # Trainable DALEC Parameters
     decomposition_rate::FT
     f_gpp::FT
     f_fol::FT
     f_root::FT
     leaf_lifespan::FT
     tor_wood::FT
     tor_root::FT
     tor_litter::FT
     tor_som::FT
     Q10::FT
     canopy_efficiency::FT
     Bday::FT
     f_lab::FT
     clab_release_period::FT
     Fday::FT
     leaf_fall_period::FT
     LMCA::FT
     Clab::FT
     Cfol::FT
     Croot::FT
     Cwood::FT
     Clitter::FT
     Csom::FT
     IWUE::FT
     runoff_focal_point::FT
     wilting_point::FT
     initial_water::FT
     foliar_cf::FT
     ligneous_cf::FT
     dom_cf::FT
     resilience::FT
     lab_lifespan::FT
     moisture_factor::FT 
end

"""
    DALECModel{FTT, A, R, F, T, P, PA, D} <: AbstractModel{FTT}

A generic DALEC Model.
"""
struct DALECModel{FTT, A, R, F, T, P, PA, D} <: AbstractModel{FTT}
    # drivers
    atmos::A
    rad::R
    fire::F
    temporal::T
    # submodels
    photosynthesis::P
    # parameters
    parameters::PA
    # domain
    domain::D
    # DALEC ID
    id::Symbol
    # Float type
    FT::Type{FTT}
end

ClimaLSM.name(model::DALECModel) = model.id
ClimaLSM.domain(::DALECModel) = :surface


"""
    ClimaLSM.prognostic_vars(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
     R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}

Return prognostic variables for the DALEC811 model.
"""
function ClimaLSM.prognostic_vars(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
     R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}
    return (:LAI, :GPP, :ET, :temperate, :respiration_auto, :leaf_production, :labile_production, 
    :root_production, :wood_production, :lff, :lrf, :labile_release, :leaf_litter, :wood_litter,
      :root_litter, :respiration_hetero_litter, :respiration_hetero_som, :litter_to_som, :runoff,
       :labile_fire_combust, :foliar_fire_combust, :root_fire_combust, :wood_fire_combust,
        :litter_fire_combust, :som_fire_combust, :labile_fire_transfer, :foliar_fire_transfer,
         :root_fire_transfer, :wood_fire_transfer, :litter_fire_transfer, :total_fire_combust,
          :nee, :next_labile_pool, :next_foliar_pool, :next_root_pool, :next_wood_pool, :next_litter_pool,
          :next_som_pool, :next_water_pool)
end

"""
    ClimaLSM.prognostic_types(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
    R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}

Return prognostic types for the DALEC811 model.
"""
function ClimaLSM.prognostic_types(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
    R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}
    return (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT,
    FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT);
end


"""
    ClimaLSM.auxiliary_vars(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
    R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}

Return auxiliary variables for the DALEC811 model.
"""
function ClimaLSM.auxiliary_vars(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
    R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}
    return (:LAI, :GPP, :ET, :temperate, :respiration_auto, :leaf_production, :labile_production, 
    :root_production, :wood_production, :lff, :lrf, :labile_release, :leaf_litter, :wood_litter,
      :root_litter, :respiration_hetero_litter, :respiration_hetero_som, :litter_to_som, :runoff,
       :labile_fire_combust, :foliar_fire_combust, :root_fire_combust, :wood_fire_combust,
        :litter_fire_combust, :som_fire_combust, :labile_fire_transfer, :foliar_fire_transfer,
         :root_fire_transfer, :wood_fire_transfer, :litter_fire_transfer, :total_fire_combust,
          :nee, :next_labile_pool, :next_foliar_pool, :next_root_pool, :next_wood_pool, :next_litter_pool,
          :next_som_pool, :next_water_pool)
end

"""
    ClimaLSM.auxiliary_types(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
    R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}

Return auxiliary types for the DALEC811 model.
"""
function ClimaLSM.auxiliary_types(model::DALECModel{FT, A, R, F, T, P, PA, D}) where{FT, A<:DALEC811AtmosphericDrivers,
    R<:DALEC811RadiativeDrivers, F<:DALEC811FireDrivers, T<:DALEC811TemporalDrivers, P<:ACM, PA<:DALEC811Parameters, D}
    return (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT,
    FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT);
end


"""
    DALECModel(atmos::DALEC811AtmosphericDrivers, 
    rad::DALEC811RadiativeDrivers,
    fire::DALEC811FireDrivers,
    temporal::DALEC811TemporalDrivers,
    photosynthesis::ACM, 
    parameters::DALEC811Parameters, 
    domains::Union{
        ClimaLSM.Domains.Point,
        ClimaLSM.Domains.Plane,
        ClimaLSM.Domains.SphericalSurface,
    }) where {FT}

Outer constructor for the DALEC811 model. 
$(DocStringExtensions.FIELDS)
"""
function DALECModel(atmos::DALEC811AtmosphericDrivers, 
    rad::DALEC811RadiativeDrivers,
    fire::DALEC811FireDrivers,
    temporal::DALEC811TemporalDrivers,
    photosynthesis::ACM, 
    parameters::DALEC811Parameters{FT}, 
    domains::Union{
        ClimaLSM.Domains.Point,
        ClimaLSM.Domains.Plane,
        ClimaLSM.Domains.SphericalSurface,
    }) where {FT}
    return DALECModel(atmos, rad, fire, temporal, photosynthesis, parameters, domains, :dalec811, FT)
end


"""
    ClimaLSM.make_rhs(model::DALECModel{FT, A, R, F, T, P, PA, D}) where {FT, A, R, F, T, P, PA, D}

The method extends ClimaLSM.make_rhs. The ode function is discrete and is expected to run on weekly to monthly temporal resolution.
To intrage the model, the function should be cast in a DiscreteProblem and intetrated using FunctionMap{true} of the OrdinaryDiffEq module.
"""
function ClimaLSM.make_rhs(model::DALECModel{FT, A, R, F, T, P, PA, D}) where {FT, A, R, F, T, P, PA, D}
    function rhs!(dY, Y, p, t)
        
        # normalize foliar pool by LCMA (leaf carbon mass per area) to obtain LAI. 
        @. p.dalec811.LAI = Y.dalec811.next_foliar_pool / model.parameters.LMCA
        
        raw_gpp = model.photosynthesis(lat = model.atmos.LATITUDE, doy = model.temporal.DOY(t), t_max = model.atmos.T_MAX(t),
            t_min = model.atmos.T_MIN(t), lai = p.dalec811.LAI, rad = model.rad.SSRD(t),
            ca = model.atmos.ATMOSPHERIC_CO2(t), ce = model.parameters.canopy_efficiency)
        
        # compute gpp FLUXES[0]
        @.p.dalec811.GPP = raw_gpp .* min.(Y.dalec811.next_water_pool ./ model.parameters.wilting_point, 1)
        
        # compute ET_flux FLUXES[28]
        @. p.dalec811.ET = p.dalec811.GPP * model.atmos.VPD(t) / model.parameters.IWUE 
        
        # compute temperate FLUXES[1]
        @. p.dalec811.temperate = exp(model.parameters.Q10 * (FT(0.5) * (model.atmos.T_MAX(t) + model.atmos.T_MIN(t)) - model.atmos.MEAN_TEMP)) * ((model.atmos.PRECIP(t) / model.atmos.MEAN_PRECIP - FT(1)) * model.parameters.moisture_factor + FT(1))
        
        # compute autotrophic respiration FLUXES[2]
        @. p.dalec811.respiration_auto = model.parameters.f_gpp * p.dalec811.GPP
        
        # compute leaf production FLUXES[3]
        @. p.dalec811.leaf_production = (p.dalec811.GPP - p.dalec811.respiration_auto) * model.parameters.f_fol
        
        # compute labile production FLUXES[4]
        @. p.dalec811.labile_production = (p.dalec811.GPP - p.dalec811.respiration_auto - p.dalec811.leaf_production) * model.parameters.f_lab
        
        # compute root production FLUXES[5]
        @. p.dalec811.root_production = (p.dalec811.GPP - p.dalec811.respiration_auto - p.dalec811.leaf_production - p.dalec811.labile_production) * model.parameters.f_root
    
        # compute wood production FLUXES[6]
        @. p.dalec811.wood_production = p.dalec811.GPP -  p.dalec811.respiration_auto - p.dalec811.leaf_production - p.dalec811.labile_production - p.dalec811.root_production
        
        # compute leaf fall factor FLUXES[8]
        @. p.dalec811.lff = leaf_fall_factor(model.temporal.TIME(t), model.parameters.leaf_lifespan, model.parameters.leaf_fall_period, model.parameters.Fday, FT)
        
        # compute labile release factor FLUXES[15]
        @. p.dalec811.lrf = lab_release_factor(model.temporal.TIME(t), model.parameters.lab_lifespan, model.parameters.clab_release_period, model.parameters.Bday, FT)
        
        # compute labile release FLUXES[7]
        @. p.dalec811.labile_release = Y.dalec811.next_labile_pool * (FT(1) - (FT(1) - p.dalec811.lrf) ^ model.temporal.DELTA_T) / model.temporal.DELTA_T
        
        # compute leaf litter production FLUXES[9]
        @. p.dalec811.leaf_litter = Y.dalec811.next_foliar_pool * (FT(1) - (FT(1) - p.dalec811.lff) ^ model.temporal.DELTA_T) / model.temporal.DELTA_T
        
        # compute wood litter production FLUXES[10]
        @. p.dalec811.wood_litter = Y.dalec811.next_wood_pool * (FT(1) - (FT(1) - model.parameters.tor_wood)^ model.temporal.DELTA_T) / model.temporal.DELTA_T

        # compute root litter production FLUXES[11]
        @. p.dalec811.root_litter = Y.dalec811.next_root_pool * (FT(1) - (FT(1) - model.parameters.tor_root) ^ model.temporal.DELTA_T) / model.temporal.DELTA_T
        
        # compute respiration heterotrophic litter FLUXES[12]
        @. p.dalec811.respiration_hetero_litter = Y.dalec811.next_litter_pool * (FT(1) - (FT(1) - p.dalec811.temperate * model.parameters.tor_litter) ^ model.temporal.DELTA_T) / model.temporal.DELTA_T
        
        # compute respiration heterotrophic SOM FLUXES[13]
        @. p.dalec811.respiration_hetero_som = Y.dalec811.next_som_pool * (FT(1) - (FT(1) - p.dalec811.temperate * model.parameters.tor_som) ^ model.temporal.DELTA_T) / model.temporal.DELTA_T

        # compute litter to som flux
        @. p.dalec811.litter_to_som = Y.dalec811.next_litter_pool * (FT(1) - (FT(1) - p.dalec811.temperate * model.parameters.decomposition_rate) ^ model.temporal.DELTA_T) / model.temporal.DELTA_T
        
        # compute runoff flux
        @. p.dalec811.runoff = Y.dalec811.next_water_pool ^ FT(2) / model.parameters.runoff_focal_point / model.temporal.DELTA_T
                
        less_than_half_sel = map(Y.dalec811.next_water_pool) do v
            return Float64(v <= (model.parameters.runoff_focal_point ./ FT(2)))
        end
        
        greater_than_half_sel = map(Y.dalec811.next_water_pool) do v
            return Float64(v > (model.parameters.runoff_focal_point ./ FT(2)))
        end
        
        less_than_half_runoff = p.dalec811.runoff .* less_than_half_sel
        
        greater_than_half_runoff = ((Y.dalec811.next_water_pool .- model.parameters.runoff_focal_point ./ FT(4)) ./ model.temporal.DELTA_T) .* greater_than_half_sel

        @. p.dalec811.runoff = less_than_half_runoff + greater_than_half_runoff
        
        # all pools before including fire.
        @. p.dalec811.next_labile_pool = Y.dalec811.next_labile_pool + (p.dalec811.labile_production - p.dalec811.labile_release) * model.temporal.DELTA_T
        @. p.dalec811.next_foliar_pool = Y.dalec811.next_foliar_pool + (p.dalec811.leaf_production - p.dalec811.leaf_litter + p.dalec811.labile_release) * model.temporal.DELTA_T
        @. p.dalec811.next_root_pool = Y.dalec811.next_root_pool + (p.dalec811.root_production - p.dalec811.root_litter) * model.temporal.DELTA_T
        @. p.dalec811.next_wood_pool = Y.dalec811.next_wood_pool + (p.dalec811.wood_production - p.dalec811.wood_litter) * model.temporal.DELTA_T
        @. p.dalec811.next_litter_pool = Y.dalec811.next_litter_pool + (p.dalec811.leaf_litter + p.dalec811.root_litter - p.dalec811.respiration_hetero_litter - p.dalec811.litter_to_som) * model.temporal.DELTA_T
        @. p.dalec811.next_som_pool = Y.dalec811.next_som_pool + (p.dalec811.litter_to_som - p.dalec811.respiration_hetero_som + p.dalec811.wood_litter) * model.temporal.DELTA_T
        @. p.dalec811.next_water_pool = Y.dalec811.next_water_pool - p.dalec811.runoff * model.temporal.DELTA_T + model.atmos.PRECIP(t) * model.temporal.DELTA_T - p.dalec811.ET * model.temporal.DELTA_T
    
        
        @. p.dalec811.labile_fire_combust = p.dalec811.next_labile_pool * model.fire.BURNED_AREA(t) * model.parameters.ligneous_cf / model.temporal.DELTA_T
        @. p.dalec811.foliar_fire_combust = p.dalec811.next_foliar_pool * model.fire.BURNED_AREA(t) * model.parameters.foliar_cf / model.temporal.DELTA_T
        @. p.dalec811.root_fire_combust = p.dalec811.next_root_pool * model.fire.BURNED_AREA(t) * model.parameters.ligneous_cf / model.temporal.DELTA_T
        @. p.dalec811.wood_fire_combust = p.dalec811.next_wood_pool * model.fire.BURNED_AREA(t) * model.parameters.ligneous_cf / model.temporal.DELTA_T
        @. p.dalec811.litter_fire_combust = p.dalec811.next_litter_pool * model.fire.BURNED_AREA(t) * (model.parameters.ligneous_cf + model.parameters.foliar_cf) * FT(0.5) / model.temporal.DELTA_T
        @. p.dalec811.som_fire_combust = p.dalec811.next_som_pool * model.fire.BURNED_AREA(t) * model.parameters.dom_cf / model.temporal.DELTA_T
        
        @. p.dalec811.labile_fire_transfer = p.dalec811.next_labile_pool * model.fire.BURNED_AREA(t) *(FT(1) - model.parameters.ligneous_cf) * (FT(1) - model.parameters.resilience) / model.temporal.DELTA_T
        @. p.dalec811.foliar_fire_transfer = p.dalec811.next_foliar_pool * model.fire.BURNED_AREA(t) *(FT(1)- model.parameters.foliar_cf) * (FT(1)- model.parameters.resilience) / model.temporal.DELTA_T
        @. p.dalec811.root_fire_transfer = p.dalec811.next_root_pool * model.fire.BURNED_AREA(t) *(FT(1) - model.parameters.ligneous_cf) * (FT(1) - model.parameters.resilience) / model.temporal.DELTA_T
        @. p.dalec811.wood_fire_transfer = p.dalec811.next_wood_pool * model.fire.BURNED_AREA(t) * (FT(1) - model.parameters.ligneous_cf)*(FT(1) - model.parameters.resilience) / model.temporal.DELTA_T
        @. p.dalec811.litter_fire_transfer = p.dalec811.next_litter_pool * model.fire.BURNED_AREA(t) * (FT(1) -(model.parameters.ligneous_cf + model.parameters.foliar_cf) * FT(0.5)) * (FT(1) - model.parameters.resilience) / model.temporal.DELTA_T
        
        # include fire combust and fire transfer fluxes to the carbon pools/
        @. p.dalec811.next_labile_pool = p.dalec811.next_labile_pool - (p.dalec811.labile_fire_combust + p.dalec811.labile_fire_transfer) * model.temporal.DELTA_T
        @. p.dalec811.next_foliar_pool = p.dalec811.next_foliar_pool - (p.dalec811.foliar_fire_combust + p.dalec811.foliar_fire_transfer) * model.temporal.DELTA_T
        @. p.dalec811.next_root_pool = p.dalec811.next_root_pool - (p.dalec811.root_fire_combust + p.dalec811.root_fire_transfer) * model.temporal.DELTA_T
        @. p.dalec811.next_wood_pool = p.dalec811.next_wood_pool - (p.dalec811.wood_fire_combust + p.dalec811.wood_fire_transfer) * model.temporal.DELTA_T
        @. p.dalec811.next_litter_pool = p.dalec811.next_litter_pool + (p.dalec811.labile_fire_transfer + p.dalec811.foliar_fire_transfer + p.dalec811.root_fire_transfer - p.dalec811.litter_fire_combust - p.dalec811.litter_fire_transfer) * model.temporal.DELTA_T
        @. p.dalec811.next_som_pool = p.dalec811.next_som_pool + (p.dalec811.wood_fire_transfer + p.dalec811.litter_fire_transfer - p.dalec811.som_fire_combust) * model.temporal.DELTA_T
    
        @. p.dalec811.total_fire_combust = p.dalec811.labile_fire_combust + p.dalec811.foliar_fire_combust + p.dalec811.root_fire_combust + p.dalec811.wood_fire_combust + p.dalec811.litter_fire_combust + p.dalec811.som_fire_combust
    
        @. p.dalec811.nee = -p.dalec811.GPP + p.dalec811.respiration_auto + p.dalec811.respiration_hetero_litter + p.dalec811.respiration_hetero_som + p.dalec811.total_fire_combust
        
        @. dY = p - Y
    end
end

include("./util/load.jl")

end

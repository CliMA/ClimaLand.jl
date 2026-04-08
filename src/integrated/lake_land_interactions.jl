"""
    set_lake_fraction!(p, model)

Sets the lake fraction p.lake_fraction to zero if no lake is modeled,
or to the lake mask from the lake model if the model is not nothing; for 
SoilCanopyModels, no lake model is possible.
"""
set_lake_fraction!(p, model::LandModel) =
    isnothing(model.lake) ? p.lake_fraction .= 0 :
    p.lake_fraction .= model.lake.inland_water_mask
set_lake_fraction!(p, model::SoilCanopyModel) = nothing


"""
    update_lake_sediment_heat_flux!(p, lake::InlandWater.SlabLakeModel, soil)

Compute the sediment heat flux between the lake and the top soil layer.
"""
function update_lake_sediment_heat_flux!(
    p,
    lake::InlandWater.SlabLakeModel,
    soil,
)
    T_soil_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.T)
    κ_soil_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.κ)
    Δz_top = soil.domain.fields.Δz_top
    params = lake.parameters
    @. p.lake.sediment_heat_flux = InlandWater.lake_sediment_heat_flux(
        p.lake.T,
        T_soil_sfc,
        κ_soil_sfc,
        Δz_top,
        params,
    )
    return nothing
end

"""
    mask_biomass!(
        p,
        prognostic_land_components::Union{
            Val{(:canopy, :lake, :snow, :soil, :soilco2)},
            Val{(:canopy, :lake, :snow, :soil)},},
    )

Mask out areas where there is a lake by setting LAI, RAI, and SAI to zero.
Called in canopy update_aux!.

Currently treats the lake mask as 1 or 0. It is TBD how to clip LAI, SAI RAI - or if it
is needed - if the lake fraction is not binary.
"""
function Canopy.mask_biomass!(
    p,
    prognostic_land_components::Union{
        Val{(:canopy, :lake, :snow, :soil, :soilco2)},
        Val{(:canopy, :lake, :snow, :soil)},
    },
)
    canopy_mask = p.lake_fraction
    FT = eltype(p.canopy.biomass.area_index.leaf)
    @. p.canopy.biomass.area_index.leaf =
        ifelse(canopy_mask == 1, FT(0), p.canopy.biomass.area_index.leaf)
    @. p.canopy.biomass.area_index.stem =
        ifelse(canopy_mask == 1, FT(0), p.canopy.biomass.area_index.stem)
    @. p.canopy.biomass.area_index.root =
        ifelse(canopy_mask == 1, FT(0), p.canopy.biomass.area_index.root)
end

"""
    maximum_snow_cover_fraction!(
        p,
        prognostic_land_components::Union{
            Val{(:canopy, :lake, :snow, :soil, :soilco2)},
            Val{(:canopy, :lake, :snow, :soil)},},
    )

Sets the maximum snow cover fraction to be 1 - lake fraction.
"""
function Snow.maximum_snow_cover_fraction!(
    p,
    prognostic_land_components::Union{
        Val{(:canopy, :lake, :snow, :soil, :soilco2)},
        Val{(:canopy, :lake, :snow, :soil)},
    },
)
    @. p.snow.snow_cover_fraction *= (1 - p.lake_fraction)
end

"""
    lake_boundary_fluxes!(bc, prognostic_land_components, model, Y, p, t)

Integrated-mode method for lake boundary fluxes. In integrated mode,
`p.lake.R_n` is already computed by `lsm_radiant_energy_fluxes!`,
so we skip the standalone `net_radiation!` call.
"""
function InlandWater.lake_boundary_fluxes!(
    bc::InlandWater.AtmosDrivenLakeBC,
    prognostic_land_components,
    model::InlandWater.SlabLakeModel{FT},
    Y,
    p,
    t,
) where {FT}
    # Compute turbulent fluxes
    turbulent_fluxes!(p.lake.turbulent_fluxes, bc.atmos, model, Y, p, t)
    # R_n already computed by lsm_radiant_energy_fluxes! — do NOT call net_radiation!

    earth_param_set = model.parameters.earth_param_set
    @. p.lake.runoff = -(
        p.drivers.P_liq + p.drivers.P_snow + p.lake.turbulent_fluxes.vapor_flux
    )
    # p.lake.runoff is already area-weighted by lake fraction, so the runoff
    # energy flux should use it directly.
    @. p.lake.runoff_energy_flux = InlandWater.lake_runoff_energy_flux(
        p.lake.runoff,
        p.lake.T,
        p.lake.q_l,
        earth_param_set,
    )
    return nothing
end

"""
    bare_soil_fraction(p, snow, lake::SlabLakeModel)

Returns the bare soil fraction of 1-snow cover fraction-lake fraction.
"""
bare_soil_fraction(p, snow, lake::SlabLakeModel) =
    @. lazy(1 - p.snow.snow_cover_fraction - p.lake_fraction)

"""
    update_ground_albedo_PAR!(p, Y, soil, snow, lake::SlabLakeModel)

Computes the ground albedo in the PAR band when a lake is present.
"""
function update_ground_albedo_PAR!(p, Y, soil, snow, lake::SlabLakeModel)
    snow_frac = p.snow.snow_cover_fraction
    α_soil = p.soil.PAR_albedo
    α_snow = p.snow.α_snow
    f_lake = p.lake_fraction
    α_lake = p.lake.albedo
    @. p.α_ground.PAR =
        p.bare_soil_fraction * α_soil + snow_frac * α_snow + f_lake * α_lake
end

"""
    update_ground_albedo_NIR!(p, Y, soil, snow, lake::SlabLakeModel)

Computes the ground albedo in the NIR band when a lake is present.
"""
function update_ground_albedo_NIR!(p, Y, soil, snow, lake::SlabLakeModel)
    snow_frac = p.snow.snow_cover_fraction
    α_soil = p.soil.NIR_albedo
    α_snow = p.snow.α_snow
    f_lake = p.lake_fraction
    α_lake = p.lake.albedo
    @. p.α_ground.NIR =
        p.bare_soil_fraction * α_soil + snow_frac * α_snow + f_lake * α_lake
end

"""
    update_soil_heat_flux_with_lake_sediment_flux!(energy_bc, p, lake::SlabLakeModel)

Alters the energy flux boundary condition of the soil model when a sediment heat flux
from a lake is modelled.
"""
function update_soil_heat_flux_with_lake_sediment_flux!(
    energy_bc,
    p,
    lake::SlabLakeModel,
)
    @. energy_bc += p.lake_fraction * p.lake.sediment_heat_flux
end

"""
    update_lake_radiative_fluxes!(lake::SlabLakeModel, p, LW_u_lake, LW_d_canopy, _σ)
"""
function update_lake_radiative_fluxes!(
    lake::SlabLakeModel,
    p,
    LW_u_lake,
    LW_d_canopy,
    _σ,
)
    α_lake = p.lake.albedo
    R_net_lake = p.lake.R_n
    par_d = p.canopy.radiative_transfer.par_d
    nir_d = p.canopy.radiative_transfer.nir_d
    f_trans_par = p.canopy.radiative_transfer.par.trans
    f_trans_nir = p.canopy.radiative_transfer.nir.trans
    @. R_net_lake = -(
        f_trans_nir * nir_d * (1 - α_lake) + f_trans_par * par_d * (1 - α_lake)
    )

    ϵ_lake = lake.parameters.emissivity
    T_lake = p.lake.T
    @. LW_u_lake = ϵ_lake * _σ * T_lake^4 + (1 - ϵ_lake) * LW_d_canopy
    @. R_net_lake -= ϵ_lake * LW_d_canopy - ϵ_lake * _σ * T_lake^4
    return nothing
end

"""
    ground_lw_upwelling(lake::SlabLakeModel, p, LW_u_soil, LW_u_snow, LW_u_lake)

Return the area-weighted ground LW upwelling flux when a lake is modelled,
"""
function ground_lw_upwelling(
    lake::SlabLakeModel,
    p,
    LW_u_soil,
    LW_u_snow,
    LW_u_lake,
)
    snow_frac = p.snow.snow_cover_fraction
    f_lake = p.lake_fraction
    return @. lazy(
        p.bare_soil_fraction * LW_u_soil +
        snow_frac * LW_u_snow +
        f_lake * LW_u_lake,
    )
end

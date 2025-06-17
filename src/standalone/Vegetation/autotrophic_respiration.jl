export AutotrophicRespirationParameters, AutotrophicRespirationModel

abstract type AbstractAutotrophicRespirationModel{FT} <:
              AbstractCanopyComponent{FT} end

"""
    AutotrophicRespirationParameters{FT<:AbstractFloat}

The required parameters for the autrophic respiration model, which is based 
off of the JULES model.
Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AutotrophicRespirationParameters{FT <: AbstractFloat}
    "Vcmax25 to N factor (mol CO2 m-2 s-1 kg C (kg C)-1)"
    ne::FT
    "Live stem wood coefficient (kg C m-3)"
    ηsl::FT
    "Specific leaf density (kg C m-2 [leaf])"
    σl::FT
    "Ratio root nitrogen to top leaf nitrogen (-), typical value 1.0"
    μr::FT
    "Ratio stem nitrogen to top leaf nitrogen (-), typical value 0.1"
    μs::FT
    "Relative contribution or Rgrowth (-)"
    Rel::FT
end

Base.eltype(::AutotrophicRespirationParameters{FT}) where {FT} = FT

"""
    AutotrophicRespirationModel{FT, ARP <: AutotrophicRespirationParameters{FT},} <: AbstractAutotrophicRespirationModel{FT} 

The JULES autotrophic respiration model.

Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
"""
struct AutotrophicRespirationModel{
    FT,
    ARP <: AutotrophicRespirationParameters{FT},
} <: AbstractAutotrophicRespirationModel{FT}
    parameters::ARP
end

function AutotrophicRespirationModel{FT}(
    parameters::AutotrophicRespirationParameters{FT},
) where {FT <: AbstractFloat}
    return AutotrophicRespirationModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.name(model::AbstractAutotrophicRespirationModel) =
    :autotrophic_respiration
ClimaLand.auxiliary_vars(model::AutotrophicRespirationModel) = (:Ra,) # Ra = Rpm + Rg
ClimaLand.auxiliary_types(model::AutotrophicRespirationModel{FT}) where {FT} =
    (FT,)
ClimaLand.auxiliary_domain_names(::AutotrophicRespirationModel) = (:surface,)



"""
    update_autotrophic_respiration!(p, Y, autotrophic_respiration::AutotrophicRespirationModel, canopy)

Computes the autotrophic respiration rate  (mol co2 m^-2 s^-1) as the sum of the plant maintenance
and growth respirations, according to the JULES model.

Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model 
description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
"""
function update_autotrophic_respiration!(
    p,
    Y,
    autotrophic_respiration::AutotrophicRespirationModel,
    canopy,
)
    hydraulics = canopy.hydraulics
    n_stem = hydraulics.n_stem
    n_leaf = hydraulics.n_leaf
    h_canopy = hydraulics.compartment_surfaces[end]
    i_end = n_stem + n_leaf
    ψ = p.canopy.hydraulics.ψ
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf
    SAI = area_index.stem
    RAI = area_index.root
    earth_param_set = canopy.parameters.earth_param_set
    grav = LP.grav(earth_param_set)
    ρ_l = LP.ρ_cloud_liq(earth_param_set)
    (; sc, pc) = canopy.photosynthesis.parameters
    (; G_Function, Ω) = canopy.radiative_transfer.parameters
    cosθs = p.drivers.cosθs
    An = p.canopy.photosynthesis.An
    Rd = p.canopy.photosynthesis.Rd

    β = @. lazy(moisture_stress(ψ.:($$i_end) * ρ_l * grav, sc, pc))
    Vcmax25 = get_Vcmax25(p, canopy.photosynthesis)
    @. p.canopy.autotrophic_respiration.Ra = compute_autrophic_respiration(
        autotrophic_respiration,
        Vcmax25,
        LAI,
        SAI,
        RAI,
        extinction_coeff(G_Function, cosθs),
        Ω,
        An,
        Rd,
        β,
        h_canopy,
    )
end

"""
    compute_autrophic_respiration(model::AutotrophicRespirationModel,
                                  Vcmax25,
                                  LAI,
                                  SAI,
                                  RAI,
                                  K,
                                  Ω,
                                  An,
                                  Rd,
                                  β,
                                  h,
                                 )

Computes the autotrophic respiration (mol co2 m^-2 s^-1) as the sum of the plant maintenance
and growth respirations, according to the JULES model.

Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
"""
function compute_autrophic_respiration(
    model::AutotrophicRespirationModel,
    Vcmax25,
    LAI,
    SAI,
    RAI,
    K,
    Ω,
    An,
    Rd,
    β,
    h,
)

    (; ne, ηsl, σl, μr, μs, Rel) = model.parameters
    Nl, Nr, Ns =
        nitrogen_content(ne, Vcmax25, LAI, SAI, RAI, ηsl, h, σl, μr, μs)
    Rpm = plant_respiration_maintenance(Rd, β, Nl, Nr, Ns)
    Rg = plant_respiration_growth(Rel, An, Rpm)
    Ra = Rpm + Rg
    return Ra * (1 - exp(-K * LAI * Ω)) / (K * Ω) # adjust to canopy level
end

Base.broadcastable(model::AutotrophicRespirationModel) = tuple(model) # this is so that @. does not broadcast on Ref(canopy.autotrophic_respiration)

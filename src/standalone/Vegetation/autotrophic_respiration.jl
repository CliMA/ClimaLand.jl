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
    "Vcmax25 (leaf level) to N factor (mol CO2 m-2 s-1 kg C (kg C)-1)"
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
    β = @. lazy(moisture_stress(ψ.:($$i_end) * ρ_l * grav, sc, pc))
    Vcmax25_leaf = get_Vcmax25_leaf(p, canopy.photosynthesis)
    Rd_leaf = get_Rd_leaf(p, canopy.photosynthesis)
    An_leaf = get_An_leaf(p, canopy.photosynthesis)
    @. p.canopy.autotrophic_respiration.Ra = compute_autrophic_respiration(
        autotrophic_respiration,
        Vcmax25_leaf,
        LAI,
        SAI,
        RAI,
        extinction_coeff(G_Function, cosθs),
        Ω,
        An_leaf,
        Rd_leaf,
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

## For interfacing with ClimaParams

"""
    AutotrophicRespirationParameters(FT; kwargs...)
    AutotrophicRespirationParameters(toml_dict; kwargs...)

Constructors for the AutotrophicRespirationParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values.
With either constructor, you can manually override any parameter via kwargs:
```julia
AutotrophicRespirationParameters(FT; ne = 99999)
AutotrophicRespirationParameters(toml_dict; ne = 99999)
```
"""
AutotrophicRespirationParameters(
    ::Type{FT};
    kwargs...,
) where {FT <: AbstractFloat} =
    AutotrophicRespirationParameters(CP.create_toml_dict(FT); kwargs...)

function AutotrophicRespirationParameters(toml_dict::CP.ParamDict; kwargs...)
    name_map = (;
        :N_factor_Vcmax25 => :ne,
        :live_stem_wood_coeff => :ηsl,
        :specific_leaf_density => :σl,
        :root_leaf_nitrogen_ratio => :μr,
        :relative_contribution_factor => :Rel,
        :stem_leaf_nitrogen_ratio => :μs,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    AutotrophicRespirationParameters{FT}(; parameters..., kwargs...)
end



"""
    nitrogen_content(
                     ne::FT, # Mean leaf nitrogen concentration (kg N (kg C)-1)
                     Vcmax25::FT, #
                     LAI::FT, # Leaf area index
                     SAI::FT,
                     RAI::FT,
                     ηsl::FT, # live stem  wood coefficient (kg C m-3)
                     h::FT, # canopy height (m)
                     σl::FT # Specific leaf density (kg C m-2 [leaf])
                     μr::FT, # Ratio root nitrogen to top leaf nitrogen (-), typical value 1.0
                     μs::FT, # Ratio stem nitrogen to top leaf nitrogen (-), typical value 0.1
                    ) where {FT}

Computes the nitrogen content of leafs (Nl), roots (Nr) and stems (Ns).
"""
function nitrogen_content(
    ne::FT, # Mean leaf nitrogen concentration (kg N (kg C)-1)
    Vcmax25::FT, #
    LAI::FT, # Leaf area index
    SAI::FT,
    RAI::FT,
    ηsl::FT, # live stem  wood coefficient (kg C m-3)
    h::FT, # canopy height (m)
    σl::FT, # Specific leaf density (kg C m-2 [leaf])
    μr::FT, # Ratio root nitrogen to top leaf nitrogen (-), typical value 1.0
    μs::FT, # Ratio stem nitrogen to top leaf nitrogen (-), typical value 0.1
) where {FT}
    Sc = ηsl * h * LAI * ClimaLand.heaviside(SAI)
    Rc = σl * RAI
    nm = Vcmax25 / ne
    Nl = nm * σl * LAI
    Nr = μr * nm * Rc
    Ns = μs * nm * Sc
    return Nl, Nr, Ns
end

"""
    plant_respiration_maintenance(
        Rd::FT, # Dark respiration
        β::FT, # Soil moisture factor
        Nl::FT, # Nitrogen content of leafs
        Nr::FT, # Nitrogen content of roots
        Ns::FT, # Nitrogen content of stems
        ) where {FT}

Computes plant maintenance respiration as a function of dark respiration (Rd),
the nitrogen content of leafs (Nl), roots (Nr) and stems (Ns),
and the soil moisture factor (β).
"""
function plant_respiration_maintenance(
    Rd::FT, # Dark respiration
    β::FT, # Soil moisture factor
    Nl::FT, # Nitrogen content of leafs
    Nr::FT, # Nitrogen content of roots
    Ns::FT, # Nitrogen content of stems
) where {FT}
    # When LAI is zero, Nl = 0
    Rpm = Rd * (β + (Nr + Ns) / max(Nl, eps(FT)))
    return Rpm
end

"""
    plant_respiration_growth(
        Rel::FT, # Factor of relative contribution
        An::FT, # Net photosynthesis
        Rpm::FT # Plant maintenance respiration
        ) where {FT}

Computes plant growth respiration as a function of net photosynthesis (An),
plant maintenance respiration (Rpm), and a relative contribution factor, Rel.
"""
function plant_respiration_growth(Rel::FT, An::FT, Rpm::FT) where {FT}
    Rg = Rel * (An - Rpm)
    return Rg
end

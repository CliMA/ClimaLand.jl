export AutotrophicRespirationParameters, AutotrophicRespirationModel

abstract type AbstractAutotrophicRespirationModel{FT} <:
              AbstractCanopyComponent{FT} end

"""
    AutotrophicRespirationParameters{FT<:AbstractFloat}

The required parameters for the autrophic respiration model, inspired by the
JULES model.
Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AutotrophicRespirationParameters{FT <: AbstractFloat}
    "Live stem wood coefficient (kg C m-3)"
    ηsl::FT
    "Specific leaf density (kg C m-2 [leaf])"
    σl::FT
    "Ratio stem nitrogen to root nitrogen (-)"
    μs::FT
    "Relative contribution or Rgrowth (-)"
    Rel::FT
    "Q10 temperature sensitivity of maintenance respiration (-), JULES default 2.0"
    Q10::FT
    "Reference temperature for the Q10 maintenance-respiration factor (K), JULES default 298.15"
    T_ref::FT
    "Reference leaf-level dark respiration rate (mol CO2 m-2 s-1) used as the base rate for the LAI-independent root and stem maintenance respiration"
    Rd_ref::FT
end

Base.eltype(::AutotrophicRespirationParameters{FT}) where {FT} = FT

"""
    AutotrophicRespirationModel{FT, ARP <: AutotrophicRespirationParameters{FT},} <: AbstractAutotrophicRespirationModel{FT}

Autotrophic respiration model inspired by JULES.

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
and growth respirations, inspired by the JULES model.

Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model
description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
"""
function update_autotrophic_respiration!(
    p,
    Y,
    autotrophic_respiration::AutotrophicRespirationModel,
    canopy,
)
    area_index = p.canopy.biomass.area_index
    SAI = area_index.stem
    RAI = area_index.root
    Rd_canopy = get_Rd_canopy(p, canopy.photosynthesis)
    An_canopy = get_An_canopy(p, canopy.photosynthesis)
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    @. p.canopy.autotrophic_respiration.Ra = compute_autrophic_respiration(
        autotrophic_respiration,
        SAI,
        RAI,
        An_canopy,
        Rd_canopy,
        T_canopy,
    )
end

"""
    compute_autrophic_respiration(model::AutotrophicRespirationModel,
                                  SAI,
                                  RAI,
                                  An_canopy,
                                  Rd_canopy,
                                  T_canopy,
                                 )

Computes the autotrophic respiration (mol co2 m^-2 s^-1) as the sum of the plant maintenance
and growth respirations, inspired by the JULES model.

Maintenance respiration is split by tissue. The leaf term uses `Rd_canopy`
(∝ LAI) so it vanishes when leafless; the root and stem terms use a reference
rate `Rd_ref` scaled by the time-constant area indices RAI/SAI, so they persist
through winter. The sum is scaled by the Q10 factor `Q10^((T_canopy - T_ref)/10)`.

Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
"""
function compute_autrophic_respiration(
    model::AutotrophicRespirationModel,
    SAI,
    RAI,
    An_canopy,
    Rd_canopy,
    T_canopy,
)

    (; μs, Rel, Q10, T_ref, Rd_ref) = model.parameters
    FT = typeof(T_canopy)
    f_T = Q10^((T_canopy - T_ref) / FT(10))
    # Rd_canopy already carries the moisture-stress weighting, so R_leaf needs no
    # extra β factor.
    R_leaf = Rd_canopy
    R_root = Rd_ref * RAI
    R_stem = Rd_ref * μs * SAI
    Rpm = f_T * (R_leaf + R_root + R_stem)
    Rg = plant_respiration_growth(Rel, An_canopy, Rpm)
    Ra = Rpm + Rg
    return Ra # already canopy level
end

Base.broadcastable(model::AutotrophicRespirationModel) = tuple(model) # this is so that @. does not broadcast on Ref(canopy.autotrophic_respiration)

## For interfacing with ClimaParams

"""
    AutotrophicRespirationParameters(toml_dict::CP.ParamDict; kwargs...)

Constructor for the `AutotrophicRespirationParameters` struct by passing a TOML
dictionary.
You can manually override any parameter via keyword arguments:
```julia
AutotrophicRespirationParameters(toml_dict; μs = 0.1)
```
"""
function AutotrophicRespirationParameters(
    toml_dict::CP.ParamDict;
    ηsl = toml_dict["live_stem_wood_coeff"],
    σl = toml_dict["specific_leaf_density"],
    Rel = toml_dict["relative_contribution_factor"],
    μs = toml_dict["stem_leaf_nitrogen_ratio"],
    Q10 = toml_dict["autotrophic_respiration_Q10"],
    T_ref = toml_dict["autotrophic_respiration_T_ref"],
    Rd_ref = toml_dict["autotrophic_respiration_Rd_ref"],
)
    FT = CP.float_type(toml_dict)
    AutotrophicRespirationParameters{FT}(; ηsl, σl, Rel, μs, Q10, T_ref, Rd_ref)
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
    # Growth respiration cannot be negative.
    Rg = Rel * max(An - Rpm, FT(0))
    return Rg
end

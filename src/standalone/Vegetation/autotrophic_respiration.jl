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
    h_canopy = canopy.biomass.height
    area_index = p.canopy.biomass.area_index
    LAI = area_index.leaf
    SAI = area_index.stem
    RAI = area_index.root
    β = p.canopy.soil_moisture_stress.βm
    Vcmax25_canopy = get_Vcmax25_canopy(p, canopy.photosynthesis)
    Rd_canopy = get_Rd_canopy(p, canopy.photosynthesis)
    An_canopy = get_An_canopy(p, canopy.photosynthesis)
    # Drive the Q10 factor with air temperature, not the prognostic canopy
    # temperature: where LAI+SAI ≈ 0 the latter can diverge and overflow
    # Q10^((T-T_ref)/10) to Inf, giving Ra = NaN. T_air stays bounded.
    T_air = p.drivers.T
    @. p.canopy.autotrophic_respiration.Ra = compute_autrophic_respiration(
        autotrophic_respiration,
        Vcmax25_canopy,
        LAI,
        SAI,
        RAI,
        An_canopy,
        Rd_canopy,
        β,
        h_canopy,
        T_air,
    )
end

"""
    compute_autrophic_respiration(model::AutotrophicRespirationModel,
                                  Vcmax25,
                                  LAI,
                                  SAI,
                                  RAI,
                                  An_canopy,
                                  Rd_canopy,
                                  β,
                                  h,
                                  T_air,
                                 )

Computes the autotrophic respiration (mol co2 m^-2 s^-1) as the sum of the plant maintenance
and growth respirations, based on the JULES model.

Maintenance respiration is split by tissue. The leaf term uses `Rd_canopy`
(∝ LAI) so it vanishes when leafless; the root and stem terms use a reference
rate `Rd_ref` scaled by the time-constant area indices RAI/SAI, so they persist
through winter. This fixes the original JULES formulation, where all of `Rpm`
scaled with `Rd_canopy ∝ LAI` and collapsed to ~0 whenever LAI → 0.

The whole term is scaled by the Q10 factor `Q10^((T_air - T_ref)/10)`, driven by
air temperature (see `update_autotrophic_respiration!` for why not canopy temp).

Clark, D. B., et al. "The Joint UK Land Environment Simulator (JULES), model description–Part 2: carbon fluxes and vegetation dynamics." Geoscientific Model Development 4.3 (2011): 701-722.
"""
function compute_autrophic_respiration(
    model::AutotrophicRespirationModel,
    Vcmax25_canopy,
    LAI,
    SAI,
    RAI,
    An_canopy,
    Rd_canopy,
    β,
    h,
    T_air,
)

    (; ηsl, σl, μr, μs, Rel, Q10, T_ref, Rd_ref) = model.parameters
    FT = typeof(T_air)
    f_T = Q10^((T_air - T_ref) / FT(10))
    # Leaf maintenance scales with Rd_canopy (∝ LAI); root/stem use Rd_ref and
    # the time-constant RAI/SAI so they persist when LAI → 0. μr/μs are the
    # root/stem-to-leaf nitrogen ratios; ηsl*h*SAI is the live stem carbon.
    R_leaf = Rd_canopy * β
    R_root = Rd_ref * μr * RAI * β
    R_stem = Rd_ref * μs * (ηsl * h / σl) * SAI
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
AutotrophicRespirationParameters(toml_dict; ne = 99999)
```
"""
function AutotrophicRespirationParameters(
    toml_dict::CP.ParamDict;
    ne = toml_dict["N_factor_Vcmax25"],
    ηsl = toml_dict["live_stem_wood_coeff"],
    σl = toml_dict["specific_leaf_density"],
    μr = toml_dict["root_leaf_nitrogen_ratio"],
    Rel = toml_dict["relative_contribution_factor"],
    μs = toml_dict["stem_leaf_nitrogen_ratio"],
    Q10 = toml_dict["autotrophic_respiration_Q10"],
    T_ref = toml_dict["autotrophic_respiration_T_ref"],
    Rd_ref = toml_dict["autotrophic_respiration_Rd_ref"],
)
    FT = CP.float_type(toml_dict)
    AutotrophicRespirationParameters{FT}(;
        ne,
        ηsl,
        σl,
        μr,
        Rel,
        μs,
        Q10,
        T_ref,
        Rd_ref,
    )
end



"""
    nitrogen_content(
                     ne::FT, # Mean leaf nitrogen concentration (kg N (kg C)-1)
                     Vcmax25_canopy::FT, #
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
    Vcmax25_canopy::FT, #
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
    nm = Vcmax25_canopy / (ne * max(LAI, sqrt(eps(FT))))
    Nl = nm * σl
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
    # Clamp at zero so a positive winter maintenance baseline (An ≈ 0) can't
    # drive growth respiration, and hence Ra, negative.
    Rg = Rel * max(An - Rpm, FT(0))
    return Rg
end

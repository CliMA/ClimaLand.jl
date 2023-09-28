export AutotrophicRespirationParameters, AutotrophicRespirationModel

abstract type AbstractAutotrophicRespirationModel{FT} <:
              AbstractCanopyComponent{FT} end

"""
    AutotrophicRespirationParameters{FT<:AbstractFloat}

The required parameters for the autrophic respiration model.
$(DocStringExtensions.FIELDS)
"""
struct AutotrophicRespirationParameters{FT <: AbstractFloat}
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
    "Factor to convert from mol CO2 to kg C" # Definitely not a parameter. not sure what to do.
    f1::FT
    "Factor of relative contribution or Rgrowth (-)"
    f2::FT
end

"""
    function AutotrophicRespirationParameters{FT}(;
        ne = FT(8 * 1e-4),
        ηsl = FT(0.01),
        σl = FT(0.05),
        μr = FT(1.0),
        μs = FT(0.1),
        f1 = FT(0.012), 
        f2 = FT(0.25)        
) where {FT}

A constructor supplying default values for the AutotrophicRespirationParameters struct.
"""
function AutotrophicRespirationParameters{FT}(;
    ne = FT(8 * 1e-4),
    ηsl = FT(0.01),
    σl = FT(0.05),
    μr = FT(1.0),
    μs = FT(0.1),
    f1 = FT(0.012),
    f2 = FT(0.25),
) where {FT}
    return AutotrophicRespirationParameters{FT}(ne, ηsl, σl, μr, μs, f1, f2)
end

struct AutotrophicRespirationModel{FT} <:
       AbstractAutotrophicRespirationModel{FT} # we could give it a more specific name...
    parameters::AutotrophicRespirationParameters{FT}
end

ClimaLSM.name(model::AbstractAutotrophicRespirationModel) =
    :autotrophic_respiration
ClimaLSM.auxiliary_vars(model::AutotrophicRespirationModel) = (:Ra,) # Ra = Rpm + Rg
ClimaLSM.auxiliary_types(model::AutotrophicRespirationModel{FT}) where {FT} =
    (FT,)
ClimaLSM.auxiliary_domain_names(::AutotrophicRespirationModel) = (:surface,)

"""

"""
function compute_autrophic_respiration(
    model::AutotrophicRespirationModel,
    Vcmax25,
    LAI,
    RAI,
    GPP,
    Rd,
    β,
    h,
)

    (; ne, ηsl, σl, μr, μs, f1, f2) = model.parameters

    Nl, Nr, Ns = nitrogen_content(ne, Vcmax25, LAI, RAI, ηsl, h, σl, μr, μs)
    Rpm = plant_respiration_maintenance(Rd, β, Nl, Nr, Ns, f1)
    Rg = plant_respiration_growth(f2, GPP, Rpm)
    Ra = Rpm + Rg # Should this be a function in canopy_parameterizations.jl or is it ok here?
    return Ra
end

Base.broadcastable(model::AutotrophicRespirationModel) = tuple(model) # this is so that @. does not broadcast on Ref(canopy.autotrophic_respiration)

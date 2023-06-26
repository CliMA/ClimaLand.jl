export AbstractSoilHydrologyClosure, vanGenuchten, BrooksCorey

"""
    AbstractSoilHydrologyClosure{FT <: AbstractFloat} 

The abstract type of soil hydrology closure, of which
vanGenuchten{FT} and BrooksCorey{FT} are the two supported 
concrete types.

To add a new parameterization, methods are required for:
- matric_potential,
- inverse_matric_potential,
- pressure_head,
- dψdϑ,
- hydraulic_conductivity.
"""
abstract type AbstractSoilHydrologyClosure{FT <: AbstractFloat} end
Base.broadcastable(x::AbstractSoilHydrologyClosure) = tuple(x)

"""
    vanGenuchten{FT} <: AbstractSoilHydrologyClosure{FT}

The van Genuchten soil hydrology closure, chosen when the 
hydraulic conductivity and matric potential are modeled
using the van Genuchten parameterization (van Genuchten 1980;
see also Table 8.2 of G. Bonan 2019).

$(DocStringExtensions.FIELDS)
"""
struct vanGenuchten{FT} <: AbstractSoilHydrologyClosure{FT}
    "The inverse of the air entry potential (1/m)"
    α::FT
    "The van Genuchten pore-size distribution index (unitless)"
    n::FT
    "The van Genuchten parameter m = 1 - 1/n (unitless)"
    m::FT
    "A derived parameter: the critical saturation at which capillary flow no longer replenishes the surface"
    S_c::FT
    function vanGenuchten(; α::FT, n::FT) where {FT}
        m = 1 - 1 / n
        S_c = (1 + ((n - 1) / n)^(1 - 2 * n))^(-m)
        return new{FT}(α, n, m, S_c)
    end

end

"""
   BrooksCorey{FT} <: AbstractSoilHydrologyClosure{FT}

The Brooks and Corey soil hydrology closure, chosen when the 
hydraulic conductivity and matric potential are modeled
using the Brooks and Corey parameterization (Brooks and Corey,
1964, 1966; see also Table 8.2 of G. Bonan 2019).

$(DocStringExtensions.FIELDS)
"""
struct BrooksCorey{FT} <: AbstractSoilHydrologyClosure{FT}
    "The pore-size distribution index (unitless)"
    c::FT
    "The air entry matric potential, when S=1 (m)"
    ψb::FT
    "A derived parameter: the critical saturation at which capillary flow no longer replenishes the surface"
    S_c::FT
    function BrooksCorey(; c::FT, ψb::FT) where {FT}
        S_c = (1 + 1 / c)^(-c)
        return new{FT}(c, ψb, S_c)
    end
end

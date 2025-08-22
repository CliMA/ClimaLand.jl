export FarquharParameters, FarquharModel

abstract type AbstractPhotosynthesisModel{FT} <: AbstractCanopyComponent{FT} end

"""
    FarquharParameters{
        FT<:AbstractFloat,
        MECH <: Union{FT, ClimaCore.Fields.Field},
        VC <: Union{FT, ClimaCore.Fields.Field},
    }

The required parameters for the Farquhar photosynthesis model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct FarquharParameters{
    FT <: AbstractFloat,
    MECH <: Union{FT, ClimaCore.Fields.Field},
    VC <: Union{FT, ClimaCore.Fields.Field},
}
    "Vcmax at 25 °C (mol CO2/m^2/s; per leaf area)"
    Vcmax25::VC
    "Γstar at 25 °C (mol/mol)"
    Γstar25::FT
    "Michaelis-Menten parameter for CO2 at 25 °C (mol/mol)"
    Kc25::FT
    "Michaelis-Menten parameter for O2 at 25 °C (mol/mol)"
    Ko25::FT
    "Energy of activation for CO2 (J/mol)"
    ΔHkc::FT
    "Energy of activation for oxygen (J/mol)"
    ΔHko::FT
    "Energy of activation for Vcmax (J/mol)"
    ΔHVcmax::FT
    "Energy of activation for Γstar (J/mol)"
    ΔHΓstar::FT
    "Energy of activation for Jmax (J/mol)"
    ΔHJmax::FT
    "Energy of activation for Rd (J/mol)"
    ΔHRd::FT
    "Reference temperature equal to 25 degrees Celsius (K)"
    To::FT
    "Intercelluar O2 concentration (mol/mol); taken to be constant"
    oi::FT
    "Quantum yield of photosystem II (Bernacchi, 2003; unitless)"
    ϕ::FT
    "Curvature parameter, a fitting constant to compute J, unitless"
    θj::FT
    "Constant factor appearing the dark respiration term for C3, equal to 0.015."
    fC3::FT
    "Constant factor appearing the dark respiration term for C4, equal to 0.025."
    fC4::FT
    "Sensitivity to low water pressure, in the moisture stress factor, (Pa^{-1}) [Tuzet et al. (2003)]"
    sc::FT
    "Reference water pressure for the moisture stress factor (Pa) [Tuzet et al. (2003)]"
    pc::FT
    "Q10 temperature parameter for Vcmax and Rd for C4 photosynthesis; unitless"
    Q10::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K^{-1}"
    s1::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K"
    s2::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K^{-1}"
    s3::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K"
    s4::FT
    "Parameter appearing in temperature dependence of C4 Rd; K^{-1}"
    s5::FT
    "Parameter appearing in temperature dependence of C4 Rd; K"
    s6::FT
    "Quantum yield for C4 photosynthesis; mol/mol"
    E::FT
    "Photosynthesis mechanism: 1.0 indicates C3, 0.0 indicates C4"
    is_c3::MECH
end

Base.eltype(::FarquharParameters{FT}) where {FT} = FT

struct FarquharModel{FT, FP <: FarquharParameters{FT}} <:
       AbstractPhotosynthesisModel{FT}
    parameters::FP
end

function FarquharModel{FT}(
    parameters::FarquharParameters{FT},
) where {FT <: AbstractFloat}
    return FarquharModel{eltype(parameters), typeof(parameters)}(parameters)
end

ClimaLand.name(model::AbstractPhotosynthesisModel) = :photosynthesis
ClimaLand.auxiliary_vars(model::FarquharModel) = (:An, :GPP, :Rd)
ClimaLand.auxiliary_types(model::FarquharModel{FT}) where {FT} = (FT, FT, FT)
ClimaLand.auxiliary_domain_names(::FarquharModel) =
    (:surface, :surface, :surface)


function gross_leaf_photosynthesis_at_a_point_Farquhar(
    T,
    APAR_leaf_moles,
    c_co2,
    medlyn_factor,
    R,
    Vcmax25_leaf,
    is_c3,
    Γstar25,
    ΔHJmax,
    ΔHVcmax,
    ΔHΓstar,
    fC3,
    fC4,
    ΔHRd,
    To,
    θj,
    ϕ,
    oi,
    Kc25,
    Ko25,
    ΔHkc,
    ΔHko,
    Q10,
    s1,
    s2,
    s3,
    s4,
    s5,
    s6,
    E,
)
    Jmax = max_electron_transport(Vcmax25_leaf, ΔHJmax, T, To, R)
    J = electron_transport(APAR_leaf_moles, Jmax, θj, ϕ)
    Vcmax_leaf = compute_Vcmax(
        is_c3,
        Vcmax25_leaf,
        T,
        R,
        To,
        ΔHVcmax,
        Q10,
        s1,
        s2,
        s3,
        s4,
    )
    Γstar = co2_compensation(Γstar25, ΔHΓstar, T, To, R)
    ci = intercellular_co2(c_co2, Γstar, medlyn_factor)
    Aj = light_assimilation(is_c3, J, ci, Γstar, APAR_leaf_moles, E)
    Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko(Ko25, ΔHko, T, To, R)
    Ac = rubisco_assimilation(is_c3, Vcmax_leaf, ci, Γstar, Kc, Ko, oi)
    return gross_photosynthesis(Ac, Aj)
end

"""
    update_photosynthesis!(
        p,
        Y,
        model::FarquharModel,
        canopy
)

Computes the net leaf photosynthesis rate `An` (mol CO2/m^2/s) for the Farquhar model, along with the
dark respiration `Rd` (mol CO2/m^2/s) and gross primary productivity `GPP` (mol CO2/m^2/s), 
and updates them in place.
"""
function update_photosynthesis!(p, Y, model::FarquharModel, canopy)
    (;
        Vcmax25,
        is_c3,
        Γstar25,
        ΔHJmax,
        ΔHVcmax,
        ΔHΓstar,
        fC3,
        fC4,
        ΔHRd,
        To,
        θj,
        ϕ,
        oi,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        Q10,
        s1,
        s2,
        s3,
        s4,
        s5,
        s6,
        E,
    ) = model.parameters

    # unpack a bunch of stuff from p and params
    Rd = p.canopy.photosynthesis.Rd
    An = p.canopy.photosynthesis.An
    GPP = p.canopy.photosynthesis.GPP
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    f_abs = p.canopy.radiative_transfer.par.abs
    ψ = p.canopy.hydraulics.ψ
    c_co2_air = p.drivers.c_co2
    cosθs = p.drivers.cosθs
    P_air = p.drivers.P
    T_air = p.drivers.T
    q_air = p.drivers.q
    earth_param_set = canopy.parameters.earth_param_set
    lightspeed = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    grav = LP.grav(earth_param_set)
    ρ_l = LP.ρ_cloud_liq(earth_param_set)
    R = LP.gas_constant(earth_param_set)
    thermo_params = earth_param_set.thermo_params
    (; G_Function, λ_γ_PAR, Ω) = canopy.radiative_transfer.parameters
    (; sc, pc) = canopy.photosynthesis.parameters
    (; g1,) = canopy.conductance.parameters
    n_stem = canopy.hydraulics.n_stem
    n_leaf = canopy.hydraulics.n_leaf
    i_end = n_stem + n_leaf
    par_d = p.canopy.radiative_transfer.par_d
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf

    β = @. lazy(moisture_stress(ψ.:($$i_end) * ρ_l * grav, sc, pc))
    medlyn_factor = @. lazy(medlyn_term(g1, T_air, P_air, q_air, thermo_params))

    @. Rd = dark_respiration(
        is_c3,
        Vcmax25,
        β,
        T_canopy,
        R,
        To,
        fC3,
        ΔHRd,
        Q10,
        s5,
        s6,
        fC4,
    )
    APAR_leaf_moles = @. lazy(
        compute_APAR_leaf_moles(
            f_abs,
            par_d,
            λ_γ_PAR,
            lightspeed,
            planck_h,
            N_a,
            LAI,
        ),
    )
    A_leaf = @. lazy(
        gross_leaf_photosynthesis_at_a_point_Farquhar(
            T_canopy,
            APAR_leaf_moles,
            c_co2_air,
            medlyn_factor,
            R,
            Vcmax25,
            is_c3,
            Γstar25,
            ΔHJmax,
            ΔHVcmax,
            ΔHΓstar,
            fC3,
            fC4,
            ΔHRd,
            To,
            θj,
            ϕ,
            oi,
            Kc25,
            Ko25,
            ΔHkc,
            ΔHko,
            Q10,
            s1,
            s2,
            s3,
            s4,
            s5,
            s6,
            E,
        ),
    )
    # Compute net from gross
    # See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
    @. An = max(0, β * A_leaf - Rd)
    # Compute GPP: scale from leaf level A to canopy level GPP
    @. GPP = GPP_from_leaf_level_An(
        A_leaf * β,
        extinction_coeff(G_Function, cosθs),
        LAI,
        Ω,
    )

end

Base.broadcastable(m::FarquharParameters) = tuple(m)

get_Vcmax25_leaf(p, m::FarquharModel) = m.parameters.Vcmax25
get_Rd_leaf(p, m::FarquharModel) = p.canopy.photosynthesis.Rd
get_An_leaf(p, m::FarquharModel) = p.canopy.photosynthesis.An
function get_Jmax_leaf(Y, p, canopy, m::FarquharModel)
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    earth_param_set = canopy.parameters.earth_param_set
    (; Vcmax25, ΔHJmax, To) = m.parameters
    R = LP.gas_constant(canopy.parameters.earth_param_set)
    Jmax_leaf =
        @. lazy(max_electron_transport(Vcmax25, ΔHJmax, T_canopy, To, R))
    return Jmax_leaf
end


function get_J_over_Jmax(Y, p, canopy, m::FarquharModel)
    J = get_J_leaf(Y, p, canopy, m)
    Jmax = get_Jmax_leaf(Y, p, canopy, m)
    FT = eltype(p.canopy.photosynthesis.An)
    return @. lazy(J / max(Jmax, sqrt(eps(FT))))
end
"""
    max_electron_transport(Vcmax::FT) where {FT}

Computes the maximum potential rate of electron transport (`Jmax`),
in units of mol/m^2/s,
as a function of Vcmax at 25 °C (`Vcmax25`),
a constant (`ΔHJmax`), a standard temperature (`To`),
the unversal gas constant (`R`), and the temperature (`T`);
used in the Farquhar model.

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function max_electron_transport(
    Vcmax25::FT,
    ΔHJmax::FT,
    T::FT,
    To::FT,
    R::FT,
) where {FT}
    Jmax25 = Vcmax25 * FT(exp(1))
    Jmax = Jmax25 * arrhenius_function(T, To, R, ΔHJmax)
    return Jmax
end


"""
    electron_transport(APAR_leaf_moles::FT,
                       Jmax::FT,
                       θj::FT,
                       ϕ::FT) where {FT}

Computes the rate of electron transport (`J`),
in units of mol/m^2/s, as a function of
the maximum potential rate of electron transport (`Jmax`),
absorbed photosynthetically active radiation (`APAR_leaf_moles`),
an empirical "curvature parameter" (`θj`; Bonan Eqn 11.21)
and the quantum yield of photosystem II (`ϕ`);
used in the Farquhar model.

See Ch 11, G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function electron_transport(
    APAR_leaf_moles::FT,
    Jmax::FT,
    θj::FT,
    ϕ::FT,
) where {FT}
    # Light utilization of APAR (leaf level)
    IPSII = ϕ * APAR_leaf_moles / 2
    # This is a solution to a quadratic equation
    # θj *J^2 - (IPSII+Jmax)*J+IPSII*Jmax = 0, Equation 11.21
    J =
        (IPSII + Jmax - sqrt((IPSII + Jmax)^2 - 4 * θj * IPSII * Jmax)) /
        (2 * θj)
    return J
end

function get_J_leaf(Y, p, canopy, m::FarquharModel)
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    earth_param_set = canopy.parameters.earth_param_set
    f_abs_par = p.canopy.radiative_transfer.par.abs
    par_d = p.canopy.radiative_transfer.par_d
    (; λ_γ_PAR,) = canopy.radiative_transfer.parameters
    c = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf
    APAR_leaf_moles = @. lazy(
        compute_APAR_leaf_moles(
            f_abs_par,
            par_d,
            λ_γ_PAR,
            c,
            planck_h,
            N_a,
            LAI,
        ),
    )

    Jmax_leaf = get_Jmax_leaf(Y, p, canopy, m)
    (; θj, ϕ) = m.parameters
    return @. lazy(electron_transport(APAR_leaf_moles, Jmax_leaf, θj, ϕ))
end

include("./optimality_farquhar.jl")

## For interfacing with ClimaParams

"""
    function FarquharParameters(
        ::Type{FT},
        is_c3::Union{FT, ClimaCore.Fields.Field};
        Vcmax25 = FT(5e-5),
        kwargs...  # For individual parameter overrides
    )

    function FarquharParameters(
        toml_dict::CP.AbstractTOMLDict,
        is_c3::Union{AbstractFloat, ClimaCore.Fields.Field};
        Vcmax25 = FT(5e-5),
        kwargs...  # For individual parameter overrides
    )

Constructors for the FarquharParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values.Possible calls:
```julia
ClimaLand.Canopy.FarquharParameters(Float64, 1.0)
# Kwarg overrides
ClimaLand.Canopy.FarquharParameters(Float64, 1.0; Vcmax25 = 99999999, pc = 444444444)
# TOML Dictionary:
import ClimaParams as CP
toml_dict = CP.create_toml_dict(Float32);
ClimaLand.Canopy.FarquharParameters(toml_dict, 1.0f0; Vcmax25 = 99999999, pc = 444444444)
```
"""
FarquharParameters(
    ::Type{FT},
    is_c3::Union{FT, ClimaCore.Fields.Field};
    kwargs...,
) where {FT <: AbstractFloat} =
    FarquharParameters(CP.create_toml_dict(FT), is_c3; kwargs...)

function FarquharParameters(
    toml_dict::CP.AbstractTOMLDict,
    is_c3::Union{AbstractFloat, ClimaCore.Fields.Field};
    Vcmax25 = 5e-5,
    kwargs...,
)
    name_map = (;
        :Jmax_activation_energy => :ΔHJmax,
        :intercellular_O2_concentration => :oi,
        :CO2_compensation_point_25c => :Γstar25,
        :Farquhar_curvature_parameter => :θj,
        :kelvin_25C => :To,
        :photosystem_II_quantum_yield => :ϕ,
        :O2_michaelis_menten => :Ko25,
        :CO2_michaelis_menten => :Kc25,
        :dark_respiration_factor => :fC3,
        :O2_activation_energy => :ΔHko,
        :low_water_pressure_sensitivity => :sc,
        :Rd_activation_energy => :ΔHRd,
        :Vcmax_activation_energy => :ΔHVcmax,
        :Γstar_activation_energy => :ΔHΓstar,
        :CO2_activation_energy => :ΔHkc,
        :moisture_stress_ref_water_pressure => :pc,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    C4_parameters = (
        :fC4 => FT(0.025),
        :Q10 => FT(2),
        :s1 => FT(0.3),
        :s2 => FT(313.15),
        :s3 => FT(0.2),
        :s4 => FT(288.15),
        :s5 => FT(1.3),
        :s6 => FT(328.15),
        :E => FT(0.05),
    )
    # if is_c3 is a field, is_c3 may contain values between 0.0 and 1.0 after regridding
    # this deals with that possibility by rounding to the closest int
    is_c3 = max.(min.(is_c3, FT(1)), FT(0)) # placeholder
    is_c3 = round.(is_c3)
    MECH = typeof(is_c3)
    Vcmax25 = FT.(Vcmax25)
    VC = typeof(Vcmax25)
    return FarquharParameters{FT, MECH, VC}(;
        is_c3,
        Vcmax25,
        parameters...,
        C4_parameters...,
        kwargs...,
    )
end

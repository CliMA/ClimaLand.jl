export FarquharParameters, FarquharModel
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
    "Vcmax at 25 °C (mol CO2/m^2/s); leaf level"
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

"""
   FarquharModel{FT, FP <: FarquharParameters{FT}}

A photosynthesis model taking leaf-level Vcmax25 and multiple other
constants and predicting dark respiration at the leaf level, net
photosynthesis at the leaf level, and GPP at the canopy level.

The Farquhar model computes photosynthetic rates using leaf-level
properties, and so Vcmax25, Rd, An, etc are per unit area leaf. This
version also works with ci, Γstar, Ko, Ki, in units of mol/mol, and not Pa.
"""
struct FarquharModel{FT, FP <: FarquharParameters{FT}} <:
       AbstractPhotosynthesisModel{FT}
    parameters::FP
end

function FarquharModel{FT}(
    parameters::FarquharParameters{FT},
) where {FT <: AbstractFloat}
    return FarquharModel{eltype(parameters), typeof(parameters)}(parameters)
end

ClimaLand.auxiliary_vars(model::FarquharModel) = (:An, :Rd, :GPP)
ClimaLand.auxiliary_types(model::FarquharModel{FT}) where {FT} = (FT, FT, FT)
ClimaLand.auxiliary_domain_names(::FarquharModel) =
    (:surface, :surface, :surface)

"""
    rubisco_assimilation_farquhar(is_c3::AbstractFloat, args...)

Calls the correct rubisco assimilation function based on the `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_rubisco_assimilation_farquhar`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_rubisco_assimilation_farquhar`.
"""
function rubisco_assimilation_farquhar(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_rubisco_assimilation_farquhar(args...) :
    c4_rubisco_assimilation_farquhar(args...)
end

"""
    c3_rubisco_assimilation_farquhar(Vcmax::FT,
                         ci::FT,
                         Γstar::FT,
                         Kmm::FT) where {FT}

Computes the Rubisco limiting rate of photosynthesis for C3 plants (`Ac`),
in units of moles CO2/m^2/s,
as a function of the maximum rate of carboxylation of Rubisco (`Vcmax`),
the leaf internal carbon dioxide partial pressure (`ci`),
the CO2 compensation point (`Γstar`), and the effective Michaelis-Menten parameter
`K_mm`.

The units of ci, Γstar, and Kmm may be Pa or mol/mol, but they must be consistent. The
units of Vcmax must be mol CO2/m^2.s. 

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_rubisco_assimilation_farquhar(
    Vcmax::FT,
    ci::FT,
    Γstar::FT,
    Kmm::FT,
) where {FT}
    Ac = Vcmax * (ci - Γstar) / (ci + Kmm)
    return Ac
end

"""
    c4_rubisco_assimilation_farquhar(Vcmax::FT,_...) where {FT}

Computes the Rubisco limiting rate of photosynthesis for C4 plants (`Ac`)
in units of moles CO2/m^2/s,
as equal to the maximum rate of carboxylation of Rubisco (`Vcmax`).
"""
function c4_rubisco_assimilation_farquhar(Vcmax::FT, _...) where {FT}
    Ac = Vcmax
    return Ac
end

"""
    light_assimilation_farquhar(is_c3::AbstractFloat, args...)

Calls the correct light assimilation function based on the `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_light_assimilation_farquhar`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_light_assimilation_farquhar`.
"""
function light_assimilation_farquhar(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_light_assimilation_farquhar(args...) :
    c4_light_assimilation_farquhar(args...)
end
"""
    c3_light_assimilation_farquhar(
                       J::FT,
                       ci::FT,
                       Γstar::FT,
                       ::FT,
                       ::FT) where {FT}

Computes the electron transport limiting rate (`Aj`),
in units of moles CO2/m^2/s.

For C3 plants, this is a function of
the rate of electron transport (`J`), the leaf internal carbon dioxide partial pressure (`ci`),
and the CO2 compensation point (`Γstar`).
See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_light_assimilation_farquhar(
    J::FT,
    ci::FT,
    Γstar::FT,
    _...,
) where {FT}
    Aj = J * (ci - Γstar) / (4 * (ci + 2 * Γstar))
    return Aj
end

"""
    c4_light_assimilation_farquhar(::FT, ::FT, ::FT, APAR::FT, E::FT) where {FT}

Computes the electron transport limiting rate (`Aj`),
in units of moles CO2/m^2/s.

For C4 plants, this is a function of APAR and a efficiency parameter E, see Equation 11.70 of
 G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c4_light_assimilation_farquhar(
    ::FT,
    ::FT,
    ::FT,
    APAR::FT,
    E::FT,
) where {FT}
    Aj = APAR * E
    return Aj
end

"""
    gross_leaf_photosynthesis_at_a_point_Farquhar

Computes the gross photosynthesis at the leaf level, not accounting for moisture stress.
"""
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
    Jmax = max_electron_transport_farquhar(Vcmax25_leaf, ΔHJmax, T, To, R)
    J = electron_transport_farquhar(APAR_leaf_moles, Jmax, θj, ϕ)
    Vcmax_leaf = compute_Vcmax_farquhar(
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
    Γstar = co2_compensation_farquhar(Γstar25, ΔHΓstar, T, To, R)
    ci = intercellular_co2_farquhar(c_co2, Γstar, medlyn_factor)
    Aj = light_assimilation_farquhar(is_c3, J, ci, Γstar, APAR_leaf_moles, E)
    Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko(Ko25, ΔHko, T, To, R)
    Kmm = Kc * (1 + oi / Ko)
    Ac = rubisco_assimilation_farquhar(is_c3, Vcmax_leaf, ci, Γstar, Kmm)
    return gross_photosynthesis(Ac, Aj)
end

"""
    update_photosynthesis!(
        p,
        Y,
        model::FarquharModel,
        canopy
)

Computes the net leaf-level photosynthesis rate `An` (mol CO2/m^2/s) for the Farquhar
model, along with the dark leaf-level respiration `Rd` (mol CO2/m^2/s), and
canopy level gross photosynthesis (mol CO2/m^2/s).
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

    β = p.canopy.soil_moisture_stress.βm
    medlyn_factor = @. lazy(medlyn_term(g1, T_air, P_air, q_air, thermo_params))

    @. Rd = dark_respiration_farquhar(
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
    ) # has moisture stress
    A = @. lazy(
        gross_leaf_photosynthesis_at_a_point_Farquhar(
            T_canopy,
            compute_APAR_leaf_moles(
                f_abs,
                par_d,
                λ_γ_PAR,
                lightspeed,
                planck_h,
                N_a,
                LAI,
            ),
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
    ) # has moisture stress
    # Compute net assimilation:
    @. An = net_photosynthesis(A * β, Rd) # Rd has β accounted for already.
    # GPP
    @. GPP = A * β * LAI
end
Base.broadcastable(m::FarquharParameters) = tuple(m)

get_Vcmax25_leaf(p, m::FarquharModel) = m.parameters.Vcmax25
get_Rd_leaf(p, m::FarquharModel) = p.canopy.photosynthesis.Rd
get_An_leaf(p, m::FarquharModel) = p.canopy.photosynthesis.An

function get_J_over_Jmax(Y, p, canopy, m::FarquharModel)
    Jmax = compute_Jmax_leaf(Y, p, canopy, m) # lazy
    J = compute_J_leaf(Y, p, canopy, m) # lazy
    FT = eltype(canopy.parameters.earth_param_set)
    return @. lazy(J / max(Jmax, sqrt(eps(FT))))
end

function compute_Jmax_leaf(Y, p, canopy, m::FarquharModel) # used internally to farquhar; helper function
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    (; Vcmax25, ΔHJmax, To) = m.parameters
    R = LP.gas_constant(canopy.parameters.earth_param_set)
    return @. lazy(
        max_electron_transport_farquhar(Vcmax25, ΔHJmax, T_canopy, To, R),
    )
end

function compute_J_leaf(Y, p, canopy, m::FarquharModel) # used internally to farquhar; helper function
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)

    earth_param_set = canopy.parameters.earth_param_set
    f_abs_par = p.canopy.radiative_transfer.par.abs
    par_d = p.canopy.radiative_transfer.par_d
    (; λ_γ_PAR,) = canopy.radiative_transfer.parameters
    c = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    LAI = p.canopy.hydraulics.area_index.leaf
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

    Jmax_leaf = compute_Jmax_leaf(Y, p, canopy, m)
    (; θj, ϕ) = m.parameters
    return @. lazy(
        electron_transport_farquhar(APAR_leaf_moles, Jmax_leaf, θj, ϕ),
    )
end

## For interfacing with ClimaParams

"""
    function FarquharParameters(
        ::Type{FT},
        is_c3::Union{FT, ClimaCore.Fields.Field};
        Vcmax25 = FT(5e-5),
        kwargs...  # For individual parameter overrides
    )

    function FarquharParameters(
        toml_dict::CP.ParamDict,
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
    toml_dict::CP.ParamDict,
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

"""
    max_electron_transport_farquhar(
        Vcmax25::FT,
        ΔHJmax::FT,
        T::FT,
        To::FT,
        R::FT,
    ) where {FT}

Computes the maximum potential rate of electron transport (`Jmax`),
in units of mol/m^2/s,
as a function of Vcmax at 25 °C (`Vcmax25`),
a constant (`ΔHJmax`), a standard temperature (`To`),
the unversal gas constant (`R`), and the temperature (`T`).

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function max_electron_transport_farquhar(
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
    electron_transport_farquhar(APAR::FT,
                       Jmax::FT,
                       θj::FT,
                       ϕ::FT) where {FT}

Computes the rate of electron transport (`J`),
in units of mol/m^2/s, as a function of
the maximum potential rate of electron transport (`Jmax`),
absorbed photosynthetically active radiation (`APAR`),
an empirical "curvature parameter" (`θj`; Bonan Eqn 11.21)
and the quantum yield of photosystem II (`ϕ`).

In this function, all fluxes are per unit leaf area, including APAR.

See Ch 11, G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function electron_transport_farquhar(
    APAR::FT,
    Jmax::FT,
    θj::FT,
    ϕ::FT,
) where {FT}
    # Light utilization of APAR
    IPSII = ϕ * APAR / 2
    # This is a solution to a quadratic equation
    # θj *J^2 - (IPSII+Jmax)*J+IPSII*Jmax = 0, Equation 11.21
    J =
        (IPSII + Jmax - sqrt((IPSII + Jmax)^2 - 4 * θj * IPSII * Jmax)) /
        (2 * θj)
    return J
end

"""
    compute_Vcmax_farquhar(is_c3::AbstractFloat, args...)

Calls the correct Vcmax function based on `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_compute_Vcmax`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_compute_Vcmax`.
"""
function compute_Vcmax_farquhar(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_compute_Vcmax_farquhar(args...) :
    c4_compute_Vcmax_farquhar(args...)
end

"""
    c4_compute_Vcmax_farquhar(Vcmax25::FT, T::FT, R::FT, To::FT, ::FT, Q10::FT, s1::FT, s2::FT,
                     s3::FT, s4::FT) where {FT}

Computes the maximum rate of carboxylation of Rubisco (`Vcmax`),
in units of mol/m^2/s,
as a function of temperature (`T`), the universal
gas constant `R`, Vcmax25,  and other parameters.

For C4 photosynthesis, this uses Equation 11.73 from G. Bonan
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c4_compute_Vcmax_farquhar(
    Vcmax25::FT,
    T::FT,
    R::FT,
    To::FT,
    ::FT,
    Q10::FT,
    s1::FT,
    s2::FT,
    s3::FT,
    s4::FT,
) where {FT}
    Vcmax =
        Vcmax25 * Q10^((T - To) / 10) / (1 + exp(s1 * (T - s2))) /
        (1 + exp(s3 * (s4 - T)))
    return Vcmax
end

"""
    c3_compute_Vcmax_farquhar(Vcmax25::FT, T::FT, R::FT, To::FT, ΔHVcmax::FT) where {FT}

Computes the maximum rate of carboxylation of Rubisco (`Vcmax`),
in units of mol/m^2/s,
as a function of temperature (`T`), the universal
gas constant `R`, Vcmax25, and other parameters.

For C3 photosynthesis, this uses Table 11.5 from G. Bonan:
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_compute_Vcmax_farquhar(
    Vcmax25::FT,
    T::FT,
    R::FT,
    To::FT,
    ΔHVcmax::FT,
    args...,
) where {FT}
    Vcmax = Vcmax25 * arrhenius_function(T, To, R, ΔHVcmax)
    return Vcmax
end


"""
    dark_respiration_farquhar(is_c3::AbstractFloat, args...)

Calls the correct dark respiration function based on `is_c3`.

A `is_c3` value of 1.0 corresponds to C3 photosynthesis and calls
`c3_dark_respiration`, while 0.0 corresponds to C4 photsynthesis and calls
`c4_dark_respiration`.

TODO: Generalize to p-model.
"""
function dark_respiration_farquhar(is_c3::AbstractFloat, args...)
    is_c3 > 0.5 ? c3_dark_respiration_farquhar(args...) :
    c4_dark_respiration_farquhar(args...)
end
"""
    c4_dark_respiration_farquhar(VCmax25::FT,
                        β::FT,
                        T::FT,
                        R::FT
                        To::FT,
                        ::FT,
                        ::FT,
                        Q10::FT,
                        s5::FT,
                        s6::FT,
                        fC4::FT) where {FT}

Computes dark respiration (`Rd`),
in units of mol CO2/m^2/s, as a function of
 the moisture stress factor (`β`),
the unversal gas constant (`R`), the temperature (`T`),
Vcmax25, and
other parameters.

See Equation 11.73 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c4_dark_respiration_farquhar(
    Vcmax25::FT,
    β::FT,
    T::FT,
    R::FT,
    To::FT,
    ::FT,
    ::FT,
    Q10::FT,
    s5::FT,
    s6::FT,
    fC4::FT,
) where {FT}
    Rd = fC4 * Vcmax25 * β * Q10^((T - To) / 10) / (1 + exp(s5 * (T - s6)))
    return Rd
end

"""
    c3_dark_respiration_farquhar(Vcmax25::FT, β::FT,
                        T::FT,
                        R::FT,
                        To::FT,
                        fC3::FT,
                        ΔHRd::FT,) where {FT}

Computes dark respiration (`Rd`),
in units of mol CO2/m^2/s, as a function of
 the moisture stress factor (`β`),
the unversal gas constant (`R`), and the temperature (`T`),
Vcmax25,  and
other parameters.

See Table 11.5 of G. Bonan's textbook,
Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function c3_dark_respiration_farquhar(
    Vcmax25::FT,
    β::FT,
    T::FT,
    R::FT,
    To::FT,
    fC3::FT,
    ΔHRd::FT,
    _...,
) where {FT}
    Rd = fC3 * Vcmax25 * β * arrhenius_function(T, To, R, ΔHRd)
    return Rd
end


"""
    intercellular_co2_farquhar(ca::FT, Γstar::FT, medlyn_factor::FT) where{FT}

Computes the intercellular CO2 concentration (mol/mol) given the
atmospheric concentration (`ca`, mol/mol), the CO2 compensation (`Γstar`,
 mol/mol), and the Medlyn factor (unitless).
"""
function intercellular_co2_farquhar(
    ca::FT,
    Γstar::FT,
    medlyn_term::FT,
) where {FT}
    c_i = max(ca * (1 - 1 / medlyn_term), Γstar)
    return c_i
end

"""
    co2_compensation_farquhar(Γstar25::FT,
                     ΔHΓstar::FT,
                     T::FT,
                     To::FT,
                     R::FT) where {FT}

Computes the CO2 compensation point (`Γstar`),
in units of mol/mol,
as a function of its value at 25 °C (`Γstar25`),
a constant energy of activation (`ΔHΓstar`), a standard temperature (`To`),
the unversal gas constant (`R`), and the temperature (`T`).

See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
function co2_compensation_farquhar(
    Γstar25::FT,
    ΔHΓstar::FT,
    T::FT,
    To::FT,
    R::FT,
) where {FT}
    Γstar = Γstar25 * arrhenius_function(T, To, R, ΔHΓstar)
    return Γstar
end

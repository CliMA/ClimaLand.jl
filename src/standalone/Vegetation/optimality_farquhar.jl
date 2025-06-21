export OptimalityFarquharParameters, OptimalityFarquharModel


"""
    OptimalityFarquharParameters{FT<:AbstractFloat}

The required parameters for the optimality Farquhar photosynthesis model.
Currently, only C3 photosynthesis is supported.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OptimalityFarquharParameters{
    FT <: AbstractFloat,
    MECH <: Union{FT, ClimaCore.Fields.Field},
}
    "Photosynthesis mechanism: C3 only"
    is_c3::MECH
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
    "Intercellular O2 concentration (mol/mol); taken to be constant"
    oi::FT
    "Quantum yield of photosystem II (Bernacchi, 2003; unitless)"
    ϕ::FT
    "Curvature parameter, a fitting constant to compute J, unitless"
    θj::FT
    "Constant factor appearing the dark respiration term for C3 plants, equal to 0.015."
    fC3::FT
    "Fitting constant to compute the moisture stress factor (Pa^{-1})"
    sc::FT
    "Fitting constant to compute the moisture stress factor (Pa)"
    pc::FT
    "Constant describing cost of maintaining electron transport (unitless)"
    c::FT
end

Base.eltype(::OptimalityFarquharParameters{FT}) where {FT} = FT

"""
    OptimalityFarquharModel{FT,
                            OPFT <: OptimalityFarquharParameters{FT}
                            } <: AbstractPhotosynthesisModel{FT}

Optimality model of Smith et al. (2019) for estimating Vcmax, based on the assumption that Aj = Ac.
Smith et al. (2019). Global photosynthetic capacity is optimized to the environment. Ecology Letters, 22(3), 506–517. https://doi.org/10.1111/ele.13210
"""
struct OptimalityFarquharModel{FT, OPFT <: OptimalityFarquharParameters{FT}} <:
       AbstractPhotosynthesisModel{FT}
    "Required parameters for the Optimality based Farquhar model of Smith et al. (2019)"
    parameters::OPFT
end

function OptimalityFarquharModel{FT}(
    parameters::OptimalityFarquharParameters{FT},
) where {FT <: AbstractFloat}
    return OptimalityFarquharModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.auxiliary_vars(model::OptimalityFarquharModel) =
    (:An, :GPP, :Rd, :Vcmax25, :Jmax25)
ClimaLand.auxiliary_types(model::OptimalityFarquharModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::OptimalityFarquharModel) =
    (:surface, :surface, :surface, :surface, :surface)

"""
    update_photosynthesis!(p, Y, model::OptimalityFarquharModel,canopy)

Computes the net photosynthesis rate `An` (mol CO2/m^2/s)for the Optimality Farquhar model, along with the
dark respiration `Rd` (mol CO2/m^2/s), the value of `Vcmax25`(mol CO2/m^2/s) , and the gross primary
productivity `GPP` (mol CO2/m^2/s), and updates them in place.
"""
function update_photosynthesis!(p, Y, model::OptimalityFarquharModel, canopy)
    # unpack a bunch of stuff from p and params
    Rd = p.canopy.photosynthesis.Rd
    An = p.canopy.photosynthesis.An
    GPP = p.canopy.photosynthesis.GPP
    Vcmax25 = p.canopy.photosynthesis.Vcmax25
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
    energy_per_mole_photon_par = planck_h * lightspeed / λ_γ_PAR * N_a
    (; sc, pc) = canopy.photosynthesis.parameters
    (; g1,) = canopy.conductance.parameters
    n_stem = canopy.hydraulics.n_stem
    n_leaf = canopy.hydraulics.n_leaf
    i_end = n_stem + n_leaf
    par_d = p.canopy.radiative_transfer.par_d
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf

    (;
        Γstar25,
        ΔHVcmax,
        ΔHΓstar,
        ΔHRd,
        fC3,
        To,
        θj,
        ϕ,
        is_c3,
        oi,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        c,
    ) = model.parameters

    β = p.canopy.soil_moisture_stress.βm
    medlyn_factor = @. lazy(medlyn_term(g1, T_air, P_air, q_air, thermo_params))
    Γstar = @. lazy(co2_compensation(Γstar25, ΔHΓstar, T, To, R))
    ci = @. lazy(intercellular_co2(c_co2_air, Γstar, medlyn_factor))# may change?
    rates = @. lazy(
        optimality_max_photosynthetic_rates(
            f_abs * par_d / energy_per_mole_photon_par,
            θj,
            ϕ,
            oi,
            ci,
            Γstar,
            MM_Kc(Kc25, ΔHkc, T_canopy, To, R),
            MM_Ko(Ko25, ΔHko, T_canopy, To, R),
            c,
        ),
    )
    Jmax = rates.:1
    Vcmax = rates.:2
    J = @. lazy(
        electron_transport(
            f_abs * par_d / energy_per_mole_photon_par,
            Jmax,
            θj,
            ϕ,
        ),
    )
    Aj = @. lazy(light_assimilation(is_c3, J, ci, Γstar))
    Ac = @. lazy(
        rubisco_assimilation(
            is_c3,
            Vcmax,
            ci,
            Γstar,
            MM_Kc(Kc25, ΔHkc, T_canopy, To, R),
            MM_Ko(Ko25, ΔHko, T_canopy, To, R),
            oi,
        ),
    )

    @. Vcmax25 = Vcmax / arrhenius_function(T_canopy, To, R, ΔHVcmax)
    @. Jmax25 = Jmax / arrhenius_function(T_canopy, To, R, ΔHJmax)
    @. Rd = dark_respiration(is_c3, Vcmax25, β, T_canopy, R, To, fC3, ΔHRd)
    @. An = net_photosynthesis(Ac, Aj, Rd, β)
    # Compute GPP: TODO - move to diagnostics only
    @. GPP = compute_GPP(An, extinction_coeff(G_Function, cosθs), LAI, Ω)
end

get_Vcmax25(p, m::OptimalityFarquharModel) = p.canopy.photosynthesis.Vcmax25

function get_Jmax(Y, p, canopy, m::OptimalityFarquharModel)
    (; To, ΔHJmax) = m.parameters
    R = LP.get_default_parameter(eltype(m.parameters), :gas_constant)
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    Jmax25 = p.canopy.photosynthesis.Jmax25
    return @. lazy(Jmax25 * arrhenius_function(T_canopy, To, R, ΔHJmax))
end

function get_electron_transport(Y, p, canopy, m::OptimalityFarquharModel)
    (; θj, ϕ) = m.parameters
    earth_param_set = canopy.parameters.earth_param_set
    f_abs_par = p.canopy.radiative_transfer.par.abs
    par_d = p.canopy.radiative_transfer.par_d
    (; λ_γ_PAR,) = canopy.radiative_transfer.parameters
    c = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    APAR = @. lazy(compute_APAR(f_abs_par, par_d, λ_γ_PAR, c, planck_h, N_a))
    Jmax = get_Jmax(Y, p, canopy, m)
    return @. lazy(electron_transport(APAR, Jmax, θj, ϕ))
end

Base.broadcastable(m::OptimalityFarquharParameters) = tuple(m)

## For interfacing with ClimaParams

"""
    function OptimalityFarquharParameters(
        ::Type{FT},
        kwargs...  # For individual parameter overrides
    )

    function OptimalityFarquharParameters(
        toml_dict::CP.AbstractTOMLDict,
        kwargs...  # For individual parameter overrides
    )

Constructors for the OptimalityFarquharParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values.Possible calls:
```julia
ClimaLand.Canopy.OptimalityFarquharParameters(Float64)
# Kwarg overrides
ClimaLand.Canopy.OptimalityFarquharParameters(Float64; pc = 444444444)
# Toml Dictionary:
import ClimaParams as CP
toml_dict = CP.create_toml_dict(Float32);
ClimaLand.Canopy.OptimalityFarquharParameters(toml_dict; pc = 444444444)
```
"""
OptimalityFarquharParameters(
    ::Type{FT};
    kwargs...,
) where {FT <: AbstractFloat} =
    OptimalityFarquharParameters(CP.create_toml_dict(FT); kwargs...)

function OptimalityFarquharParameters(toml_dict; kwargs...)
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
        :electron_transport_maintenance => :c,
        :Γstar_activation_energy => :ΔHΓstar,
        :CO2_activation_energy => :ΔHkc,
        :moisture_stress_ref_water_pressure => :pc,
    )

    params = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    is_c3 = FT(1)
    return OptimalityFarquharParameters{FT, FT}(; params..., kwargs..., is_c3)
end

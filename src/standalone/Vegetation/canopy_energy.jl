export PrescribedCanopyTempModel,
    BigLeafEnergyModel,
    AbstractCanopyEnergyModel,
    canopy_temperature,
    root_energy_flux_per_ground_area!,
    BigLeafEnergyParameters,
    total_energy_per_area!

"""
    AbstractCanopyEnergyModel{FT}

An abstract struct for the Canopy Energy Models. Both PrescribedCanopyTempModel and
BigLeafEnergyModel are subtypes of this abstract type.
"""
abstract type AbstractCanopyEnergyModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(model::AbstractCanopyEnergyModel) = :energy


ClimaLand.auxiliary_vars(model::AbstractCanopyEnergyModel) =
    (:fa_energy_roots, :∂LW_n∂Tc, :∂qc∂Tc)
ClimaLand.auxiliary_types(model::AbstractCanopyEnergyModel{FT}) where {FT} =
    (FT, FT, FT)
ClimaLand.auxiliary_domain_names(model::AbstractCanopyEnergyModel) =
    (:surface, :surface, :surface)

"""
    PrescribedCanopyTempModel{FT} <: AbstractCanopyEnergyModel{FT}

A model for the energy of the canopy which assumes the canopy temperature
is the same as the atmosphere temperature prescribed in the
`PrescribedAtmos` struct.

No equation for the energy of the canopy is solved.
"""
struct PrescribedCanopyTempModel{FT} <: AbstractCanopyEnergyModel{FT} end

"""
    canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p)

Returns the canopy temperature under the `PrescribedCanopyTemp` model,
where the canopy temperature is assumed to be the same as the atmosphere
temperature.
"""
function canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p)
    p.drivers.T
end

## Prognostic Canopy Temperature
"""
    BigLeafEnergyParameters{FT <: AbstractFloat}

$(DocStringExtensions.FIELDS)

"""
struct BigLeafEnergyParameters{FT <: AbstractFloat}
    "Specific heat per emitting area [J/m^2/K]"
    ac_canopy::FT
end

Base.eltype(::BigLeafEnergyParameters{FT}) where {FT} = FT

"""
    BigLeafEnergyParameters{FT}(; ac_canopy = FT(2e3)) where {FT}

Construct `BigLeafEnergyParameters` with the default `ac_canopy = 2e3`.
"""
function BigLeafEnergyParameters{FT}(; ac_canopy = FT(2e3)) where {FT}
    return BigLeafEnergyParameters{FT}(ac_canopy)
end

"""
    BigLeafEnergyParameters(toml_dict::CP.ParamDict;
                            ac_canopy = toml_dict["ac_canopy"],
                        )

Construct `BigLeafEnergyParameters` from a TOML dict.
"""
function BigLeafEnergyParameters(
    toml_dict::CP.ParamDict;
    ac_canopy = toml_dict["ac_canopy"],
)
    FT = CP.float_type(toml_dict)
    return BigLeafEnergyParameters{FT}(ac_canopy)
end


"""
    BigLeafEnergyModel{FT} <: AbstractCanopyEnergyModel{FT}
"""
struct BigLeafEnergyModel{FT, BEP <: BigLeafEnergyParameters{FT}} <:
       AbstractCanopyEnergyModel{FT}
    parameters::BEP
end

function BigLeafEnergyModel{FT}(
    parameters::BigLeafEnergyParameters{FT},
) where {FT <: AbstractFloat}
    return BigLeafEnergyModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.prognostic_vars(model::BigLeafEnergyModel) = (:T,)
ClimaLand.prognostic_types(model::BigLeafEnergyModel{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(model::BigLeafEnergyModel) = (:surface,)

"""
    canopy_temperature(model::BigLeafEnergyModel, canopy, Y, p)

Returns the canopy temperature under the `BigLeafEnergyModel` model,
where the canopy temperature is modeled prognostically.
"""
canopy_temperature(model::BigLeafEnergyModel, canopy, Y, p) = Y.canopy.energy.T

function make_compute_imp_tendency(
    model::BigLeafEnergyModel{FT},
    canopy,
) where {FT}
    function compute_imp_tendency!(dY, Y, p, t)
        area_index = p.canopy.hydraulics.area_index
        ac_canopy = model.parameters.ac_canopy
        # Energy Equation:
        # (ρc_canopy h_canopy AI) ∂T∂t = -∑F
        # or( ac_canopy AI)∂T∂t = -∑F
        # where ∑F = F_sfc - F_bot, and both F_sfc and F_bot are per
        # unit area ground [W/m^2].
        # Because they are per unit area ground, we need the factor of
        # area index on the LHF, as ac_canopy [J/m^2/K]
        # is per unit area plant.

        # d(energy.T) = - net_energy_flux / specific_heat_capacity
        # To prevent dividing by zero, change AI" to
        # "max(AI, eps(FT))"

        @. dY.canopy.energy.T =
            -(
                -p.canopy.radiative_transfer.LW_n -
                p.canopy.radiative_transfer.SW_n +
                p.canopy.turbulent_fluxes.shf +
                p.canopy.turbulent_fluxes.lhf - p.canopy.energy.fa_energy_roots
            ) / (ac_canopy * max(area_index.leaf + area_index.stem, eps(FT)))
    end
    return compute_imp_tendency!
end

"""
    root_energy_flux_per_ground_area!(
        fa_energy::ClimaCore.Fields.Field,
        ground::PrescribedGroundConditions{FT},
        model::AbstractCanopyEnergyModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    ) where {FT}


A method which updates the ClimaCore.Fields.Field `fa_energy` in place
with  the energy flux associated with the root-soil
water flux for the `CanopyModel` run in standalone mode,
with a `PrescribedGroundConditions`.This value is ignored and set to zero
in this case.

Background information: This energy
flux is not typically included in land surface
models. We account for it when the soil model is prognostic because
the soil model includes the energy in the soil water in its energy
balance; therefore, in order to conserve energy, the canopy model
must account for it as well.
"""
function root_energy_flux_per_ground_area!(
    fa_energy::ClimaCore.Fields.Field,
    ground::PrescribedGroundConditions{FT},
    model::AbstractCanopyEnergyModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
) where {FT}
    fa_energy .= FT(0)
end

function ClimaLand.make_compute_jacobian(
    model::BigLeafEnergyModel{FT},
    canopy,
) where {FT}
    function compute_jacobian!(
        jacobian::MatrixFields.FieldMatrixWithSolver,
        Y,
        p,
        dtγ,
        t,
    )
        (; matrix) = jacobian

        # The derivative of the residual with respect to the prognostic variable
        ∂Tres∂T = matrix[@name(canopy.energy.T), @name(canopy.energy.T)]
        ∂LHF∂qc = p.canopy.turbulent_fluxes.∂LHF∂qc
        ∂SHF∂Tc = p.canopy.turbulent_fluxes.∂SHF∂Tc
        ∂LW_n∂Tc = p.canopy.energy.∂LW_n∂Tc
        ∂qc∂Tc = p.canopy.energy.∂qc∂Tc
        ϵ_c = p.canopy.radiative_transfer.ϵ
        area_index = p.canopy.hydraulics.area_index
        ac_canopy = model.parameters.ac_canopy
        earth_param_set = canopy.parameters.earth_param_set
        _T_freeze = LP.T_freeze(earth_param_set)
        _σ = LP.Stefan(earth_param_set)
        @. ∂LW_n∂Tc = -2 * 4 * _σ * ϵ_c * Y.canopy.energy.T^3 # ≈ ϵ_ground = 1
        @. ∂qc∂Tc = partial_q_sat_partial_T_liq(
            p.drivers.P,
            Y.canopy.energy.T - _T_freeze,
        )# use atmos air pressure as approximation for surface air pressure
        @. ∂Tres∂T =
            float(dtγ) * MatrixFields.DiagonalMatrixRow(
                (∂LW_n∂Tc - ∂SHF∂Tc - ∂LHF∂qc * ∂qc∂Tc) /
                (ac_canopy * max(area_index.leaf + area_index.stem, eps(FT))),
            ) - (I,)
    end
    return compute_jacobian!
end

"""
    partial_q_sat_partial_T_liq(P::FT, T::FT) where {FT}

Computes the quantity ∂q_sat∂T at temperature T and pressure P,
over liquid water. The temperature must be in Celsius.

Uses the polynomial approximation from Flatau et al. (1992).
"""
function partial_q_sat_partial_T_liq(P::FT, T::FT) where {FT}
    esat = FT(
        6.11213476e2 +
        4.44007856e1 * T +
        1.43064234 * T^2 +
        2.64461437e-2 * T^3 +
        3.05903558e-4 * T^4 +
        1.96237241e-6 * T^5 +
        8.92344772e-9 * T^6 - 3.73208410e-11 * T^7 + 2.09339997e-14 * T^8,
    )
    desatdT = FT(
        4.44017302e1 +
        2.86064092 * T +
        7.94683137e-2 * T^2 +
        1.21211669e-3 * T^3 +
        1.03354611e-5 * T^4 +
        4.04125005e-8 * T^5 - 7.88037859e-11 * T^6 - 1.14596802e-12 * T^7 +
        3.81294516e-15 * T^8,
    )

    return FT(0.622) * P / (P - FT(0.378) * esat)^2 * desatdT
end

"""
    ClimaLand.total_energy_per_area!(
        surface_field,
        model::BigLeafEnergyModel,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value of
the big leaf model's energy.

Note that this assumes that there is at most a single canopy and stem
component, and that the area index for them refers to the integrated
area index (in height) - not the value per layer.
"""
function ClimaLand.total_energy_per_area!(
    surface_field,
    model::BigLeafEnergyModel,
    Y,
    p,
    t,
)
    area_index = p.canopy.hydraulics.area_index
    @. surface_field .=
        model.parameters.ac_canopy *
        (area_index.stem + area_index.leaf) *
        Y.canopy.energy.T
    return nothing
end

"""
    ClimaLand.total_energy_per_area!(
        surface_field,
        model::AbstractCanopyEnergyModel,
        Y,
        p,
        t,
)

A default function which errors for generic energy models for the
canopy.

Note that we only support two models for canopy energy. The `BigLeafEnergyModel` has a special method for this, and the other has the temperature
prescribed and does not conserve energy.
"""
function ClimaLand.total_energy_per_area!(
    surface_field,
    model::AbstractCanopyEnergyModel,
    Y,
    p,
    t,
)
    @error("This method has not been implemented.")
end

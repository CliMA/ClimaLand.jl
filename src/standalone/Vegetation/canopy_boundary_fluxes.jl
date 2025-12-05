using SurfaceFluxes
using Thermodynamics
using StaticArrays
import SurfaceFluxes.Parameters as SFP

import ClimaLand: turbulent_fluxes!, AbstractBC, get_earth_param_set

function get_earth_param_set(model::CanopyModel)
    return model.earth_param_set
end


"""
    AbstractCanopyBC <: ClimaLand.AbstractBC

An abstract type for boundary conditions for the canopy model.
"""
abstract type AbstractCanopyBC <: ClimaLand.AbstractBC end
"""
    AtmosDrivenCanopyBC{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
        G <: AbstractGroundConditions,
        R <: AbstractCanopyFluxParameterization,
        C::Tuple
    } <: AbstractCanopyBC

A struct used to specify the canopy fluxes, referred
to as "boundary conditions", at the surface and
bottom of the canopy, for water and energy.

These fluxes include turbulent surface fluxes
computed with Monin-Obukhov theory, radiative fluxes,
and root extraction.
$(DocStringExtensions.FIELDS)
"""
struct AtmosDrivenCanopyBC{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
    G <: AbstractGroundConditions,
    R <: AbstractCanopyFluxParameterization,
    C <: Tuple,
} <: AbstractCanopyBC
    "The atmospheric conditions driving the model"
    atmos::A
    "The radiative fluxes driving the model"
    radiation::B
    "Ground conditions"
    ground::G
    "Turbulent flux (latent, sensible, vapor, and momentum) parameterization"
    turbulent_flux_parameterization::R
    "Prognostic land components present"
    prognostic_land_components::C
end

"""
    AtmosDrivenCanopyBC(
        atmos,
        radiation,
        ground,
        turbulent_flux_parameterization;
        prognostic_land_components = (:canopy,),
    )

An outer constructor for `AtmosDrivenCanopyBC` which is
intended for use as a default when running canopy
models.

This is also checks the logic that:
- If the `ground` field is Prescribed, :soil should not be a prognostic_land_component
- If the `ground` field is not Prescribed, :soil should be modeled prognostically.
"""
function AtmosDrivenCanopyBC(
    atmos,
    radiation,
    ground,
    turbulent_flux_parameterization;
    prognostic_land_components = (:canopy,),
)
    if typeof(ground) <: PrescribedGroundConditions
        @assert !(:soil ∈ prognostic_land_components)
    else
        @assert :soil ∈ prognostic_land_components
    end

    args = (
        atmos,
        radiation,
        ground,
        turbulent_flux_parameterization,
        prognostic_land_components,
    )
    return AtmosDrivenCanopyBC(args...)
end

function ClimaLand.get_drivers(bc::AtmosDrivenCanopyBC)
    if typeof(bc.ground) <: PrescribedGroundConditions
        return (bc.atmos, bc.radiation, bc.ground)
    else
        return (bc.atmos, bc.radiation)
    end
end


function make_update_boundary_fluxes(canopy::CanopyModel)
    function update_boundary_fluxes!(p, Y, t)
        canopy_boundary_fluxes!(p, canopy, canopy.energy, Y, t)
    end
    return update_boundary_fluxes!
end

"""
    canopy_boundary_fluxes!(p::NamedTuple,
                            canopy::CanopyModel,
                            Y::ClimaCore.Fields.FieldVector,
                            t,
                            )

Computes the boundary fluxes for the canopy prognostic
equations; updates the specific fields in the auxiliary
state `p` which hold these variables. This function is called
within the explicit tendency of the canopy model.

- `p.canopy.turbulent_fluxes`: Canopy SHF, LHF, transpiration, derivatives of these with respect to T,q
- `p.canopy.hydraulics.fa[end]`: Transpiration
- `p.canopy.hydraulics.fa_roots`: Root water flux
- `p.canopy.radiative_transfer.LW_n`: net long wave radiation
- `p.canopy.radiative_transfer.SW_n`: net short wave radiation
"""
function canopy_boundary_fluxes!(
    p::NamedTuple,
    canopy::CanopyModel,
    energy::AbstractCanopyEnergyModel,
    Y::ClimaCore.Fields.FieldVector,
    t,
)
    bc = canopy.boundary_conditions
    radiation = bc.radiation
    atmos = bc.atmos
    sf_parameterization = bc.turbulent_flux_parameterization
    root_water_flux = p.canopy.hydraulics.fa_roots
    root_energy_flux = p.canopy.energy.fa_energy_roots
    fa = p.canopy.hydraulics.fa
    LAI = p.canopy.biomass.area_index.leaf
    SAI = p.canopy.biomass.area_index.stem
    canopy_tf = p.canopy.turbulent_fluxes
    i_end = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf
    # Compute transpiration, SHF, LHF
    ClimaLand.turbulent_fluxes!(
        canopy_tf,
        atmos,
        sf_parameterization,
        canopy,
        Y,
        p,
        t,
    )
    # Transpiration is per unit ground area, not leaf area (mult by LAI)
    fa.:($i_end) .= PlantHydraulics.transpiration_per_ground_area(
        canopy.hydraulics.transpiration,
        Y,
        p,
        t,
    )
    # Note that in the three functions below,
    # we dispatch off of the ground conditions `bc.ground`
    # to handle standalone canopy simulations vs integrated ones

    # Update the root flux of water per unit ground area in place
    root_water_flux_per_ground_area!(
        root_water_flux,
        bc.ground,
        canopy.hydraulics,
        canopy,
        Y,
        p,
        t,
    )
    # Update the root flux of energy per unit ground area in place
    root_energy_flux_per_ground_area!(
        root_energy_flux,
        bc.ground,
        canopy.energy,
        canopy,
        Y,
        p,
        t,
    )

    # Update the canopy radiation
    canopy_sw_energy_fluxes!(
        p,
        bc.ground,
        canopy,
        bc.radiation,
        canopy.earth_param_set,
        Y,
        t,
    )

    canopy_lw_energy_fluxes!(
        p,
        bc.ground,
        canopy,
        bc.radiation,
        canopy.earth_param_set,
        Y,
        t,
    )

end


function flux_balance_equation(T::FT, ϵ_c::FT, ϵ_g::FT, T_g::FT, σ::FT, LW_d::FT, SW_n::FT, root_energy_flux::FT, ts_in, u::FT, h::FT, r_stomata_canopy::FT, d_sfc::FT, z_0m::FT, z_0b::FT, Cd::FT, LAI::FT, SAI::FT, earth_param_set, gustiness::FT) where {FT}
    T = max(T, ts_in.T- FT(5))
    T = min(T, ts_in.T+FT(5))
    LW_n = canopy_net_longwave(T, ϵ_c, ϵ_g, T_g, σ, LW_d)
    R_n = SW_n + LW_n
    turb_fluxes = canopy_compute_turbulent_fluxes_at_a_point(ts_in, u, h, gustiness,T, r_stomata_canopy, d_sfc, z_0m, z_0b, Cd, LAI, SAI, earth_param_set)
    return FT(-R_n + (turb_fluxes.shf + turb_fluxes.lhf) - root_energy_flux)
end
function root_solve(ϵ_c, ϵ_g, T_g, σ, LW_d, SW_n, root_energy_flux, ts_in, u, h, r_stomata_canopy, d_sfc, z_0m, z_0b, Cd, LAI, SAI, earth_param_set, gustiness)
    f(T) = flux_balance_equation(T, ϵ_c, ϵ_g, T_g, σ, LW_d, SW_n, root_energy_flux, ts_in, u, h, r_stomata_canopy, d_sfc, z_0m, z_0b, Cd, LAI, SAI, earth_param_set, gustiness)
    FT = typeof(ϵ_c)
    soln = RootSolvers.find_zero(f, SecantMethod(ts_in.T - FT(3), ts_in.T +FT(3)))
    return soln.root
end
                   
function solve_for_T!(p, Y, bc, atmos, radiation, earth_param_set)
    T = p.canopy.energy.T
    root_energy_flux = p.canopy.energy.fa_energy_roots
    LW_d = p.drivers.LW_d
    SW_n = p.canopy.radiative_transfer.SW_n
    ts_in = p.drivers.thermal_state
    u = p.drivers.u
    h = atmos.h
    gustiness = atmos.gustiness
    r_stomata_canopy = p.canopy.conductance.r_stomata_canopy
    sf_parameterization = bc.turbulent_flux_parameterization
    d_sfc = sf_parameterization.displ
    z_0m = sf_parameterization.z_0m
    z_0b = sf_parameterization.z_0b
    Cd = sf_parameterization.Cd
    LAI = p.canopy.biomass.area_index.leaf
    SAI = p.canopy.biomass.area_index.stem
    ϵ_c = p.canopy.radiative_transfer.ϵ # this takes into account LAI/SAI
    # Long wave: use ground conditions from the ground driver
    T_g = p.drivers.T_ground
    ϵ_g = bc.ground.ϵ
    FT = eltype(T)
    _σ = FT(LP.Stefan(earth_param_set))
    T .= root_solve.(ϵ_c, ϵ_g, T_g, _σ, LW_d, SW_n, root_energy_flux, ts_in, u, h, r_stomata_canopy, d_sfc, z_0m, z_0b, Cd, LAI, SAI, earth_param_set, gustiness)
    if parent(T)[1] < 0
        @show ϵ_c, ϵ_g, T_g, _σ, LW_d, SW_n, root_energy_flux, ts_in, u, h, r_stomata_canopy, d_sfc, z_0m, z_0b, Cd, LAI, SAI, earth_param_set, gustiness
    end
end

function canopy_boundary_fluxes!(
    p::NamedTuple,
    canopy::CanopyModel,
    energy::SteadyStateModel,
    Y::ClimaCore.Fields.FieldVector,
    t,
)
    bc = canopy.boundary_conditions
    radiation = bc.radiation
    atmos = bc.atmos
    sf_parameterization = bc.turbulent_flux_parameterization
    root_water_flux = p.canopy.hydraulics.fa_roots
    root_energy_flux = p.canopy.energy.fa_energy_roots
    fa = p.canopy.hydraulics.fa
    LAI = p.canopy.biomass.area_index.leaf
    SAI = p.canopy.biomass.area_index.stem
    canopy_tf = p.canopy.turbulent_fluxes
    i_end = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf

    #Explicit terms
    root_water_flux_per_ground_area!(
        root_water_flux,
        bc.ground,
        canopy.hydraulics,
        canopy,
        Y,
        p,
        t,
    )
    # Update the root flux of energy per unit ground area in place
    root_energy_flux_per_ground_area!(
        root_energy_flux,
        bc.ground,
        canopy.energy,
        canopy,
        Y,
        p,
        t,
    )
    canopy_sw_energy_fluxes!(
        p,
        bc.ground,
        canopy,
        bc.radiation,
        canopy.earth_param_set,
        Y,
        t,
    )

    # Solve for T
    solve_for_T!(p, Y, bc, atmos, radiation, canopy.earth_param_set)
    # compute fluxes using this T

    
    # Compute transpiration, SHF, LHF
    ClimaLand.turbulent_fluxes!(
        canopy_tf,
        atmos,
        sf_parameterization,
        canopy,
        Y,
        p,
        t,
    )
    # Transpiration is per unit ground area, not leaf area (mult by LAI)
    fa.:($i_end) .= PlantHydraulics.transpiration_per_ground_area(
        canopy.hydraulics.transpiration,
        Y,
        p,
        t,
    )

    canopy_lw_energy_fluxes!(
        p,
        bc.ground,
        canopy,
        bc.radiation,
        canopy.earth_param_set,
        Y,
        t,
    )

end

"""
    function ClimaLand.turbulent_fluxes!(
        dest,
        atmos::AbstractAtmosphericDrivers,
        sf_parameterization::MoninObukhovCanopyFluxes,
        model::CanopyModel,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

A canopy specific function for compute turbulent fluxes with the atmosphere;
returns the latent heat flux, sensible heat flux, vapor flux, and aerodynamic resistance.

We cannot use the default version in src/shared_utilities/drivers.jl
because the canopy requires a different resistance for vapor and sensible heat
fluxes, and the resistances depend on ustar, which we must compute using
SurfaceFluxes before adjusting to account for these resistances.
"""
function ClimaLand.turbulent_fluxes!(
    dest,
    atmos::AbstractAtmosphericDrivers,
    sf_parameterization::MoninObukhovCanopyFluxes,
    model::CanopyModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    @static if pkgversion(ClimaLand) < v"1.1" &&
               typeof(atmos) <: CoupledAtmosphere
        return nothing # Coupler has already computed the fluxes and updated `dest` in place
    end
    T_sfc = canopy_temperature(model.energy, model, Y, p)
    r_stomata_canopy = p.canopy.conductance.r_stomata_canopy
    d_sfc = sf_parameterization.displ
    z_0m = sf_parameterization.z_0m
    z_0b = sf_parameterization.z_0b
    Cd = sf_parameterization.Cd

    momentum_fluxes = Val(return_momentum_fluxes(atmos))
    dest .=
        canopy_turbulent_fluxes_at_a_point.(
            momentum_fluxes, # return_extra_fluxes
            p.drivers.thermal_state,
            p.drivers.u,
            atmos.h,
            atmos.gustiness,
            T_sfc,
            r_stomata_canopy,
            d_sfc,
            z_0m,
            z_0b,
            Cd,
            p.canopy.biomass.area_index.leaf,
            p.canopy.biomass.area_index.stem,
            Ref(model.earth_param_set),
        )
    return nothing
end


"""
    canopy_turbulent_fluxes_at_a_point(return_extra_fluxes, args...)

This is a wrapper function that allows us to dispatch on the type of `return_extra_fluxes`
as we compute the canopy turbulent fluxes pointwise. This is needed because space for the
extra fluxes is only allocated in the cache when running with a `CoupledAtmosphere`.
The function `canopy_compute_turbulent_fluxes_at_a_point` does the actual flux computation.

The `return_extra_fluxes` argument indicates whether to return the following:
- momentum fluxes (`ρτxz`, `ρτyz`)
- buoyancy flux (`buoy_flux`)
"""
function canopy_turbulent_fluxes_at_a_point(
    return_extra_fluxes::Val{false},
    args...,
)
    (LH, SH, Ẽ, r_ae, ∂LHF∂qc, ∂SHF∂Tc, _, _, _) =
        canopy_compute_turbulent_fluxes_at_a_point(args...)
    return (
        lhf = LH,
        shf = SH,
        transpiration = Ẽ,
        r_ae = r_ae,
        ∂LHF∂qc = ∂LHF∂qc,
        ∂SHF∂Tc = ∂SHF∂Tc,
    )
end
function canopy_turbulent_fluxes_at_a_point(
    return_extra_fluxes::Val{true},
    args...,
)
    (LH, SH, Ẽ, r_ae, ∂LHF∂qc, ∂SHF∂Tc, ρτxz, ρτyz, buoy_flux) =
        canopy_compute_turbulent_fluxes_at_a_point(args...)
    return (
        lhf = LH,
        shf = SH,
        transpiration = Ẽ,
        r_ae = r_ae,
        ∂LHF∂qc = ∂LHF∂qc,
        ∂SHF∂Tc = ∂SHF∂Tc,
        ρτxz = ρτxz,
        ρτyz = ρτyz,
        buoy_flux = buoy_flux,
    )
end


"""
    function canopy_compute_turbulent_fluxes_at_a_point(
        ts_in,
        u::Union{FT, SVector{2, FT}},
        h::FT,    
        gustiness::FT,
        T_sfc::FT,
        r_stomata_canopy::FT,
        d_sfc::FT,
        z_0m::FT,
        z_0b::FT,
        Cd::FT,
        LAI::FT,
        SAI::FT,
        earth_param_set::EP;
    ) where {FT <: AbstractFloat, EP, F}

Computes the turbulent surface fluxes for the canopy at a point
and returns the fluxes in a named tuple.

Note that an additiontal resistance is used in computing both
evaporation and sensible heat flux, and this modifies the output
of `SurfaceFluxes.surface_conditions`.
"""
function canopy_compute_turbulent_fluxes_at_a_point(
    ts_in,
    u::Union{FT, SVector{2, FT}},
    h::FT,
    gustiness::FT,
    T_sfc::FT,
    r_stomata_canopy::FT,
    d_sfc::FT,
    z_0m::FT,
    z_0b::FT,
    Cd::FT,
    LAI::FT,
    SAI::FT,
    earth_param_set::EP;
) where {FT <: AbstractFloat, EP}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    if u isa FT
        u = SVector{2, FT}(u, 0)
    end
    state_in = SurfaceFluxes.StateValues(h - d_sfc, u, ts_in)

    ρ_sfc = ClimaLand.compute_ρ_sfc(thermo_params, ts_in, T_sfc)
    q_sfc = Thermodynamics.q_vap_saturation_generic(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Liquid(),
    )
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)
    state_sfc = SurfaceFluxes.StateValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)

    # State containers
    states = SurfaceFluxes.ValuesOnly(
        state_in,
        state_sfc,
        z_0m,
        z_0b,
        beta = FT(1),
        gustiness = gustiness,
    )
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    scheme = SurfaceFluxes.PointValueScheme()
    conditions =
        SurfaceFluxes.surface_conditions(surface_flux_params, states, scheme)
    _LH_v0::FT = LP.LH_v0(earth_param_set)
    _ρ_liq::FT = LP.ρ_cloud_liq(earth_param_set)
    cp_m_sfc::FT = Thermodynamics.cp_m(thermo_params, ts_sfc)
    r_ae::FT = 1 / (conditions.Ch * SurfaceFluxes.windspeed(states))
    ustar::FT = conditions.ustar
    ρ_sfc = Thermodynamics.air_density(thermo_params, ts_sfc)
    T_int = Thermodynamics.air_temperature(thermo_params, ts_in)
    Rm_int = Thermodynamics.gas_constant_air(thermo_params, ts_in)
    ρ_air = Thermodynamics.air_density(thermo_params, ts_in)
    r_b_leaf::FT = 1 / (Cd * max(ustar, gustiness))
    r_b_canopy_lai = r_b_leaf / LAI
    r_b_canopy_total = r_b_leaf / (LAI + SAI)

    E0::FT =
        SurfaceFluxes.evaporation(surface_flux_params, states, conditions.Ch)
    E = E0 * r_ae / (r_b_canopy_lai + r_stomata_canopy + r_ae) # CLM 5, tech note Equation 5.101, and fig 5.2b, assuming all sunlit, f_wet = 0
    Ẽ = E / _ρ_liq

    SH =
        SurfaceFluxes.sensible_heat_flux(
            surface_flux_params,
            conditions.Ch,
            states,
            scheme,
        ) * r_ae / (r_b_canopy_total + r_ae)
    # The above follows from CLM 5, tech note Equation 5.88, setting H_v = SH and solving to remove T_s, ignoring difference between cp in atmos and above canopy.
    LH = _LH_v0 * E

    # Derivatives
    # We ignore ∂r_ae/∂T_sfc, ∂u*/∂T_sfc, ∂r_stomata∂Tc
    ∂ρsfc∂Tc =
        ρ_air *
        (Thermodynamics.cv_m(thermo_params, ts_in) / Rm_int) *
        (T_sfc / T_int)^(Thermodynamics.cv_m(thermo_params, ts_in) / Rm_int - 1) /
        T_int
    ∂cp_m_sfc∂Tc = 0 # Possibly can address at a later date

    ∂LHF∂qc =
        ρ_sfc * _LH_v0 / (r_b_canopy_lai + r_stomata_canopy + r_ae) +
        LH / ρ_sfc * ∂ρsfc∂Tc

    ∂SHF∂Tc =
        ρ_sfc * cp_m_sfc / (r_b_canopy_total + r_ae) +
        SH / ρ_sfc * ∂ρsfc∂Tc +
        SH / cp_m_sfc * ∂cp_m_sfc∂Tc

    return (
        lhf = LH,
        shf = SH,
        transpiration = Ẽ,
        r_ae = r_ae,
        ∂LHF∂qc = ∂LHF∂qc,
        ∂SHF∂Tc = ∂SHF∂Tc,
        ρτxz = conditions.ρτxz,
        ρτyz = conditions.ρτyz,
        conditions.buoy_flux,
    )
end

"""
    boundary_vars(bc, ::ClimaLand.TopBoundary)
    boundary_var_domain_names(bc, ::ClimaLand.TopBoundary)
    boundary_var_types(::AbstractCanopyEnergyModel, bc, ::ClimaLand.TopBoundary)

Fallbacks for the boundary conditions methods which add the turbulent
fluxes to the auxiliary variables.
"""
boundary_vars(bc, ::ClimaLand.TopBoundary) = (:turbulent_fluxes,)
boundary_var_domain_names(bc, ::ClimaLand.TopBoundary) = (:surface,)
boundary_var_types(::CanopyModel{FT}, bc, ::ClimaLand.TopBoundary) where {FT} =
    (
        NamedTuple{
            (:lhf, :shf, :transpiration, :r_ae, :∂LHF∂qc, :∂SHF∂Tc),
            Tuple{FT, FT, FT, FT, FT, FT},
        },
    )

"""
    boundary_var_types(
        ::CanopyModel{FT},
        ::AtmosDrivenCanopyBC{<:CoupledAtmosphere, <:CoupledRadiativeFluxes},
        ::ClimaLand.TopBoundary,
    ) where {FT}

An extension of the `boundary_var_types` method for AtmosDrivenCanopyBC. This
specifies the type of the additional variables.

This method includes additional flux-related properties needed by the atmosphere:
momentum fluxes (`ρτxz`, `ρτyz`) and the buoyancy flux (`buoy_flux`).
These are updated in place when the coupler computes turbulent fluxes,
rather than in `canopy_boundary_fluxes!`.

Note that we currently store these in the land model because the coupler
computes turbulent land/atmosphere fluxes using ClimaLand functions, and
the land model needs to be able to store the fluxes as an intermediary.
Once we compute fluxes entirely within the coupler, we can remove this.
"""
boundary_var_types(
    ::CanopyModel{FT},
    ::AtmosDrivenCanopyBC{<:CoupledAtmosphere, <:CoupledRadiativeFluxes},
    ::ClimaLand.TopBoundary,
) where {FT} = (
    NamedTuple{
        (
            :lhf,
            :shf,
            :transpiration,
            :r_ae,
            :∂LHF∂qc,
            :∂SHF∂Tc,
            :ρτxz,
            :ρτyz,
            :buoy_flux,
        ),
        Tuple{FT, FT, FT, FT, FT, FT, FT, FT, FT},
    },
)

"""
    SurfaceFluxes

## Interface
  - [`surface_conditions`](@ref) computes
    - Monin-Obukhov length
    - Potential temperature flux (if not given) using Monin-Obukhov theory
    - transport fluxes using Monin-Obukhov theory
    - friction velocity/temperature scale/tracer scales
    - exchange coefficients

## References
 - [Nishizawa2018](@cite)
 - [Byun1990](@cite)

"""
module SurfaceFluxes


import RootSolvers
const RS = RootSolvers

using DocStringExtensions
const DSE = DocStringExtensions

import Thermodynamics
const TD = Thermodynamics

include("UniversalFunctions.jl")
import .UniversalFunctions
const UF = UniversalFunctions

include("Parameters.jl")
import .Parameters

const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

abstract type SolverScheme end
struct LayerAverageScheme <: SolverScheme end
struct PointValueScheme <: SolverScheme end

"""
    SurfaceFluxConditions

Surface flux conditions, returned from `surface_conditions`.

# Fields

$(DSE.FIELDS)
"""
struct SurfaceFluxConditions{FT <: Real}
    L_MO::FT
    shf::FT
    lhf::FT
    buoy_flux::FT
    ρτxz::FT
    ρτyz::FT
    ustar::FT
    Cd::FT
    Ch::FT
    evaporation::FT
end

SurfaceFluxConditions(L_MO, shf, lhf, buoy_flux, ρτxz, ρτyz, ustar, Cd, Ch, E) =
    SurfaceFluxConditions(
        promote(L_MO, shf, lhf, buoy_flux, ρτxz, ρτyz, ustar, Cd, Ch, E)...,
    )

function Base.show(io::IO, sfc::SurfaceFluxConditions)
    println(io, "----------------------- SurfaceFluxConditions")
    println(io, "L_MO                   = ", sfc.L_MO)
    println(io, "Sensible Heat Flux     = ", sfc.shf)
    println(io, "Latent Heat Flux       = ", sfc.lhf)
    println(io, "Buoyancy Flux          = ", sfc.buoy_flux)
    println(io, "Friction velocity u⋆   = ", sfc.ustar)
    println(io, "C_drag                 = ", sfc.Cd)
    println(io, "C_heat                 = ", sfc.Ch)
    println(io, "evaporation            = ", sfc.evaporation)
    println(io, "-----------------------")
end

"""
   StateValues

Input container for state variables at either first / interior nodes.

# Fields

$(DSE.FIELDS)
"""
struct StateValues{FT <: Real, A, TS <: TD.ThermodynamicState}
    z::FT
    u::A
    ts::TS
end

abstract type AbstractSurfaceConditions{
    FT <: Real,
    SVA <: StateValues,
    SVB <: StateValues,
} end

"""
    Fluxes

Input container for state variables, latent and sensible heat fluxes roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct Fluxes{FT, SVA, SVB} <: AbstractSurfaceConditions{FT, SVA, SVB}
    state_in::SVA
    state_sfc::SVB
    shf::FT
    lhf::FT
    z0m::FT
    z0b::FT
    gustiness::FT
end

function Fluxes(
    state_in::SVA,
    state_sfc::SVB,
    shf::FT,
    lhf::FT,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1),
) where {SVA, SVB, FT}
    return Fluxes{FT, SVA, SVB}(
        state_in,
        state_sfc,
        shf,
        lhf,
        z0m,
        z0b,
        gustiness,
    )
end


"""
    FluxesAndFrictionVelocity

Input container, given surface state variables, latent and sensible heat fluxes,
and the friction velocity, roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct FluxesAndFrictionVelocity{FT, SVA, SVB} <:
       AbstractSurfaceConditions{FT, SVA, SVB}
    state_in::SVA
    state_sfc::SVB
    shf::FT
    lhf::FT
    ustar::FT
    z0m::FT
    z0b::FT
    gustiness::FT
end

function FluxesAndFrictionVelocity(
    state_in::SVA,
    state_sfc::SVB,
    shf::FT,
    lhf::FT,
    ustar::FT,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1),
) where {SVA, SVB, FT}
    return FluxesAndFrictionVelocity{FT, SVA, SVB}(
        state_in,
        state_sfc,
        shf,
        lhf,
        ustar,
        z0m,
        z0b,
        gustiness,
    )
end

"""
    Coefficients

Input container, given surface state variables, and exchange coefficients,roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct Coefficients{FT, SVA, SVB} <: AbstractSurfaceConditions{FT, SVA, SVB}
    state_in::SVA
    state_sfc::SVB
    Cd::FT
    Ch::FT
    gustiness::FT
    beta::FT
end

function Coefficients(
    state_in::SVA,
    state_sfc::SVB,
    Cd::FT,
    Ch::FT;
    gustiness::FT = FT(1),
    beta::FT = FT(1),
) where {SVA, SVB, FT}
    return Coefficients{FT, SVA, SVB}(
        state_in,
        state_sfc,
        Cd,
        Ch,
        gustiness,
        beta,
    )
end


"""
    ValuesOnly

Input container, given only surface state variables, roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct ValuesOnly{FT, SVA, SVB} <: AbstractSurfaceConditions{FT, SVA, SVB}
    state_in::SVA
    state_sfc::SVB
    z0m::FT
    z0b::FT
    gustiness::FT
    beta::FT
end

function ValuesOnly(
    state_in::SVA,
    state_sfc::SVB,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1),
    beta::FT = FT(1),
) where {SVA, SVB, FT}
    return ValuesOnly{FT, SVA, SVB}(
        state_in,
        state_sfc,
        z0m,
        z0b,
        gustiness,
        beta,
    )
end

ts_in(sc::AbstractSurfaceConditions) = sc.state_in.ts
ts_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.ts

z_in(sc::AbstractSurfaceConditions) = sc.state_in.z
z_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.z
Δz(sc::AbstractSurfaceConditions) = z_in(sc) - z_sfc(sc)

z0(sc::AbstractSurfaceConditions, ::UF.HeatTransport) = sc.z0b
z0(sc::AbstractSurfaceConditions, ::UF.MomentumTransport) = sc.z0m

Δu1(sc::AbstractSurfaceConditions) = sc.state_in.u[1] - sc.state_sfc.u[1]
Δu2(sc::AbstractSurfaceConditions) = sc.state_in.u[2] - sc.state_sfc.u[2]

qt_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_in(sc))
qt_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_sfc(sc))
Δqt(param_set::APS, sc::AbstractSurfaceConditions) =
    qt_in(param_set, sc) - qt_sfc(param_set, sc)

u_in(sc::AbstractSurfaceConditions) = sc.state_in.u
u_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.u

function windspeed(sc::AbstractSurfaceConditions)
    return max(hypot(Δu1(sc), Δu2(sc)), sc.gustiness)
end

"""
    surface_conditions(
        param_set::AbstractSurfaceFluxesParameters,
        sc::SurfaceFluxes.AbstractSurfaceConditions,
        scheme::SurfaceFluxes.SolverScheme = PointValueScheme();
        tol_neutral = SFP.cp_d(param_set) / 100,
        tol::RS.AbstractTolerance = RS.RelativeSolutionTolerance(FT(0.01)),
        maxiter::Int = 10,
        soltype::RS.SolutionType = RS.CompactSolution(),
        noniterative_stable_sol::Bool=true,
    )

The main user facing function of the module.
It computes the surface conditions
based on the Monin-Obukhov similarity functions. Requires
information about thermodynamic parameters (`param_set`)
the surface state `sc`, the universal function type and
the discretisation `scheme`. Default tolerance for
Monin-Obukhov length is absolute (i.e. has units [m]).
Returns the RootSolvers `CompactSolution` by default.

Result struct of type SurfaceFluxConditions contains:
  - L_MO:   Monin-Obukhov lengthscale
  - shf:    Sensible Heat Flux
  - lhf:    Latent Heat Flux
  - ρτxz:   Momentum Flux (Eastward component)
  - ρτyz:   Momentum Flux (Northward component)
  - ustar:  Friction velocity
  - Cd:     Momentum Exchange Coefficient
  - Ch:     Thermal Exchange Coefficient
"""
function surface_conditions(
    param_set::APS{FT},
    sc::AbstractSurfaceConditions,
    scheme::SolverScheme = PointValueScheme();
    tol_neutral = SFP.cp_d(param_set) / 100,
    tol::RS.AbstractTolerance = RS.RelativeSolutionTolerance(FT(0.01)),
    maxiter::Int = 10,
    soltype::RS.SolutionType = RS.CompactSolution(),
    noniterative_stable_sol::Bool = true,
) where {FT}
    uft = SFP.universal_func_type(param_set)
    L_MO = obukhov_length(
        param_set,
        sc,
        uft,
        scheme,
        tol,
        tol_neutral,
        maxiter,
        soltype,
        noniterative_stable_sol,
    )
    ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
    Cd = momentum_exchange_coefficient(
        param_set,
        L_MO,
        sc,
        uft,
        scheme,
        tol_neutral,
    )
    Ch =
        heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)
    shf = sensible_heat_flux(param_set, Ch, sc, scheme)
    lhf = latent_heat_flux(param_set, Ch, sc, scheme)
    buoy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    ρτxz, ρτyz = momentum_fluxes(param_set, Cd, sc, scheme)
    E = evaporation(param_set, sc, Ch)
    return SurfaceFluxConditions(
        L_MO,
        shf,
        lhf,
        buoy_flux,
        ρτxz,
        ρτyz,
        ustar,
        Cd,
        Ch,
        E,
    )
end

"""
    obukhov_length(sfc::SurfaceFluxConditions)

    obukhov_length( # internal method
        param_set::AbstractSurfaceFluxesParameters,
        sc::AbstractSurfaceConditions,
        uft,
        scheme,
        tol,
        tol_neutral,
        maxiter,
        soltype,
    )

Compute and return the Monin-Obukhov lengthscale (LMO).

The internal method for computing LMO depends on the
particular surface condition `sc`, of which there are
several options:

 - `FluxesAndFrictionVelocity`
 - `Coefficients`

## `AbstractSurfaceConditions` (fallback)

The Monin-Obukhov length is computed by solving a non-linear
equation given a tolerance `tol` and maximum iterations `maxiter`.

## `FluxesAndFrictionVelocity`

Surface fluxes and friction velocity are known.
Iterations are not needed to determine LMO.

## `Coefficients`

Exchange coefficients are known.
Iterations are not needed to determine LMO.
"""
function obukhov_length end

obukhov_length(sfc::SurfaceFluxConditions) = sfc.L_MO

function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end

function compute_richardson_number(
    sc::AbstractSurfaceConditions,
    DSEᵥ_in,
    DSEᵥ_sfc,
    grav,
)
    return (grav * Δz(sc) * (DSEᵥ_in - DSEᵥ_sfc)) /
           (DSEᵥ_in * (windspeed(sc))^2)
end

function compute_∂Ri∂ζ(param_set, sc::AbstractSurfaceConditions, uft, scheme, ζ)
    # In this design, this ∂Ri∂ζ function is intended to be an
    # internal function to support the Newton iteration scheme
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    dz = Δz(sc)
    z0m = z0(sc, UF.MomentumTransport())
    z0h = z0(sc, UF.HeatTransport())
    ufₛ = UF.universal_func(uft, Δz(sc) / ζ, SFP.uf_params(param_set))
    ψₘ = UF.psi(ufₛ, ζ, UF.MomentumTransport())
    ψₕ = UF.psi(ufₛ, ζ, UF.HeatTransport())
    ψₘ₀ = UF.psi(ufₛ, z0m * ζ / dz, UF.MomentumTransport())
    ψₕ₀ = UF.psi(ufₛ, z0h * ζ / dz, UF.HeatTransport())
    ϕₘ = UF.phi(ufₛ, ζ, UF.MomentumTransport())
    ϕₕ = UF.phi(ufₛ, ζ, UF.HeatTransport())
    ϕₘ₀ = UF.phi(ufₛ, z0m * ζ / dz, UF.MomentumTransport())
    ϕₕ₀ = UF.phi(ufₛ, z0h * ζ / dz, UF.HeatTransport())
    F_m = log(dz / z0m) - ψₘ + ψₘ₀
    F_h = log(dz / z0h) - ψₕ + ψₕ₀
    ∂Ri∂ζ =
        compute_Ri_b(param_set, sc, uft, scheme, ζ) / ζ *
        (1 + 1 / F_h * (ϕₕ - ϕₕ₀) - 2 * (ϕₘ - ϕₘ₀) * (1 / F_m))
end

function compute_Ri_b(param_set, sc::AbstractSurfaceConditions, uft, scheme, ζ)
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    ufₛ = UF.universal_func(uft, Δz(sc) / ζ, SFP.uf_params(param_set))
    F_m = compute_Fₘₕ(sc, ufₛ, ζ, UF.MomentumTransport())
    F_h = compute_Fₘₕ(sc, ufₛ, ζ, UF.HeatTransport())
    return ζ * F_h / F_m^2
end

function compute_Fₘₕ(sc, ufₛ, ζ, transport)
    ψ = UF.psi(ufₛ, ζ, transport)
    ψ₀ = UF.psi(ufₛ, z0(sc, transport) * ζ / Δz(sc), transport)
    return log(Δz(sc) / z0(sc, transport)) - ψ + ψ₀
end

function obukhov_length(
    param_set::APS{FT},
    sc::Union{Fluxes, ValuesOnly},
    uft::UF.AUFT,
    scheme,
    tol,
    tol_neutral,
    maxiter,
    soltype,
    noniterative_stable_sol,
) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    grav = SFP.grav(param_set)
    DSEᵥ_in =
        TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    DSEᵥ_sfc = TD.virtual_dry_static_energy(
        thermo_params,
        ts_sfc(sc),
        grav * z_sfc(sc),
    )
    ΔDSEᵥ = DSEᵥ_in - DSEᵥ_sfc
    if ΔDSEᵥ >= 0 && noniterative_stable_sol == true # Stable Layer
        ### Analytical Solution
        ### Gryanik et al. (2021)
        ### DOI: 10.1029/2021MS002590)
        Ri_b = compute_richardson_number(sc, DSEᵥ_in, DSEᵥ_sfc, grav)
        ϵₘ = Δz(sc) / z0(sc, UF.MomentumTransport())
        ϵₕ = Δz(sc) / z0(sc, UF.HeatTransport())
        C = (log(ϵₘ))^2 / (log(ϵₕ))
        ζₐ = ufparams.ζ_a
        ufₐ = UF.universal_func(uft, Δz(sc) / ζₐ, SFP.uf_params(param_set))
        γ = ufparams.γ
        Pr₀ = UF.Pr_0(ufₐ)
        ψ_ma = UF.psi(ufₐ, ζₐ, UF.MomentumTransport())
        ψ_ha = UF.psi(ufₐ, ζₐ, UF.HeatTransport())
        A =
            (log(ϵₘ) - ψ_ma)^(2 * (γ - 1)) /
            (ζₐ^(γ - 1) * (log(ϵₕ) - ψ_ha)^(γ - 1)) *
            ((log(ϵₘ) - ψ_ma)^2 / (log(ϵₕ) - ψ_ha) - C)
        ζₛ = C * Ri_b + A * Ri_b^γ
        ufₛ = UF.universal_func(uft, Δz(sc) / ζₛ, SFP.uf_params(param_set))
        ψₘ = UF.psi(ufₛ, ζₛ, UF.MomentumTransport())
        ψₕ = UF.psi(ufₛ, ζₛ, UF.HeatTransport())
        # Compute exchange coefficients
        κ = SFP.von_karman_const(param_set)
        Cd = κ^2 / (log(ϵₘ)^2) * (1 - ψₘ / log(ϵₘ))^(-2)
        Ch =
            κ^2 / (Pr₀ * log(ϵₘ) * log(ϵₕ)) *
            (1 - ψₘ / log(ϵₘ))^(-1) *
            (1 - ψₕ / Pr₀ / log(ϵₕ))^(-1)
        return non_zero(Δz(sc) / ζₛ)
    elseif tol_neutral >= abs(ΔDSEᵥ) # Neutral Layer
        # Large L_MO -> virtual dry static energy suggests neutral boundary layer
        # Return ζ->0 in the neutral boundary layer case, where ζ = z / L_MO
        return L_MO = FT(Inf) * sign(non_zero(ΔDSEᵥ))
    else
        function root_ζ(ζ)
            f =
                compute_richardson_number(sc, DSEᵥ_in, DSEᵥ_sfc, grav) -
                compute_Ri_b(param_set, sc, uft, scheme, ζ)
            return f
        end
        function root_and_deriv_ζ(ζ)
            f = root_ζ(ζ)
            f′ = -compute_∂Ri∂ζ(param_set, sc, uft, scheme, ζ)
            return (f, f′)
        end
        ζ₀ = sign(ΔDSEᵥ)
        sol = RS.find_zero(
            root_and_deriv_ζ,
            RS.NewtonsMethod(ζ₀),
            soltype,
            tol,
            maxiter,
        )
        L_MO = Δz(sc) / sol.root
        return non_zero(L_MO)
    end
end

function obukhov_length(
    param_set,
    sc::FluxesAndFrictionVelocity,
    uft::UF.AUFT,
    scheme,
    args...,
)
    return -sc.ustar^3 / SFP.von_karman_const(param_set) /
           non_zero(compute_buoyancy_flux(param_set, sc, scheme))
end

function obukhov_length(
    param_set,
    sc::Coefficients,
    uft::UF.AUFT,
    scheme,
    args...,
)
    lhf = latent_heat_flux(param_set, sc.Ch, sc, scheme)
    shf = sensible_heat_flux(param_set, sc.Ch, sc, scheme)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    buoyancy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    return -ustar^3 / SFP.von_karman_const(param_set) / non_zero(buoyancy_flux)
end

"""
    compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)

Returns the buoyancy flux when the surface fluxes are known.
"""
function compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    ε_vd = SFP.Rv_over_Rd(param_set)
    cp_m = TD.cp_m(thermo_params, ts_in)
    L_v = TD.latent_heat_vapor(thermo_params, ts_in)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc)
    T_in = TD.air_temperature(thermo_params, ts_in)
    return grav / ρ_sfc * (shf / cp_m / T_in + (ε_vd - 1) * lhf / L_v)
end

function compute_buoyancy_flux(
    param_set,
    sc::Union{FluxesAndFrictionVelocity, Fluxes},
    scheme,
)
    return compute_buoyancy_flux(
        param_set,
        sc.shf,
        sc.lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
end

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Returns buoyancy star based on known air densities.
"""
function compute_bstar(
    param_set,
    L_MO,
    sc::AbstractSurfaceConditions,
    uft::UF.AUFT,
    scheme,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    grav = SFP.grav(param_set)
    DSEᵥ_sfc = TD.virtual_dry_static_energy(
        thermo_params,
        ts_sfc(sc),
        grav * z_sfc(sc),
    )
    DSEᵥ_in =
        TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    DSEᵥ_star =
        compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            UF.HeatTransport(),
            uft,
            scheme,
        ) * (DSEᵥ_in - DSEᵥ_sfc)
    return grav * DSEᵥ_star / DSEᵥ_in
end

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Returns buoyancy star based on known friction velocity  and fluxes.
"""
compute_bstar(
    param_set,
    L_MO,
    sc::FluxesAndFrictionVelocity,
    uft::UF.AUFT,
    scheme,
) = -compute_buoyancy_flux(param_set, sc, scheme) / sc.ustar

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Return buoyancy star based on known fluxes.
"""
compute_bstar(param_set, L_MO, sc::Fluxes, uft::UF.AUFT, scheme) =
    -compute_buoyancy_flux(param_set, sc, scheme) /
    compute_ustar(param_set, L_MO, sc, uft, scheme)


"""
    compute_ustar(
        param_set::AbstractSurfaceFluxesParameters,
        L_MO,
        sc::AbstractSurfaceCondition,
        uft,
        scheme
    )

Return the friction velocity. This method is dispatched
by the surface condition:

## `sc::FluxesAndFrictionVelocity`

Friction velocity is known.

## `sc::Fluxes`

Compute given the Monin-Obukhov lengthscale.

## `sc::Coefficients`

Compute given the exchange coefficients.

## `sc::ValuesOnly`
Compute given the Monin-Obukhov lengthscale.
"""
function compute_ustar end

compute_ustar(param_set, L_MO, sc::FluxesAndFrictionVelocity, uft, scheme) =
    sc.ustar

compute_ustar(param_set, L_MO, sc::Fluxes, uft, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        UF.MomentumTransport(),
        uft,
        scheme,
    )

compute_ustar(param_set, L_MO, sc::Coefficients, uft, scheme) =
    sqrt(sc.Cd) * (windspeed(sc))

compute_ustar(param_set, L_MO, sc::ValuesOnly, uft, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        UF.MomentumTransport(),
        uft,
        scheme,
    )

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Compute and return Cd, the momentum exchange coefficient, given the
Monin-Obukhov lengthscale.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    uft::UF.AUFT,
    scheme,
    tol_neutral,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    κ = SFP.von_karman_const(param_set)
    grav = SFP.grav(param_set)
    transport = UF.MomentumTransport()
    DSEᵥ_sfc = TD.virtual_dry_static_energy(
        thermo_params,
        ts_sfc(sc),
        grav * z_sfc(sc),
    )
    DSEᵥ_in =
        TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    ΔDSEᵥ = DSEᵥ_in - DSEᵥ_sfc
    if abs(ΔDSEᵥ) <= tol_neutral
        Cd = (κ / log(Δz(sc) / z0(sc, transport)))^2
    else
        ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
        Cd = ustar^2 / windspeed(sc)^2
    end
    return Cd
end

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)

Return Cd, the momentum exchange coefficient.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    sc::Coefficients,
    uft,
    scheme,
    tol_neutral,
)
    return sc.Cd
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)

Compute and return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    uft,
    scheme,
    tol_neutral,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    transport = UF.HeatTransport()
    κ = SFP.von_karman_const(param_set)
    grav = SFP.grav(param_set)
    DSEᵥ_sfc = TD.virtual_dry_static_energy(
        thermo_params,
        ts_sfc(sc),
        grav * z_sfc(sc),
    )
    DSEᵥ_in =
        TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    ΔDSEᵥ = DSEᵥ_in - DSEᵥ_sfc
    z0_b = z0(sc, UF.HeatTransport())
    z0_m = z0(sc, UF.MomentumTransport())
    if abs(ΔDSEᵥ) <= tol_neutral
        Ch = κ^2 / (log(Δz(sc) / z0_b) * log(Δz(sc) / z0_m))
    else
        ϕ_heat = compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            transport,
            uft,
            scheme,
        )
        ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
        Ch = ustar * ϕ_heat / windspeed(sc)
    end
    return Ch
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    sc::Coefficients,
    uft,
    scheme,
    tol_neutral,
)
    return sc.Ch
end

"""
    momentum_fluxes(param_set, Cd, sc, scheme)

Compute and return the momentum fluxes
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Cd: Momentum exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function momentum_fluxes(param_set, Cd, sc::AbstractSurfaceConditions, scheme)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    ρτxz = -ρ_sfc * Cd * Δu1(sc) * windspeed(sc)
    ρτyz = -ρ_sfc * Cd * Δu2(sc) * windspeed(sc)
    return (ρτxz, ρτyz)
end

"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

Compute and return the sensible heat fluxes
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Ch: Thermal exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly, Coefficients},
    scheme,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    cp_d = SFP.cp_d(param_set)
    R_d = SFP.R_d(param_set)
    T_0 = SFP.T_0(param_set)
    cp_m_in = TD.cp_m(thermo_params, ts_in(sc))
    cp_m_sfc = TD.cp_m(thermo_params, ts_sfc(sc))
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    T_in = TD.air_temperature(thermo_params, ts_in(sc))
    T_sfc = TD.air_temperature(thermo_params, ts_sfc(sc))
    ΔΦ = grav * Δz(sc)
    ΔDSE = cp_m_in * (T_in - T_0) - cp_m_sfc * (T_sfc - T_0) + ΔΦ
    return -ρ_sfc * Ch * windspeed(sc) * ΔDSE
end

"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known sensible heat flux.
"""
function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    scheme,
)
    return sc.shf
end

"""
    latent_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known latent heat flux.
"""
function latent_heat_flux(
    param_set,
    L_MO,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    scheme,
)
    return sc.lhf
end

"""
    latent_heat_flux(param_set, Ch, sc, scheme)

Compute and return the latent heat flux
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Ch: Thermal exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function latent_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly, Coefficients},
    scheme,
)
    Lv_0 = SFP.LH_v0(param_set)
    E = evaporation(param_set, sc, Ch)
    lhf = Lv_0 * E
    return lhf
end

"""
    evaporation(param_set, sc, Ch)

Compute and return the evaporation. When `sc` is `Fluxes` or `FluxesAndFrictionVelocity`,
evaporation is directly calculated from the latent heat flux.
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - Ch: Thermal exchange coefficient
"""
function evaporation(
    param_set,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    Ch,
)
    Lv_0 = SFP.LH_v0(param_set)
    return sc.lhf / Lv_0
end

"""
    evaporation(param_set, sc, Ch)

Compute and return the evaporation. When `sc` is `ValuesOnly` or `Coefficients`, a `beta` factor
is used to represent the resistance of the surface.
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - Ch: Thermal exchange coefficient
"""
function evaporation(param_set, sc::Union{ValuesOnly, Coefficients}, Ch)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    return -ρ_sfc * Ch * windspeed(sc) * Δqt(param_set, sc) * sc.beta
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, ::LayerAverageScheme)

Computes the coefficient for the physical scale of a variable based on Nishizawa(2018)
for the FV scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g. Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set::APS,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    transport,
    uft,
    ::LayerAverageScheme,
)
    von_karman_const = SFP.von_karman_const(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    π_group = UF.π_group(uf, transport)
    R_z0 = 1 - z0(sc, transport) / Δz(sc)
    denom1 = log(Δz(sc) / z0(sc, transport))
    denom2 = -UF.Psi(uf, Δz(sc) / uf.L, transport)
    denom3 =
        z0(sc, transport) / Δz(sc) *
        UF.Psi(uf, z0(sc, transport) / uf.L, transport)
    denom4 = R_z0 * (UF.psi(uf, z0(sc, transport) / uf.L, transport) - 1)
    Σterms = denom1 + denom2 + denom3 + denom4
    return von_karman_const / (π_group * Σterms)
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, ::PointValueScheme)

Computes the coefficient for the physical scale of a variable based on Byun (1990)
for the Finite Differences scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g. Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    transport,
    uft,
    ::PointValueScheme,
)
    von_karman_const = SFP.von_karman_const(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    π_group = UF.π_group(uf, transport)
    denom1 = log(Δz(sc) / z0(sc, transport))
    denom2 = -UF.psi(uf, Δz(sc) / uf.L, transport)
    denom3 = UF.psi(uf, z0(sc, transport) / uf.L, transport)
    Σterms = denom1 + denom2 + denom3
    return von_karman_const / (π_group * Σterms)
end

"""
    recover_profile(param_set, sc, L_MO, Z, X_in, X_sfc, transport, scheme)

Recover profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_star: Scale parameter for variable X
  - X_sfc: For variable X, values at surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretization scheme (currently supports FD and FV)

# TODO: add tests
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO,
    Z,
    X_star,
    X_sfc,
    transport,
    scheme::Union{LayerAverageScheme, PointValueScheme},
)
    ufp = SFP.uf_params(param_set)
    uft = UF.universal_func_type(typeof(ufp))
    uf = UF.universal_func(uft, L_MO, ufp)
    von_karman_const = SFP.von_karman_const(param_set)
    num1 = log(Z / z0(sc, transport))
    num2 = -UF.psi(uf, Z / L_MO, transport)
    num3 = UF.psi(uf, z0(sc, transport) / L_MO, transport)
    Σnum = num1 + num2 + num3
    return Σnum * X_star / von_karman_const + X_sfc
end

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # SurfaceFluxes module

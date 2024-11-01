module Snow

using DocStringExtensions
import ...Parameters as LP
using ClimaCore
using Thermodynamics
using ClimaLand
using ClimaLand:
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    turbulent_fluxes,
    net_radiation,
    AbstractModel,
    heaviside

import ClimaLand:
    AbstractBC,
    make_update_aux,
    make_compute_exp_tendency,
    make_update_boundary_fluxes,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    prognostic_domain_names,
    auxiliary_types,
    auxiliary_domain_names,
    surface_temperature,
    surface_specific_humidity,
    surface_height,
    surface_albedo,
    surface_emissivity,
    get_drivers
export SnowParameters, SnowModel, AtmosDrivenSnowBC, snow_boundary_fluxes!

"""
    AbstractSnowModel{FT} <: ClimaLand.AbstractExpModel{FT}

Defines a new type of abstract explicit model for snow modeling.
Currently, the only supported concrete example is called `SnowModel`
and is used as a bulk snow model.
"""
abstract type AbstractSnowModel{FT} <: ClimaLand.AbstractExpModel{FT} end

"""
    AbstractDensityModel{FT}
Defines the model type for density and depth parameterizations
for use within an `AbstractSnowModel` type. Current examples include the
`ConstantDensityModel` and the `Anderson1976` models.
"""
abstract type AbstractDensityModel{FT <: AbstractFloat} end

"""
    ConstantDensityModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density
is always treated as a constant (type FT).
"""
struct ConstantDensityModel{FT} <: AbstractDensityModel{FT}
    ρ_snow::FT
end

"""
    ConstantDensityModel{FT}(ρ::FT)
An outer constructor for the `ConstantDensityModel` density parameterization for usage in a snow model.
"""
function ConstantDensityModel{FT}(ρ::FT) where {FT}
    return ConstantDensityModel{FT}(ρ)
end

"""
    ConstantDensityModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density
compacts according to the seminal works of Anderson in the 1970s (and specifically the
numerical implementation documented in the Snow17 model).
"""
struct Anderson1976{FT} <: AbstractDensityModel{FT}
    c1::FT
    c2::FT
    c3::FT
    c4::FT
    c5::FT
    cx::FT
    ρ_d::FT
end

"""
    Anderson1976(FT; c1, c2, c3, c4, c5, ρ_d, cx)
An outer constructor for the `Anderson1976` density parameterization for usage in a snow model.
Uses the default constants defined by Anderson but allows for redefinition of the model constants.
These include:
- c1: fractional increase in density (0.026 cm/hr)
- c2: compression constant estimated by Kojima in 1967 (21 cm^3/g)
- c3: fractional settling rate at 0 degrees C for densities less than the critical density ρ_d (0.005 1/hr)
- c4: constant (0.1 1/degrees C)
- c5: scaling of the settling rate when no water is present in the snowpack (0, it is 2 when water is present)
- cx: destructive metamorphism decay factor for densities greater than the critical density ρ_d (23)
- ρ_d: the critical density (unitless as a fraction of the density of liquid water, 0.15)
"""
function Anderson1976(
    FT::DataType;
    c1 = 0.026,
    c2 = 21.0,
    c3 = 0.005,
    c4 = 0.10,
    c5 = 0.0, #Snow17 code says 1.0, paper says 0.0?
    ρ_d = 0.15,
    cx = 23.0,
)
    return Anderson1976{FT}(
        FT(c1),
        FT(c2),
        FT(c3),
        FT(c4),
        FT(c5),
        FT(ρ_d),
        FT(cx),
    )
end


"""
    SnowParameters{FT <: AbstractFloat, PSE}

A struct for storing parameters of the `SnowModel`.

Note that in our current implementation of runoff, a physical
timescale is required and computed using Ksat and the depth
of the snow. For shallow snowpacks, this will fall below the timestep
of the model. For that reason, we pass the timestep of the model as 
a parameter, and take the larger of the timestep and the physical timescale
as the value used in the model. Future implementations will revisit this.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SnowParameters{
    FT <: AbstractFloat,
    DM <: AbstractDensityModel,
    PSE,
}
    "Choice of parameterization for snow density"
    density::DM
    "Roughness length over snow for momentum (m)"
    z_0m::FT
    "Roughness length over snow for scalars (m)"
    z_0b::FT
    "Albedo of snow (unitless)"
    α_snow::FT
    "Emissivity of snow (unitless)"
    ϵ_snow::FT
    "Volumetric holding capacity of water in snow (unitless)"
    θ_r::FT
    "Hydraulic conductivity of wet snow (m/s)"
    Ksat::FT
    "Thermal conductivity of ice (W/m/K)"
    κ_ice::FT
    "Timestep of the model (s)"
    Δt::FT
    "Areal specific heat of ground interacting with snow (J/m^2/K)"
    ρcD_g::FT
    "Clima-wide parameters"
    earth_param_set::PSE
end

"""
   SnowParameters{FT}(Δt;
                      density = ConstantDensityModel(200),
                      z_0m = FT(0.0024),
                      z_0b = FT(0.00024),
                      α_snow = FT(0.8),
                      ϵ_snow = FT(0.99),
                      θ_r = FT(0.08),
                      Ksat = FT(1e-3),
                      κ_ice = FT(2.21),
                      ρcD_g = FT(3.553e5),
                      earth_param_set::PSE) where {FT, PSE}

An outer constructor for `SnowParameters` which supplies defaults for
all arguments but `earth_param_set`.
"""
function SnowParameters{FT}(
    Δt;
    density::DM = ConstantDensityModel(FT(200)),
    z_0m = FT(0.0024),
    z_0b = FT(0.00024),
    α_snow = FT(0.8),
    ϵ_snow = FT(0.99),
    θ_r = FT(0.08),
    Ksat = FT(1e-3),
    κ_ice = FT(2.21),
    ρcD_g = FT(3.553e5),
    earth_param_set::PSE,
) where {FT <: AbstractFloat, DM <: AbstractDensityModel, PSE}
    return SnowParameters{FT, DM, PSE}(
        density,
        z_0m,
        z_0b,
        α_snow,
        ϵ_snow,
        θ_r,
        Ksat,
        κ_ice,
        Δt,
        ρcD_g,
        earth_param_set,
    )
end

Base.broadcastable(ps::SnowParameters) = tuple(ps)

"""
    struct SnowModel{
        FT,
        PS <: SnowParameters{FT},
        BC,
        D,
    } <: AbstractSnowModel{FT}

A container/type for the bulk snow model, based on the UEB snow model
of Tarboton et al. (1995) and Tarboton and Luce (1996).
"""
struct SnowModel{FT, PS <: SnowParameters{FT}, BC, D} <: AbstractSnowModel{FT}
    "Parameters required by the snow model"
    parameters::PS
    "Boundary conditions"
    boundary_conditions::BC
    "The domain of the model"
    domain::D
end

function SnowModel(;
    parameters::SnowParameters{FT, DM, PSE},
    domain::ClimaLand.Domains.AbstractDomain,
    boundary_conditions::BC,
) where {FT, DM, PSE, BC}
    args = (parameters, boundary_conditions, domain)
    SnowModel{FT, typeof.(args)...}(args...)
end

"""
    prognostic_vars(::SnowModel)

Returns the prognostic variable names of the snow model.

For this model, we track the snow water equivalent S in meters (liquid
water volume per ground area) and
the energy per unit ground area U [J/m^2] prognostically,
as well as the snow depth Z [m] for some types of density models.
"""
prognostic_vars(m::SnowModel) =
    (:S, :U, density_prog_vars(m.parameters.density)...)
density_prog_vars(m::ConstantDensityModel) = ()
density_prog_vars(m::Anderson1976) = (:Z)

"""
    prognostic_types(::SnowModel{FT})

Returns the prognostic variable types of the snow model;
both snow water equivalent and energy per unit area
are scalars.
"""
prognostic_types(m::SnowModel{FT}) where {FT} =
    (FT, FT, density_prog_types(m.parameters.density)...)
density_prog_types(m::ConstantDensityModel{FT}) where {FT} = ()
density_prog_types(m::Anderson1976{FT}) where {FT} = (FT)

"""
    prognostic_domain_names(::SnowModel)

Returns the prognostic variable domain names of the snow model;
both snow water equivalent and energy per unit area
are modeling only as a function of (x,y), and not as a function
of depth. Therefore their domain name is ":surface".
"""
prognostic_domain_names(m::SnowModel) =
    (:surface, :surface, density_prog_names(m.parameters.density)...)
density_prog_names(m::ConstantDensityModel) = ()
density_prog_names(m::Anderson1976) = (:surface)

"""
    auxiliary_vars(::SnowModel)

Returns the auxiliary variable names for the snow model. These
include the mass fraction in liquid water (`q_l`, unitless),
the thermal conductivity (`κ`, W/m/K),
the bulk temperature (`T`, K), the surface temperature (`T_sfc`, K), the bulk snow density (`ρ_snow`, kg/m^3)
the SHF, LHF, and vapor flux (`turbulent_fluxes.shf`, etc),
the net radiation (`R_n, J/m^2/s)`, the energy flux in liquid water runoff
(`energy_runoff`, J/m^2/s), the water volume in runoff (`water_runoff`, m/s), and the total energy and water fluxes applied to the snowpack.

Since the snow can melt completely in one timestep, we clip the water and energy fluxes
such that SWE cannot become negative and U cannot become unphysical. The
clipped values are what are actually applied as boundary fluxes, and are stored in
`applied_` fluxes.
"""
auxiliary_vars(m::SnowModel) = (
    :q_l,
    :κ,
    :T,
    :T_sfc,
    :ρ_snow,
    :turbulent_fluxes,
    :R_n,
    :energy_runoff,
    :water_runoff,
    :total_energy_flux,
    :total_water_flux,
    :applied_energy_flux,
    :applied_water_flux,
    :snow_cover_fraction,
    density_aux_vars(m.parameters.density)...,
)
density_aux_vars(m::Union{Anderson1976, ConstantDensityModel}) = () #do we want to include Z whenever it isn't a prognostic variable, if we already have SWE and ρ_snow?

auxiliary_types(m::SnowModel{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    NamedTuple{(:lhf, :shf, :vapor_flux, :r_ae), Tuple{FT, FT, FT, FT}},
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    density_aux_types(m.parameters.density)...,
)
density_aux_types(
    m::Union{ConstantDensityModel{FT}, Anderson1976{FT}},
) where {FT} = ()

auxiliary_domain_names(m::SnowModel) = (
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    density_aux_names(m.parameters.density)...,
)
density_aux_names(m::Union{ConstantDensityModel, Anderson1976}) = ()

ClimaLand.name(::SnowModel) = :snow

function ClimaLand.make_update_aux(model::SnowModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        parameters = model.parameters

        update_density!(parameters.density, parameters, Y, p)

        @. p.snow.κ = snow_thermal_conductivity(p.snow.ρ_snow, parameters)

        @. p.snow.q_l =
            snow_liquid_mass_fraction(Y.snow.U, Y.snow.S, parameters)

        @. p.snow.T =
            snow_bulk_temperature(Y.snow.U, Y.snow.S, p.snow.q_l, parameters)

        @. p.snow.T_sfc = snow_surface_temperature(p.snow.T)

        @. p.snow.water_runoff = compute_water_runoff(
            Y.snow.S,
            p.snow.q_l,
            p.snow.T,
            p.snow.ρ_snow,
            parameters,
        ) #do we want to dispatch this off the density model type so that it can take Y.snow.Z insetad of p.snow.ρ_snow when it exists to avoid a redundant computation?

        @. p.snow.energy_runoff =
            p.snow.water_runoff * volumetric_internal_energy_liq(FT, parameters)
        @. p.snow.snow_cover_fraction = snow_cover_fraction(Y.snow.S)
    end
end

function ClimaLand.make_update_boundary_fluxes(model::SnowModel{FT}) where {FT}
    function update_boundary_fluxes!(p, Y, t)
        # First compute the boundary fluxes
        snow_boundary_fluxes!(model.boundary_conditions, model, Y, p, t)
        # Next, clip them in case the snow will melt in this timestep
        @. p.snow.applied_water_flux = clip_total_snow_water_flux(
            Y.snow.S,
            p.snow.total_water_flux,
            model.parameters.Δt,
        )
        @. p.snow.applied_energy_flux = clip_total_snow_energy_flux(
            Y.snow.U,
            Y.snow.S,
            p.snow.total_energy_flux,
            p.snow.total_water_flux,
            model.parameters.Δt,
        )
    end
end

function ClimaLand.make_compute_exp_tendency(model::SnowModel{FT}) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        # positive fluxes are TOWARDS atmos; negative fluxes increase quantity in snow
        @. dY.snow.S = -p.snow.applied_water_flux
        @. dY.snow.U = -p.snow.applied_energy_flux
        update_density_prog!(model.parameters.density, model, dY, Y, p)
    end
    return compute_exp_tendency!
end

"""
    clip_total_snow_water_flux(S, total_water_flux, Δt)

A helper function which clips the total water flux so that
snow water equivalent S will not become negative in a timestep Δt.
"""
function clip_total_snow_water_flux(S, total_water_flux, Δt)
    if S - total_water_flux * Δt < 0
        return S / Δt
    else
        return total_water_flux
    end
end

"""
     clip_total_snow_energy_flux(U, S, total_energy_flux, total_water_flux, Δt)

A helper function which clips the total energy flux such that 
snow energy per unit ground area U will not become positive, and 
which ensures that if the snow water equivalent S goes to zero in a step, 
U will too.
"""
function clip_total_snow_energy_flux(
    U,
    S,
    total_energy_flux,
    total_water_flux,
    Δt,
)
    if (U - total_energy_flux * Δt) > 0
        return U / Δt
    elseif S - total_water_flux * Δt < 0
        return U / Δt
    else
        return total_energy_flux
    end
end

"""
    ClimaLand.get_drivers(model::SnowModel)

Returns the driver variable symbols for the SnowModel.
"""
function ClimaLand.get_drivers(model::SnowModel)
    return (
        model.boundary_conditions.atmos,
        model.boundary_conditions.radiation,
    )
end

include("./snow_parameterizations.jl")
include("./boundary_fluxes.jl")
end

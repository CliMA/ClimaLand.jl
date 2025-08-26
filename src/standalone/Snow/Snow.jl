module Snow
import ClimaParams as CP
using DocStringExtensions
import ...Parameters as LP
using ClimaCore
using LazyBroadcast: lazy
using Thermodynamics
using ClimaLand
using ClimaLand:
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    turbulent_fluxes!,
    net_radiation!,
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
    get_drivers,
    total_energy_per_area!,
    total_liq_water_vol_per_area!
export SnowParameters,
    SnowModel,
    AtmosDrivenSnowBC,
    snow_boundary_fluxes!,
    ConstantAlbedoModel,
    ZenithAngleAlbedoModel,
    WuWuSnowCoverFractionModel

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
for use within an `AbstractSnowModel` type. Currently we support a
`MinimumDensityModel` model and the `NeuralDepthModel`.

Since depth and bulk density are related via SWE, and SWE is a prognostic variable
of the snow model, only depth or bulk density can be independently modeled at one time.
This is why they are treated as a single "density model" even if the parameterization is actually
a model for snow depth.
The snow depth/density model can be diagnostic or introduce additional prognostic variables.

"""
abstract type AbstractDensityModel{FT <: AbstractFloat} end

"""
    MinimumDensityModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}

Establishes the density parameterization where snow density is estimated
using a linear interpolation between the minimum density and the density
of liquid water, based on the mass fraction of liquid water.
"""
struct MinimumDensityModel{FT} <: AbstractDensityModel{FT}
    "Minimum density of snow (kg/m^3)"
    ρ_min::FT
end


"""
    AbstractAlbedoModel{FT}

Defines the model type for albedo parameterization
for use within an `AbstractSnowModel` type. 

These parameterizations are stored in parameters.α_snow, and 
are used to update the value of p.snow.α_snow (the broadband
albedo of the snow at a point).
stored 
"""
abstract type AbstractAlbedoModel{FT <: AbstractFloat} end

"""
    ConstantAlbedoModel{FT <: AbstractFloat} <: AbstractAlbedoModel{FT}

Establishes the albedo parameterization where albedo is treated as a
constant spatially and temporally.
"""
struct ConstantAlbedoModel{FT} <: AbstractAlbedoModel{FT}
    "Albedo of snow (unitless)"
    α::FT
end

"""
    ZenithAngleAlbedoModel{FT <: AbstractFloat} <: AbstractAlbedoModel{FT}

Establishes the albedo parameterization where albedo
depends on the cosine of the zenith angle of the sun, as
    α = f(x) * [α_0 + Δα*exp(-k*cos(θs))],

where cos θs is the cosine of the zenith angle, α_0, Δα, and k 
are free parameters. The factor out front is a function of 
x = ρ_snow/ρ_liq, of the form f(x) = min(1 - β(x-x0), 1). The parameters
x0 ∈ [0,1] and β ∈ [0,1] are free. Choose β = 0 to remove this dependence on snow density.


Note: If this choice is used, the field cosθs must appear in the cache
p.drivers. This is available through the PrescribedRadiativeFluxes object.
"""
struct ZenithAngleAlbedoModel{FT} <: AbstractAlbedoModel{FT}
    "Free parameter controlling the minimum snow albedo"
    α_0::FT
    "Free parameter controlling the snow albedo when θs = 90∘"
    Δα::FT
    "Rate at which albedo drops to its minimum value with zenith angle"
    k::FT
    "Rate governing how snow albedo changes with snow density, a proxy for grain size and liquid water content, ∈ [0,1]"
    β::FT
    "Value of relative snow density ρ_snow/ρ_liq at which snow density begins to decrease albedo, ∈ [0,1]"
    x0::FT
end

function ZenithAngleAlbedoModel(
    α_0::FT,
    Δα::FT,
    k::FT;
    β = FT(0),
    x0 = FT(0.2),
) where {FT}
    @assert 0 ≤ x0 ≤ 1
    @assert 0 ≤ β ≤ 1
    ZenithAngleAlbedoModel(α_0, Δα, k, β, x0)
end

function ZenithAngleAlbedoModel(
    toml_dict::CP.AbstractTOMLDict;
    α_0 = toml_dict["alpha_0"],
    Δα = toml_dict["delta_alpha"],
    k = toml_dict["k"],
    β = toml_dict["beta"],
    x0 = toml_dict["x0"],
)
    return ZenithAngleAlbedoModel(α_0, Δα, k, β = β, x0 = x0)
end

"""
    AbstractSnowCoverFractionModel{FT}

Defines the model type for snow cover parameterization
for use within an `AbstractSnowModel` type. 

These parameterizations are stored in parameters.scf, and 
are used to update the value of p.snow.snow_cover_fraction.
"""
abstract type AbstractSnowCoverFractionModel{FT <: AbstractFloat} end

"""
    WuWuSnowCoverFractionModel{FT <: AbstractFloat} <: AbstractSnowCoverFractionModel{FT}

Establishes the snow cover parameterization of Wu, Tongwen, and 
Guoxiong Wu. "An empirical formula to compute
snow cover fraction in GCMs." Advances in Atmospheric Sciences
21 (2004): 529-535,
    scf = min(β_scf * z̃ / (z̃ + 1), 1),

where z̃ = snow depth per ground area / 0.106 m, and β_scf
is computed using a resolution dependent formula:
β_scf = max(β0 - γ(horz_degree_res - 1.5), β_min), where horz_degree_res is the
horizontal resolution of the simulation, in degrees, and β0, β_min and γ
are unitless. It is correct to think of β0, β_min, γ, and z0 as the free
parameters, while horz_degree_res is provided and β_scf is determined.

β0, β_min, γ, and β_scf must be > 0. 

From Wu and Wu et al, β0 ∼ 1.77 and γ ∼ 0.08, over a range of 1.5-4.5∘
"""
struct WuWuSnowCoverFractionModel{FT} <: AbstractSnowCoverFractionModel{FT}
    "The value used to normalize snow depth when computing snow cover fraction (m)"
    z0::FT
    "The value of β_scf (unitless; computed from other parameters)"
    β_scf::FT
    "Free parameter controlling the snow cover scaling change with resolution (1/degrees)"
    γ::FT
    "The value of β_scf at 1.5∘ horizontal resolution (unitless)"
    β0::FT
    "The minimum of β_scf as horizontal resolution gets coarser (unitless)"
    β_min::FT
    "The horizontal resolution of the model, in degrees"
    horz_degree_res::FT
    function WuWuSnowCoverFractionModel(
        z0::FT,
        β_scf::FT,
        γ::FT,
        β0::FT,
        β_min::FT,
        horz_degree_res::FT,
    ) where {FT}
        expected_β_scf = max(β0 - γ * (horz_degree_res - FT(1.5)), β_min)
        @assert z0 > eps(FT)
        @assert expected_β_scf ≈ β_scf
        @assert β0 > eps(FT)
        @assert β_min > eps(FT)
        @assert γ > eps(FT)
        @assert horz_degree_res > eps(FT)
        new{FT}(z0, β_scf, γ, β0, β_min, horz_degree_res)
    end
end

function WuWuSnowCoverFractionModel(
    γ::FT,
    β0::FT,
    β_min::FT,
    horz_degree_res::FT;
    z0 = FT(0.106),
) where {FT}
    @assert β_min > eps(FT)
    @assert β0 > eps(FT)
    @assert γ > eps(FT)
    @assert z0 > eps(FT)
    @assert horz_degree_res > eps(FT)
    β_scf = max(β0 - γ * (horz_degree_res - FT(1.5)), β_min)
    return WuWuSnowCoverFractionModel(z0, β_scf, γ, β0, β_min, horz_degree_res)
end

function WuWuSnowCoverFractionModel(
    toml_dict::CP.AbstractTOMLDict,
    horz_degree_res;
    γ = toml_dict["gamma"],
    β0 = toml_dict["beta_0"],
    β_min = toml_dict["beta_min"],
    z0 = toml_dict["z0"],
)
    return WuWuSnowCoverFractionModel(γ, β0, β_min, horz_degree_res; z0 = z0)
end

"""
    SnowParameters{FT <: AbstractFloat, PSE}

A struct for storing parameters of the `SnowModel`.

Note that in our current implementation of runoff, a physical timescale is
required and computed using Ksat and the depth of the snow. For shallow
snowpacks, this will fall below the timestep of the model. For that reason, we
pass the timestep of the model as a parameter, and take the larger of the
timestep and the physical timescale as the value used in the model. Future
implementations will revisit this.

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SnowParameters{
    FT <: AbstractFloat,
    DM <: AbstractDensityModel,
    AM <: AbstractAlbedoModel,
    SCFM <: AbstractSnowCoverFractionModel,
    PSE,
}
    "Choice of parameterization for snow density"
    density::DM
    "Roughness length over snow for momentum (m)"
    z_0m::FT
    "Roughness length over snow for scalars (m)"
    z_0b::FT
    "Albedo parameterization for snow"
    α_snow::AM
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
    "Parameter to prevent dividing by zero when computing snow temperature (m)"
    ΔS::FT
    "Snow cover fraction parameterization"
    scf::SCFM
    "Clima-wide parameters"
    earth_param_set::PSE
end

"""
   SnowParameters{FT}(Δt;
                      density = MinimumDensityModel(200),
                      z_0m = FT(0.0024),
                      z_0b = FT(0.00024),
                      α_snow = ConstantAlbedoModel(FT(0.8)),
                      ϵ_snow = FT(0.99),
                      θ_r = FT(0.08),
                      Ksat = FT(1e-3),
                      κ_ice = FT(2.21),
                      ΔS = FT(0.1),
                      scf = WuWuSnowCoverFractionModel(FT(0.106),FT(1.81), FT(0.08), FT(1.77), FT(1.0), FT(1),),
                      earth_param_set::PSE) where {FT, PSE}

An outer constructor for `SnowParameters` which supplies defaults for
all arguments but `earth_param_set`.
"""
function SnowParameters{FT}(
    Δt;
    density::DM = MinimumDensityModel(FT(200)),
    z_0m = FT(0.0024),
    z_0b = FT(0.00024),
    α_snow::AM = ConstantAlbedoModel(FT(0.8)),
    ϵ_snow = FT(0.99),
    θ_r = FT(0.08),
    Ksat = FT(1e-3),
    κ_ice = FT(2.21),
    ΔS = FT(0.1),
    scf::SCFM = WuWuSnowCoverFractionModel(
        FT(0.106),
        FT(1.81),
        FT(0.08),
        FT(1.77),
        FT(1),
        FT(1),
    ),
    earth_param_set::PSE,
) where {
    FT <: AbstractFloat,
    DM <: AbstractDensityModel,
    AM <: AbstractAlbedoModel,
    SCFM <: AbstractSnowCoverFractionModel,
    PSE,
}
    return SnowParameters{FT, DM, AM, SCFM, PSE}(
        density,
        z_0m,
        z_0b,
        α_snow,
        ϵ_snow,
        θ_r,
        Ksat,
        κ_ice,
        float(Δt),
        ΔS,
        scf,
        earth_param_set,
    )
end

## For interfacing with ClimaParams
"""
    function SnowParameters(
        FT,
        Δt;
        kwargs...  # For individual parameter overrides
    )

    function SnowParameters(
        toml_dict::CP.AbstractTOMLDict,
        Δt;
        kwargs...  # For individual parameter overrides
    )

Constructors for the SnowParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values. Possible calls:
```julia
Δt = 450.0
ClimaLand.Canopy.SnowParameters(Float64, Δt) # use the default values for all parameters
# Kwarg overrides
ClimaLand.Canopy.SnowParameters(Float64, Δt; ϵ_snow = 0.99)
# TOML Dictionary:
import ClimaParams as CP
toml_dict = CP.create_toml_dict(Float32);
ClimaLand.Canopy.SnowParameters(toml_dict, Δt; ϵ_snow = Float32(0.99), Ksat = Float32(1e-4))
```
"""
SnowParameters(::Type{FT}, Δt; kwargs...) where {FT <: AbstractFloat} =
    SnowParameters(CP.create_toml_dict(FT), Δt; kwargs...)

function SnowParameters(toml_dict::CP.AbstractTOMLDict, Δt; kwargs...)
    name_map = (;
        :snow_momentum_roughness_length => :z_0m,
        :snow_scalar_roughness_length => :z_0b,
        :thermal_conductivity_of_water_ice => :κ_ice,
        :snow_emissivity => :ϵ_snow,
        :holding_capacity_of_water_in_snow => :θ_r,
        :wet_snow_hydraulic_conductivity => :Ksat,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    name_map2 = (; :snow_density => :ρ_snow, :snow_albedo => :α_snow)
    parameters2 = CP.get_parameter_values(toml_dict, name_map2, "Land")

    density = MinimumDensityModel(parameters2.ρ_snow)
    α_snow = ConstantAlbedoModel(parameters2.α_snow)
    FT = CP.float_type(toml_dict)
    earth_param_set = LP.LandParameters(toml_dict)
    PSE = typeof(earth_param_set)
    return SnowParameters{FT}(
        Δt;
        earth_param_set,
        parameters...,
        density,
        α_snow,
        kwargs...,
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
    SnowModel(
        FT,
        domain,
        forcing,
        toml_dict::CP.AbstractTOMLDict,
        Δt;
        prognostic_land_components = (:snow,),
        z_0m = toml_dict["snow_momentum_roughness_length"],
        z_0b = toml_dict["snow_scalar_roughness_length"],
        ϵ_snow = toml_dict["snow_emissivity"],
        α_snow = ConstantAlbedoModel(toml_dict["snow_albedo"]),
        density = MinimumDensityModel(toml_dict["snow_density"]),
        scf = WuWuSnowCoverFractionModel(toml_dict, FT(1)),
        θ_r = toml_dict["holding_capacity_of_water_in_snow"],
        Ksat = toml_dict["wet_snow_hydraulic_conductivity"],
        ΔS = toml_dict["delta_S"],
    )

Creates a SnowModel model with the given float type FT, domain, toml_dict, forcing, and prognostic land components.

When running the snow model in standalone mode, provide `prognostic_land_components = (:snow,)`, while for running integrated 
land models, this should be a list of the component models. This value of this argument must be the same across all 
components in the integrated land model.

Default parameterizations and parameters can be overwritten using keyword arguments.
"""
function SnowModel(
    FT,
    domain,
    forcing,
    toml_dict::CP.AbstractTOMLDict,
    Δt;
    prognostic_land_components = (:snow,),
    z_0m = toml_dict["snow_momentum_roughness_length"],
    z_0b = toml_dict["snow_scalar_roughness_length"],
    ϵ_snow = toml_dict["snow_emissivity"],
    α_snow = ConstantAlbedoModel(toml_dict["snow_albedo"]),
    density = MinimumDensityModel(toml_dict["snow_density"]),
    scf = WuWuSnowCoverFractionModel(toml_dict, FT(1)),
    θ_r = toml_dict["holding_capacity_of_water_in_snow"],
    Ksat = toml_dict["wet_snow_hydraulic_conductivity"],
    ΔS = toml_dict["delta_S"],
)
    earth_param_set = LP.LandParameters(toml_dict)
    parameters = SnowParameters{FT}(
        Δt;
        earth_param_set,
        scf,
        α_snow,
        ϵ_snow,
        density,
        z_0m,
        z_0b,
        θ_r,
        Ksat,
        ΔS,
    )
    boundary_conditions = AtmosDrivenSnowBC(
        forcing.atmos,
        forcing.radiation;
        prognostic_land_components,
    )
    return SnowModel(; boundary_conditions, domain, parameters)
end

"""
    prognostic_vars(::SnowModel)

Returns the prognostic variable names of the snow model.

For this model, we track the snow water equivalent S in meters (liquid
water volume per ground area) and
the energy per unit ground area U [J/m^2] prognostically.
"""
prognostic_vars(m::SnowModel) =
    (:S, :S_l, :U, density_prog_vars(m.parameters.density)...)

"""
    density_prog_vars(::AbstractDensityModel)

A default method for adding prognostic variables to the snow model as required
by the density model choice.
"""
density_prog_vars(::AbstractDensityModel) = ()

"""
    prognostic_types(::SnowModel{FT})

Returns the prognostic variable types of the snow model;
both snow water equivalent and energy per unit area
are scalars.
"""
prognostic_types(m::SnowModel{FT}) where {FT} =
    (FT, FT, FT, density_prog_types(m.parameters.density)...)

"""
    density_prog_vars(::AbstractDensityModel)

A default method for specifying variable types of the prognostic variables required
by the density model choice, similar to `prognostic_types()`.
"""
density_prog_types(::AbstractDensityModel{FT}) where {FT} = ()

"""
    prognostic_domain_names(::SnowModel)

Returns the prognostic variable domain names of the snow model;
both snow water equivalent and energy per unit area
are modeling only as a function of (x,y), and not as a function
of depth. Therefore their domain name is ":surface".
"""
prognostic_domain_names(m::SnowModel) =
    (:surface, :surface, :surface, density_prog_names(m.parameters.density)...)

"""
    density_prog_vars(::AbstractDensityModel)

A default method for specifying variable domain names of the prognostic variables required
by the density model choice, similar to `prognostic_domain_names()`.
"""
density_prog_names(::AbstractDensityModel) = ()

"""
    auxiliary_vars(::SnowModel)

Returns the auxiliary variable names for the snow model. These
include the specific humidity at the surface of the snow `(`q_sfc`, unitless),
the mass fraction in liquid water (`q_l`, unitless),
the thermal conductivity (`κ`, W/m/K),
the bulk temperature (`T`, K), the surface temperature (`T_sfc`, K), the snow depth (`z_snow`, m),
the bulk snow density (`ρ_snow`, kg/m^3)
the SHF, LHF, and vapor flux (`turbulent_fluxes.shf`, etc),
the net radiation (`R_n, J/m^2/s)`, the energy flux in liquid water runoff
(`energy_runoff`, J/m^2/s), the water volume in runoff (`water_runoff`, m/s), and the total energy and water fluxes applied to the snowpack.

Since the snow can melt completely in one timestep, we clip the water and energy fluxes
such that SWE cannot become negative and U cannot become unphysical. The
clipped values are what are actually applied as boundary fluxes, and are stored in
`applied_` fluxes.
"""
auxiliary_vars(snow::SnowModel) = (
    :q_sfc,
    :q_l,
    :κ,
    :T,
    :T_sfc,
    :z_snow,
    :α_snow,
    :ρ_snow,
    :R_n,
    :phase_change_flux,
    :energy_runoff,
    :water_runoff,
    :liquid_water_flux,
    :total_energy_flux,
    :total_water_flux,
    :applied_energy_flux,
    :applied_water_flux,
    :snow_cover_fraction,
    boundary_vars(snow.boundary_conditions, ClimaLand.TopBoundary())...,
)

auxiliary_types(snow::SnowModel{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    boundary_var_types(
        snow,
        snow.boundary_conditions,
        ClimaLand.TopBoundary(),
    )...,
)

auxiliary_domain_names(snow::SnowModel) = (
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
    :surface,
    :surface,
    :surface,
    :surface,
    boundary_var_domain_names(
        snow.boundary_conditions,
        ClimaLand.TopBoundary(),
    )...,
)


ClimaLand.name(::SnowModel) = :snow

function ClimaLand.make_update_aux(model::SnowModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        parameters = model.parameters
        # The ordering is important here
        @. p.snow.q_l = liquid_mass_fraction(Y.snow.S, Y.snow.S_l)

        update_density_and_depth!(
            p.snow.ρ_snow,
            p.snow.z_snow,
            parameters.density,
            Y,
            p,
            parameters,
        )

        update_snow_albedo!(
            p.snow.α_snow,
            parameters.α_snow,
            Y,
            p,
            t,
            parameters.earth_param_set,
        ) # This could depend on ρ_snow

        @. p.snow.κ = snow_thermal_conductivity(p.snow.ρ_snow, parameters)

        @. p.snow.T =
            snow_bulk_temperature(Y.snow.U, Y.snow.S, p.snow.q_l, parameters)

        @. p.snow.T_sfc = snow_surface_temperature(p.snow.T)

        @. p.snow.water_runoff = compute_water_runoff(
            Y.snow.S,
            Y.snow.S_l,
            p.snow.T,
            p.snow.ρ_snow,
            p.snow.z_snow,
            parameters,
        )

        @. p.snow.energy_runoff =
            p.snow.water_runoff *
            volumetric_internal_energy_liq(p.snow.T, parameters)
        update_snow_cover_fraction!(
            p.snow.snow_cover_fraction,
            parameters.scf,
            Y,
            p,
            t,
            parameters.earth_param_set,
        ) # This depends on z_snow
    end
end

function ClimaLand.make_update_boundary_fluxes(model::SnowModel{FT}) where {FT}
    function update_boundary_fluxes!(p, Y, t)
        # First compute the boundary fluxes
        snow_boundary_fluxes!(model.boundary_conditions, model, Y, p, t)
        # Next, clip them in case the snow will melt in this timestep
        @. p.snow.applied_water_flux = clip_water_flux(
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
        # We now estimate the phase change flux: if the applied energy flux is such that T > T_f on the next step, use the residual after warming to T_f to melt snow
        # This estimate uses the current S and q_l.
        @. p.snow.phase_change_flux = phase_change_flux(
            Y.snow.U,
            Y.snow.S,
            p.snow.q_l,
            p.snow.applied_energy_flux,
            model.parameters,
        )
        @. p.snow.liquid_water_flux +=
            p.snow.phase_change_flux * p.snow.snow_cover_fraction
        @. p.snow.liquid_water_flux = clip_liquid_water_flux(
            Y.snow.S_l,
            Y.snow.S,
            p.snow.liquid_water_flux,
            p.snow.applied_water_flux,
            model.parameters.Δt,
        )
    end
end

function ClimaLand.make_compute_exp_tendency(model::SnowModel{FT}) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        # positive fluxes are TOWARDS atmos; negative fluxes increase quantity in snow
        @. dY.snow.S = -p.snow.applied_water_flux
        @. dY.snow.S_l = -p.snow.liquid_water_flux
        @. dY.snow.U = -p.snow.applied_energy_flux
        update_density_prog!(model.parameters.density, model, dY, Y, p)
    end
    return compute_exp_tendency!
end

"""
    clip_water_flux(S, total_water_flux, Δt)

A helper function which clips the total water flux so that
snow water equivalent S will not become negative in a timestep Δt.
"""
function clip_water_flux(S::FT, total_water_flux::FT, Δt::FT) where {FT}
    if S - total_water_flux * Δt < 0
        return S / Δt
    else
        return total_water_flux
    end
end

"""
    clip_liquid_water_flux(S_l::FT, S::FT, liquid_water_flux::FT, applied_water_flux::FT, Δt::FT) where {FT}

A helper function which clips the liquid water flux so that
snow liquid water S_l will not become negative or exceed S in a timestep Δt.
"""
function clip_liquid_water_flux(
    S_l::FT,
    S::FT,
    liquid_water_flux::FT,
    applied_water_flux::FT,
    Δt::FT,
) where {FT}
    predicted_S = S - applied_water_flux * Δt
    predicted_S_l = S_l - liquid_water_flux * Δt
    if predicted_S_l < 0
        return S_l / Δt
    elseif predicted_S_l > predicted_S
        return (S_l - predicted_S) / Δt
    else
        return liquid_water_flux
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

"""
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model::SnowModel,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value for
the total liquid water volume per unit ground area for the `SnowModel`.

This has already accounted for the area fraction of snow in the definition
of S; it also accounts for both liquid and frozen water present in the snow,
as the snow water equivalent is already the total liquid water volume present in the snow
if all the snow melted, per unit ground area.
"""
function ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::SnowModel,
    Y,
    p,
    t,
)
    surface_field .= Y.snow.S
    return nothing
end

"""
    ClimaLand.total_energy_per_area!(
        surface_field,
        model::SnowModel,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value for
the total energy per unit ground area for the `SnowModel`.

This has already accounted for the area fraction of snow in the definition
of S.
"""
function ClimaLand.total_energy_per_area!(
    surface_field,
    model::SnowModel,
    Y,
    p,
    t,
)
    surface_field .= Y.snow.U
    return nothing
end

end

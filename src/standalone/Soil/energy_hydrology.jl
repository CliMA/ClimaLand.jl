"""
    EnergyHydrologyParameters{FT <: AbstractFloat}

A parameter structure for the integrated soil water and energy
 equation system. In this simplest
form, we assume the conductivity and volumetric heat capacity
of the soil are constant.
$(DocStringExtensions.FIELDS)
"""
struct EnergyHydrologyParameters{
    FT <: AbstractFloat,
    C <: AbstractSoilHydrologyClosure,
    PSE,
}
    "The dry soil thermal conductivity, W/m/K"
    κ_dry::FT
    "The saturated thermal conductivity of frozen soil, W/m/K"
    κ_sat_frozen::FT
    "The saturated thermal conductivity of unfrozen soil, W/m/K"
    κ_sat_unfrozen::FT
    "The volumetric heat capacity of dry soil, J/m^3/K"
    ρc_ds::FT
    "The porosity of the soil (m^3/m^3)"
    ν::FT
    "The volumetric fraction of the soil solids in organic matter (m^3/m^3)"
    ν_ss_om::FT
    "The volumetric fraction of the soil solids in quartz (m^3/m^3)"
    ν_ss_quartz::FT
    "The volumetric fraction of the soil solids in gravel (m^3/m^3)"
    ν_ss_gravel::FT
    "The parameter α used in computing Kersten number, unitless"
    α::FT
    "The parameter β used in computing Kersten number, unitless"
    β::FT
    "The soil hydrology closure model: van Genucthen or Brooks and Corey"
    hydrology_cm::C
    "The saturated hydraulic conductivity (m/s)"
    K_sat::FT
    "The specific storativity (1/m)"
    S_s::FT
    "The residual water fraction (m^3/m^3"
    θ_r::FT
    "Ice impedance factor for the hydraulic conductivity"
    Ω::FT
    "Coefficient of viscosity factor for the hydraulic conductivity"
    γ::FT
    "Reference temperature for the viscosity factor"
    γT_ref::FT
    "Soil Albedo"
    albedo::FT
    "Soil Emissivity"
    emissivity::FT
    "Roughness length for momentum"
    z_0m::FT
    "Roughness length for scalars"
    z_0b::FT
    "Maximum dry soil layer thickness under evaporation (m)"
    d_ds::FT
    "Physical constants and clima-wide parameters"
    earth_param_set::PSE
end

function EnergyHydrologyParameters{FT}(;
    κ_dry::FT,
    κ_sat_frozen::FT,
    κ_sat_unfrozen::FT,
    ρc_ds::FT,
    ν::FT,
    ν_ss_om::FT,
    ν_ss_quartz::FT,
    ν_ss_gravel::FT,
    hydrology_cm::C,
    K_sat::FT,
    S_s::FT,
    θ_r::FT,
    albedo = FT(0.2),
    emissivity = FT(1),
    z_0m = FT(0.01),
    z_0b = FT(0.01),
    d_ds = FT(0.015),
    earth_param_set::PSE,
) where {FT <: AbstractFloat, PSE, C <: AbstractSoilHydrologyClosure}
    # These were determined in the Balland and Arp paper, from 2003.
    α = FT(0.24)
    β = FT(18.3)
    # Lundin paper
    Ω = FT(7)
    # Unclear where these are from (design doc)
    γ = FT(2.64e-2)
    γT_ref = FT(288)
    return EnergyHydrologyParameters{FT, C, PSE}(
        κ_dry,
        κ_sat_frozen,
        κ_sat_unfrozen,
        ρc_ds,
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        α,
        β,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        Ω,
        γ,
        γT_ref,
        albedo,
        emissivity,
        z_0m,
        z_0b,
        d_ds,
        earth_param_set,
    )
end

Base.broadcastable(ps::EnergyHydrologyParameters) = tuple(ps)

"""
    EnergyHydrology <: AbstractSoilModel

A model for simulating the flow of water and heat
in a porous medium by solving the Richardson-Richards equation
and the heat equation, including terms for phase change.

A variety of boundary condition types are supported, including
FluxBC, MoistureStateBC/TemperatureStateBC, 
FreeDrainage (only for the bottom of the domain),
and an AtmosDrivenFluxBC (under which radiative fluxes and
turbulent surface fluxes are computed and used as boundary conditions).
Please see the documentation for this boundary condition type for more
details.

$(DocStringExtensions.FIELDS)
"""
struct EnergyHydrology{FT, PS, D, BCRH, S} <: AbstractSoilModel{FT}
    "The parameter sets"
    parameters::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions for RRE and heat, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BCRH
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
    "A boolean flag which, when false, turns off the horizontal flow of water and heat"
    lateral_flow::Bool
end

"""
    EnergyHydrology{FT}(;
        parameters::PS
        domain::D,
        boundary_conditions::NamedTuple,
        sources::Tuple,
        lateral_flow::Bool = true
    ) where {FT, D, PS}

A constructor for a `EnergyHydrology` model, which sets the default value
of the `lateral_flow` flag to true.
"""
function EnergyHydrology{FT}(;
    parameters::EnergyHydrologyParameters{FT, PSE},
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
    lateral_flow::Bool = true,
) where {FT, D, PSE}
    top_bc = boundary_conditions.top
    if typeof(top_bc) <: AtmosDrivenFluxBC
        # If the top BC indicates atmospheric conditions are driving the model
        # add baseflow as a sink term
        subsurface_source = subsurface_runoff_source(top_bc.runoff)
        sources = append_source(subsurface_source, sources)
    end
    args = (parameters, domain, boundary_conditions, sources)
    EnergyHydrology{FT, typeof.(args)...}(args..., lateral_flow)
end

"""
    make_compute_exp_tendency(model::EnergyHydrology)

An extension of the function `make_compute_exp_tendency`, for the integrated
soil energy and heat equations, including phase change.

This function creates and returns a function which computes the entire
right hand side of the PDE for `Y.soil.ϑ_l, Y.soil.θ_i, Y.soil.ρe_int`,
and updates `dY.soil` in place with those values.
All of these quantities will be stepped explicitly.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_compute_exp_tendency(
    model::EnergyHydrology{FT},
) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        Δz_top, Δz_bottom = get_Δz(z)

        # Convert all boundary conditions to FluxBCs
        rre_top_flux_bc, heat_top_flux_bc = soil_boundary_fluxes(
            model.boundary_conditions.top,
            TopBoundary(),
            model,
            Δz_top,
            Y,
            p,
            t,
        )
        rre_bot_flux_bc, heat_bot_flux_bc = soil_boundary_fluxes(
            model.boundary_conditions.bottom,
            BottomBoundary(),
            model,
            Δz_bottom,
            Y,
            p,
            t,
        )

        interpc2f = Operators.InterpolateC2F()
        gradc2f = Operators.GradientC2F()

        # Without topography only
        # In Cartesian coordinates, W (z^) = Cov3 (z^)= Contra3 (n^ = z^)
        # In spherical coordinates, W (r^) = Cov3 (r^) = Contra3 (n^ = r^)
        # Passing WVector to gradient BC is passing a normal flux.


        # Richards-Richardson RHS
        divf2c_rre = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(rre_top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(rre_bot_flux_bc)),
        )
        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ϑ_l =
            -(divf2c_rre(-interpc2f(p.soil.K) * gradc2f(p.soil.ψ + z)))
        dY.soil.θ_i .= ClimaCore.Fields.zeros(FT, axes(Y.soil.θ_i))

        # Heat equation RHS
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(heat_top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(heat_bot_flux_bc)),
        )
        ρe_int_l = volumetric_internal_energy_liq.(p.soil.T, model.parameters)

        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ρe_int =
            -divf2c_heat(
                -interpc2f(p.soil.κ) * gradc2f(p.soil.T) -
                interpc2f(ρe_int_l * p.soil.K) * gradc2f(p.soil.ψ + z),
            )
        # Horizontal contributions
        horizontal_components!(
            dY,
            model.domain,
            Val(model.lateral_flow),
            model,
            p,
            z,
        )

        for src in model.sources

            ClimaLSM.source!(dY, src, Y, p, model)
        end
    end
    return compute_exp_tendency!
end

"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Union{HybridBox, SphericalShell},
                          lateral_flow::Val{true},
                          model::EnergyHydrology,
                          p::NamedTuple)

Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators for the `EnergyHydrology` model,
 in the case of a hybrid box or 
spherical shell domain with the model
`lateral_flag` set to true.

The horizontal contributions are
computed using the WeakDivergence and Gradient operators.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
    lateral_flow::Val{true},
    model::EnergyHydrology,
    p::NamedTuple,
    z::ClimaCore.Fields.Field,
)
    hdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    # The flux is already covariant, from hgrad, so no need to convert.
    @. dY.soil.ϑ_l += -hdiv(-p.soil.K * hgrad(p.soil.ψ + z))
    ρe_int_l = volumetric_internal_energy_liq.(p.soil.T, model.parameters)
    @. dY.soil.ρe_int +=
        -hdiv(
            -p.soil.κ * hgrad(p.soil.T) -
            p.soil.K * ρe_int_l * hgrad(p.soil.ψ + z),
        )
end

"""
    prognostic_vars(soil::EnergyHydrology)

A function which returns the names of the prognostic variables
of `EnergyHydrology`.
"""
ClimaLSM.prognostic_vars(soil::EnergyHydrology) = (:ϑ_l, :θ_i, :ρe_int)

"""
    prognostic_types(soil::EnergyHydrology{FT}) where {FT}

A function which returns the types of the prognostic variables
of `EnergyHydrology`.
"""
ClimaLSM.prognostic_types(soil::EnergyHydrology{FT}) where {FT} = (FT, FT, FT)
"""
    auxiliary_vars(soil::EnergyHydrology)

A function which returns the names of the auxiliary variables
of `EnergyHydrology`.
"""
ClimaLSM.auxiliary_vars(soil::EnergyHydrology) = (:K, :ψ, :θ_l, :T, :κ)

"""
    auxiliary_types(soil::EnergyHydrology{FT}) where {FT}

A function which returns the types of the auxiliary variables
of `EnergyHydrology`.
"""
ClimaLSM.auxiliary_types(soil::EnergyHydrology{FT}) where {FT} =
    (FT, FT, FT, FT, FT)

"""
    make_update_aux(model::EnergyHydrology)

An extension of the function `make_update_aux`, for the integrated
soil hydrology and energy model.

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::EnergyHydrology)
    function update_aux!(p, Y, t)
        (;
            ν,
            hydrology_cm,
            K_sat,
            S_s,
            θ_r,
            Ω,
            γ,
            γT_ref,
            κ_sat_frozen,
            κ_sat_unfrozen,
        ) = model.parameters

        @. p.soil.θ_l =
            volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i, θ_r)

        @. p.soil.κ = thermal_conductivity(
            model.parameters.κ_dry,
            kersten_number(
                Y.soil.θ_i,
                relative_saturation(p.soil.θ_l, Y.soil.θ_i, ν),
                model.parameters,
            ),
            κ_sat(p.soil.θ_l, Y.soil.θ_i, κ_sat_unfrozen, κ_sat_frozen),
        )

        @. p.soil.T = temperature_from_ρe_int(
            Y.soil.ρe_int,
            Y.soil.θ_i,
            volumetric_heat_capacity(p.soil.θ_l, Y.soil.θ_i, model.parameters),
            model.parameters,
        )

        @. p.soil.K =
            impedance_factor(Y.soil.θ_i / (p.soil.θ_l + Y.soil.θ_i - θ_r), Ω) *
            viscosity_factor(p.soil.T, γ, γT_ref) *
            hydraulic_conductivity(
                hydrology_cm,
                K_sat,
                effective_saturation(ν, Y.soil.ϑ_l, θ_r),
            )
        @. p.soil.ψ =
            pressure_head(hydrology_cm, θ_r, Y.soil.ϑ_l, ν - Y.soil.θ_i, S_s)
    end
    return update_aux!
end

"""
    PhaseChange{FT} <: AbstractSoilSource{FT}

PhaseChange source type.
"""
struct PhaseChange{FT} <: AbstractSoilSource{FT}
    Δz::FT
end


"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::PhaseChange{FT},
             Y::ClimaCore.Fields.FieldVector,
             p::NamedTuple,
             model
             )

Computes the source terms for phase change.

"""
function ClimaLSM.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::PhaseChange{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
) where {FT}
    params = model.parameters
    (; ν, earth_param_set) = params
    _ρ_l = FT(LSMP.ρ_cloud_liq(earth_param_set))
    _ρ_i = FT(LSMP.ρ_cloud_ice(earth_param_set))
    ρc = volumetric_heat_capacity.(p.soil.θ_l, Y.soil.θ_i, params)
    τ = thermal_time.(ρc, src.Δz, p.soil.κ)

    liquid_source =
        phase_change_source.(p.soil.θ_l, Y.soil.θ_i, p.soil.T, τ, params)
    @. dY.soil.ϑ_l += -liquid_source
    @. dY.soil.θ_i += (_ρ_l / _ρ_i) * liquid_source
end



## The functions below are required to be defined
## for use with the `AtmosDrivenFluxBC` upper
## boundary conditions. They return the
## surface conditions necessary for computing
## radiative, sensible and latent heat fluxes
## as well as evaporation.


"""
    ClimaLSM.surface_temperature(
        model::EnergyHydrology{FT},
        Y,
        p,
        t,
    ) where {FT}

Returns the surface temperature field of the
`EnergyHydrology` soil model.

The assumption is that the soil surface temperature
is the same as the temperature at the center of the
first soil layer.
"""
function ClimaLSM.surface_temperature(
    model::EnergyHydrology{FT},
    Y,
    p,
    t,
) where {FT}
    return ClimaLSM.Domains.top_center_to_surface(p.soil.T)
end

"""
    ClimaLSM.surface_emissivity(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface emissivity field of the
`EnergyHydrology` soil model.
"""
function ClimaLSM.surface_emissivity(
    model::EnergyHydrology{FT},
    Y,
    p,
) where {FT}
    return model.parameters.emissivity
end

"""
    ClimaLSM.surface_albedo(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface albedo field of the
`EnergyHydrology` soil model.
"""
function ClimaLSM.surface_albedo(model::EnergyHydrology{FT}, Y, p) where {FT}
    return model.parameters.albedo
end

"""
    ClimaLSM.surface_air_density(
        atmos::PrescribedAtmosphere{FT},
        model::EnergyHydrology{FT},
        Y,
        p,
        t,
        T_sfc
    ) where {FT}

Returns the surface air density field of the
`EnergyHydrology` soil model for the
`PrescribedAtmosphere` case.

This assumes the ideal gas law and hydrostatic
balance to estimate the air density at the surface
from the values of surface temperature and the atmospheric
thermodynamic state,
because the surface air density is not a prognostic
variable of the soil model.
"""
function ClimaLSM.surface_air_density(
    atmos::PrescribedAtmosphere{FT},
    model::EnergyHydrology{FT},
    Y,
    p,
    t,
    T_sfc,
) where {FT}

    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    ts_in = construct_atmos_ts(atmos, t, thermo_params)
    return compute_ρ_sfc.(thermo_params, Ref(ts_in), T_sfc)
end

"""
    ClimaLSM.surface_specific_humidity(
        model::EnergyHydrology{FT},
        Y,
        p,
        T_sfc,
        ρ_sfc
    ) where {FT}

Returns the surface specific humidity field of the
`EnergyHydrology` soil model.

This models the surface specific humidity as
the saturated value multiplied by
the factor `exp(ψ_sfc g M_w/(RT_sfc))` in accordance
with the Clausius-Clapeyron equation,
where `ψ_sfc` is the matric potential at the surface,
`T_sfc` the surface temperature, `g` the gravitational
acceleration on the surface of the Earth, `M_w` the molar
mass of water, and `R` the universal gas constant.
"""
function ClimaLSM.surface_specific_humidity(
    model::EnergyHydrology{FT},
    Y,
    p,
    T_sfc,
    ρ_sfc,
) where {FT}
    g = LSMP.grav(model.parameters.earth_param_set)
    R = LSMP.gas_constant(model.parameters.earth_param_set)
    M_w = LSMP.molar_mass_water(model.parameters.earth_param_set)
    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    ψ_sfc = ClimaLSM.Domains.top_center_to_surface(p.soil.ψ)
    q_sat =
        Thermodynamics.q_vap_saturation_generic.(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Liquid(),
        )
    return @. (q_sat * exp(g * ψ_sfc * M_w / (R * T_sfc)))

end

"""
    ClimaLSM.surface_height(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface height of the `EnergyHydrology` model.
"""
function ClimaLSM.surface_height(model::EnergyHydrology{FT}, Y, p) where {FT}
    face_space = ClimaLSM.Domains.obtain_face_space(model.domain.space)
    N = ClimaCore.Spaces.nlevels(face_space)
    surface_space = ClimaLSM.Domains.obtain_surface_space(model.domain.space)
    z_sfc = ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(
            ClimaCore.Fields.level(
                ClimaCore.Fields.coordinate_field(face_space).z,
                ClimaCore.Utilities.PlusHalf(N - 1),
            ),
        ),
        surface_space,
    )
    return z_sfc
end

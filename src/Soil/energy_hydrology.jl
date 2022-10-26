"""
    EnergyHydrologyParameters{FT <: AbstractFloat}

A parameter structure for the integrated soil water and energy
 equation system. In this simplest
form, we assume the conductivity and volumetric heat capacity
of the soil are constant.
$(DocStringExtensions.FIELDS)
"""
struct EnergyHydrologyParameters{FT <: AbstractFloat, PSE}
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
    "The van Genuchten parameter α (1/m)"
    vg_α::FT
    "The van Genuchten parameter n"
    vg_n::FT
    "The van Genuchten parameter m"
    vg_m::FT
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
    "Physical constants and clima-wide parameters"
    earth_param_set::PSE
end

function EnergyHydrologyParameters(;
    κ_dry::FT,
    κ_sat_frozen::FT,
    κ_sat_unfrozen::FT,
    ρc_ds::FT,
    ν::FT,
    ν_ss_om::FT,
    ν_ss_quartz::FT,
    ν_ss_gravel::FT,
    vg_α::FT,
    vg_n::FT,
    K_sat::FT,
    S_s::FT,
    θ_r::FT,
    earth_param_set::PSE,
) where {FT, PSE}
    vg_m = FT(1.0 - 1.0 / vg_n)
    # These were determined in the Balland and Arp paper, from 2003.
    α = FT(0.24)
    β = FT(18.3)
    # Lundin paper
    Ω = FT(7)
    # Unclear where these are from (design doc)
    γ = FT(2.64e-2)
    γT_ref = FT(288)
    return EnergyHydrologyParameters{FT, PSE}(
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
        vg_α,
        vg_n,
        vg_m,
        K_sat,
        S_s,
        θ_r,
        Ω,
        γ,
        γT_ref,
        earth_param_set,
    )
end




"""
    EnergyHydrology <: AbstractSoilModel

A model for simulating the flow of water and heat 
in a porous medium by solving the Richardson-Richards equation
and the heat equation, including terms for phase change.

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
end

"""
    EnergyHydrology{FT}(;
        parameters::PS
        domain::D,
        rre_heat_boundary_conditions::RREHeatBoundaryConditions{FT},
        sources::Tuple,
    ) where {FT, D, PS}

A constructor for a `EnergyHydrology` model.
"""
function EnergyHydrology{FT}(;
    parameters::EnergyHydrologyParameters{FT, PSE},
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
) where {FT, D, PSE}
    args = (parameters, domain, boundary_conditions, sources)
    EnergyHydrology{FT, typeof.(args)...}(args...)
end

"""
    make_rhs(model::EnergyHydrology)

An extension of the function `make_rhs`, for the integrated soil
energy and heat equations, including phase change.

This function creates and returns a function which computes the entire
right hand side of the PDE for `Y.soil.ϑ_l, Y.soil.θ_i, Y.soil.ρe_int`, 
and updates `dY.soil` in place with those values.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_rhs(model::EnergyHydrology{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        Δz_top, Δz_bottom = get_Δz(z)

        # Convert all boundary conditions to FluxBCs
        rre_top_flux_bc = boundary_flux(
            model.boundary_conditions.water.top,
            TopBoundary(),
            Δz_top,
            p,
            t,
            model.parameters,
        )
        rre_bot_flux_bc = boundary_flux(
            model.boundary_conditions.water.bottom,
            BottomBoundary(),
            Δz_bottom,
            p,
            t,
            model.parameters,
        )

        heat_top_flux_bc = boundary_flux(
            model.boundary_conditions.heat.top,
            TopBoundary(),
            Δz_top,
            p,
            t,
        )
        heat_bot_flux_bc = boundary_flux(
            model.boundary_conditions.heat.bottom,
            BottomBoundary(),
            Δz_bottom,
            p,
            t,
        )

        interpc2f = Operators.InterpolateC2F()
        gradc2f = Operators.GradientC2F()

        # Without topography only
        # In Cartesian coordinates, W (z^) = Cov3 (z^)= Contra3 (n^ = z^)
        # In spherical coordinates, W (r^) = Cov3 (r^) = Contra3 (n^ = r^)
        # Passing WVector to gradient BC is passing a normal flux.


        # Richards-Richardson RHS:
        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
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
        ρe_int_l =
            volumetric_internal_energy_liq.(p.soil.T, Ref(model.parameters))

        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ρe_int =
            -divf2c_heat(
                -interpc2f(p.soil.κ) * gradc2f(p.soil.T) -
                interpc2f(ρe_int_l * p.soil.K) * gradc2f(p.soil.ψ + z),
            )
        # Horizontal contributions
        horizontal_components!(dY, model.domain, model, p, z)

        for src in model.sources

            ClimaLSM.source!(dY, src, Y, p, model.parameters)
        end

        # This has to come last
        dss!(dY, model.domain)
    end
    return rhs!
end

"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::HybridBox,
                          model::EnergyHydrology,
                          p::ClimaCore.Fields.FieldVector)

Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a hybrid box domain, the horizontal contributions are
computed using the WeakDivergence and Gradient operators.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
    model::EnergyHydrology,
    p::ClimaCore.Fields.FieldVector,
    z::ClimaCore.Fields.Field,
)
    hdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    # The flux is already covariant, from hgrad, so no need to convert.
    @. dY.soil.ϑ_l += -hdiv(-p.soil.K * hgrad(p.soil.ψ + z))
    ρe_int_l = volumetric_internal_energy_liq.(p.soil.T, Ref(model.parameters))
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
ClimaLSM.auxiliary_vars(soil::EnergyHydrology) = (:K, :ψ, :T, :κ)

"""
    auxiliary_types(soil::EnergyHydrology{FT}) where {FT}

A function which returns the types of the auxiliary variables
of `EnergyHydrology`.
"""
ClimaLSM.auxiliary_types(soil::EnergyHydrology{FT}) where {FT} =
    (FT, FT, FT, FT)

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
        @unpack ν,
        vg_α,
        vg_n,
        vg_m,
        K_sat,
        S_s,
        θ_r,
        Ω,
        γ,
        γT_ref,
        κ_sat_frozen,
        κ_sat_unfrozen = model.parameters

        θ_l = volumetric_liquid_fraction.(Y.soil.ϑ_l, ν .- Y.soil.θ_i)
        p.soil.κ =
            thermal_conductivity.(
                model.parameters.κ_dry,
                kersten_number.(
                    Y.soil.θ_i,
                    relative_saturation.(θ_l, Y.soil.θ_i, ν),
                    Ref(model.parameters),
                ),
                κ_sat.(θ_l, Y.soil.θ_i, κ_sat_unfrozen, κ_sat_frozen),
            )

        p.soil.T =
            temperature_from_ρe_int.(
                Y.soil.ρe_int,
                Y.soil.θ_i,
                volumetric_heat_capacity.(
                    θ_l,
                    Y.soil.θ_i,
                    Ref(model.parameters),
                ),
                Ref(model.parameters),
            )

        @. p.soil.K =
            impedance_factor(Y.soil.θ_i / (θ_l + Y.soil.θ_i - θ_r), Ω) *
            viscosity_factor(p.soil.T, γ, γT_ref) *
            hydraulic_conductivity(
                K_sat,
                vg_m,
                effective_saturation(ν, Y.soil.ϑ_l, θ_r),
            )
        @. p.soil.ψ = pressure_head(
            vg_α,
            vg_n,
            vg_m,
            θ_r,
            Y.soil.ϑ_l,
            ν - Y.soil.θ_i,
            S_s,
        )

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
             p::ClimaCore.Fields.FieldVector
             )

Computes the source terms for phase change.

"""
function ClimaLSM.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::PhaseChange{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    params,
) where {FT}
    @unpack ν, earth_param_set = params
    _ρ_l = FT(LSMP.ρ_cloud_liq(earth_param_set))
    _ρ_i = FT(LSMP.ρ_cloud_ice(earth_param_set))
    θ_l = @. volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i)
    ρc = volumetric_heat_capacity.(θ_l, Y.soil.θ_i, Ref(params))
    τ = thermal_time.(ρc, src.Δz, p.soil.κ)

    liquid_source =
        phase_change_source.(θ_l, Y.soil.θ_i, p.soil.T, τ, Ref(params))
    @. dY.soil.ϑ_l += -liquid_source
    @. dY.soil.θ_i += (_ρ_l / _ρ_i) * liquid_source
end

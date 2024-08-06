"""
    EnergyHydrologyParameters{
            FT <: AbstractFloat,
            F <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
            SF <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
            C,
            PSE,
        }

A parameter structure for the integrated soil water and energy
 equation system.

Note that we require two different parameter types F and SF; these
are for parameters that are defined on the surface only
and those defined in the interior of the soil domain:

- Surface parameters: albedo in each wavelength band (SF)
- Scalar parameters: emissivity, α, β, γ, γT_ref, Ω,
 roughness lengths z_0, d_ds ) (FT)
- Parameters defined in the interior: all else (F)

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct EnergyHydrologyParameters{
    FT <: AbstractFloat,
    F <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
    SF <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
    C,
    PSE,
}
    "The dry soil thermal conductivity, W/m/K"
    κ_dry::F
    "The saturated thermal conductivity of frozen soil, W/m/K"
    κ_sat_frozen::F
    "The saturated thermal conductivity of unfrozen soil, W/m/K"
    κ_sat_unfrozen::F
    "The volumetric heat capacity of dry soil, J/m^3/K (per volume dry soil, not per volume soil solids)"
    ρc_ds::F
    "The porosity of the soil (m^3/m^3)"
    ν::F
    "The volumetric fraction of the soil solids in organic matter (m^3/m^3)"
    ν_ss_om::F
    "The volumetric fraction of the soil solids in quartz (m^3/m^3)"
    ν_ss_quartz::F
    "The volumetric fraction of the soil solids in gravel (m^3/m^3)"
    ν_ss_gravel::F
    "The parameter α used in computing Kersten number, unitless"
    α::FT
    "The parameter β used in computing Kersten number, unitless"
    β::FT
    "The soil hydrology closure model: van Genuchten or Brooks and Corey"
    hydrology_cm::C
    "The saturated hydraulic conductivity (m/s)"
    K_sat::F
    "The specific storativity (1/m)"
    S_s::F
    "The residual water fraction (m^3/m^3"
    θ_r::F
    "Ice impedance factor for the hydraulic conductivity"
    Ω::FT
    "Coefficient of viscosity factor for the hydraulic conductivity"
    γ::FT
    "Reference temperature for the viscosity factor"
    γT_ref::FT
    "Soil PAR Albedo"
    PAR_albedo::SF
    "Soil NIR Albedo"
    NIR_albedo::SF
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
        # add baseflow as a sink term, add sublimation as a sink term
        subl_source = sublimation_source(top_bc)
        subsurface_source = subsurface_runoff_source(top_bc.runoff)
        sources = append_source(subsurface_source, sources)
        sources = append_source(subl_source, sources)
    end
    args = (parameters, domain, boundary_conditions, sources)
    EnergyHydrology{FT, typeof.(args)...}(args..., lateral_flow)
end

function make_update_boundary_fluxes(model::EnergyHydrology)
    function update_boundary_fluxes!(p, Y, t)
        Δz_top = model.domain.fields.Δz_top
        Δz_bottom = model.domain.fields.Δz_bottom
        soil_boundary_fluxes!(
            model.boundary_conditions.top,
            ClimaLand.TopBoundary(),
            model,
            Δz_top,
            Y,
            p,
            t,
        )

        soil_boundary_fluxes!(
            model.boundary_conditions.bottom,
            ClimaLand.BottomBoundary(),
            model,
            Δz_bottom,
            Y,
            p,
            t,
        )
    end
    return update_boundary_fluxes!
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
function ClimaLand.make_compute_exp_tendency(
    model::EnergyHydrology{FT},
) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        z = model.domain.fields.z

        # Don't update the prognostic variables we're stepping implicitly
        dY.soil.ϑ_l .= 0
        dY.soil.ρe_int .= 0

        # Note that soil ice content is only updated via source terms,
        # which are plus-added to dY.soil.θ_i, so we must zero it out here first.
        dY.soil.θ_i .= 0

        # Horizontal contributions
        horizontal_components!(
            dY,
            model.domain,
            Val(model.lateral_flow),
            model,
            p,
            z,
        )
        # Source terms
        for src in model.sources
            ClimaLand.source!(dY, src, Y, p, model)
        end
    end
    return compute_exp_tendency!
end

"""
    make_compute_imp_tendency(model::EnergyHydrology)

An extension of the function `make_compute_imp_tendency`, for the integrated
soil energy and heat equations, including phase change.

This version of this function computes the right hand side of the PDE for
`Y.soil.ϑ_l`, which is the only quantity we currently step implicitly.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_compute_imp_tendency(
    model::EnergyHydrology{FT},
) where {FT}
    function compute_imp_tendency!(dY, Y, p, t)
        z = model.domain.fields.z
        rre_top_flux_bc = p.soil.top_bc.water
        rre_bottom_flux_bc = p.soil.bottom_bc.water
        heat_top_flux_bc = p.soil.top_bc.heat
        heat_bottom_flux_bc = p.soil.bottom_bc.heat

        interpc2f = Operators.InterpolateC2F()
        gradc2f = Operators.GradientC2F()

        # Without topography only
        # In Cartesian coordinates, W (z^) = Cov3 (z^)= Contra3 (n^ = z^)
        # In spherical coordinates, W (r^) = Cov3 (r^) = Contra3 (n^ = r^)
        # Passing WVector to gradient BC is passing a normal flux.

        # Richards-Richardson RHS
        divf2c_rre = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(rre_top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(rre_bottom_flux_bc)),
        )
        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ϑ_l =
            -(divf2c_rre(-interpc2f(p.soil.K) * gradc2f(p.soil.ψ + z)))

        # Heat equation RHS
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(heat_top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(heat_bottom_flux_bc)),
        )

        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ρe_int =
            -divf2c_heat(
                -interpc2f(p.soil.κ) * gradc2f(p.soil.T) -
                interpc2f(
                    volumetric_internal_energy_liq(
                        p.soil.T,
                        model.parameters.earth_param_set,
                    ) * p.soil.K,
                ) * gradc2f(p.soil.ψ + z),
            )

        # Don't update the prognostic variables we're stepping explicitly
        @. dY.soil.θ_i = 0
    end
    return compute_imp_tendency!
end

"""
    ClimaLand.make_compute_jacobian(model::EnergyHydrology{FT}) where {FT}

Creates and returns the compute_jacobian! function for the EnergyHydrology model.
This updates the contribution for the soil liquid water content only.

Using this Jacobian with a backwards Euler timestepper is equivalent
to using the modified Picard scheme of Celia et al. (1990).
"""
function ClimaLand.make_compute_jacobian(model::EnergyHydrology{FT}) where {FT}
    function compute_jacobian!(jacobian::ImplicitEquationJacobian, Y, p, dtγ, t)
        (; matrix) = jacobian
        (; ν, hydrology_cm, S_s, θ_r, ρc_ds, earth_param_set) = model.parameters

        # Create divergence operator
        divf2c_op = Operators.DivergenceF2C()
        divf2c_matrix = MatrixFields.operator_matrix(divf2c_op)
        # Create gradient operator, and set gradient at boundaries to 0
        gradc2f_op = Operators.GradientC2F(
            top = Operators.SetGradient(Geometry.WVector(FT(0))),
            bottom = Operators.SetGradient(Geometry.WVector(FT(0))),
        )
        gradc2f_matrix = MatrixFields.operator_matrix(gradc2f_op)
        # Create interpolation operator
        interpc2f_op = Operators.InterpolateC2F(
            bottom = Operators.Extrapolate(),
            top = Operators.Extrapolate(),
        )

        # The derivative of the residual with respect to the prognostic variable
        ∂ϑres∂ϑ = matrix[@name(soil.ϑ_l), @name(soil.ϑ_l)]
        ∂ρeres∂ρe = matrix[@name(soil.ρe_int), @name(soil.ρe_int)]
        # If the top BC is a `MoistureStateBC`, add the term from the top BC
        #  flux before applying divergence
        if haskey(p.soil, :dfluxBCdY)
            dfluxBCdY = p.soil.dfluxBCdY
            topBC_op = Operators.SetBoundaryOperator(
                top = Operators.SetValue(dfluxBCdY),
                bottom = Operators.SetValue(
                    Geometry.Covariant3Vector(zero(FT)),
                ),
            )
            # Add term from top boundary condition before applying divergence
            # Note: need to pass 3D field on faces to `topBC_op`. Interpolating `K` to faces
            #  for this is inefficient - we should find a better solution.
            @. ∂ϑres∂ϑ =
                -dtγ * (
                    divf2c_matrix() ⋅ (
                        MatrixFields.DiagonalMatrixRow(
                            interpc2f_op(-p.soil.K),
                        ) ⋅ gradc2f_matrix() ⋅ MatrixFields.DiagonalMatrixRow(
                            ClimaLand.Soil.dψdϑ(
                                hydrology_cm,
                                Y.soil.ϑ_l,
                                ν,
                                θ_r,
                                S_s,
                            ),
                        ) + MatrixFields.LowerDiagonalMatrixRow(
                            topBC_op(
                                Geometry.Covariant3Vector(
                                    zero(interpc2f_op(p.soil.K)),
                                ),
                            ),
                        )
                    )
                ) - (I,)
        else
            @. ∂ϑres∂ϑ =
                -dtγ * (
                    divf2c_matrix() ⋅
                    MatrixFields.DiagonalMatrixRow(interpc2f_op(-p.soil.K)) ⋅
                    gradc2f_matrix() ⋅ MatrixFields.DiagonalMatrixRow(
                        ClimaLand.Soil.dψdϑ(
                            hydrology_cm,
                            Y.soil.ϑ_l,
                            ν,
                            θ_r,
                            S_s,
                        ),
                    )
                ) - (I,)
        end
        @. ∂ρeres∂ρe =
            -dtγ * (
                divf2c_matrix() ⋅
                MatrixFields.DiagonalMatrixRow(interpc2f_op(-p.soil.κ)) ⋅
                gradc2f_matrix() ⋅ MatrixFields.DiagonalMatrixRow(
                    1 / ClimaLand.Soil.volumetric_heat_capacity(
                        p.soil.θ_l,
                        Y.soil.θ_i,
                        ρc_ds,
                        earth_param_set,
                    ),
                )
            ) - (I,)
    end
    return compute_jacobian!
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
    @. dY.soil.ρe_int +=
        -hdiv(
            -p.soil.κ * hgrad(p.soil.T) -
            p.soil.K *
            volumetric_internal_energy_liq(
                p.soil.T,
                model.parameters.earth_param_set,
            ) *
            hgrad(p.soil.ψ + z),
        )
end

"""
    prognostic_vars(soil::EnergyHydrology)

A function which returns the names of the prognostic variables
of `EnergyHydrology`.
"""
ClimaLand.prognostic_vars(soil::EnergyHydrology) = (:ϑ_l, :θ_i, :ρe_int)

"""
    prognostic_types(soil::EnergyHydrology{FT}) where {FT}

A function which returns the types of the prognostic variables
of `EnergyHydrology`.
"""
ClimaLand.prognostic_types(soil::EnergyHydrology{FT}) where {FT} = (FT, FT, FT)

ClimaLand.prognostic_domain_names(soil::EnergyHydrology) =
    (:subsurface, :subsurface, :subsurface)
"""
    auxiliary_vars(soil::EnergyHydrology)

A function which returns the names of the auxiliary variables
of `EnergyHydrology`.
"""
ClimaLand.auxiliary_vars(soil::EnergyHydrology) = (
    :K,
    :ψ,
    :θ_l,
    :T,
    :κ,
    boundary_vars(soil.boundary_conditions.top, ClimaLand.TopBoundary())...,
    boundary_vars(
        soil.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)

"""
    auxiliary_types(soil::EnergyHydrology{FT}) where {FT}

A function which returns the types of the auxiliary variables
of `EnergyHydrology`.
"""
ClimaLand.auxiliary_types(soil::EnergyHydrology{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    boundary_var_types(
        soil,
        soil.boundary_conditions.top,
        ClimaLand.TopBoundary(),
    )...,
    boundary_var_types(
        soil,
        soil.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)

ClimaLand.auxiliary_domain_names(soil::EnergyHydrology) = (
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    boundary_var_domain_names(
        soil.boundary_conditions.top,
        ClimaLand.TopBoundary(),
    )...,
    boundary_var_domain_names(
        soil.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)
"""
    make_update_aux(model::EnergyHydrology)

An extension of the function `make_update_aux`, for the integrated
soil hydrology and energy model.

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_update_aux(model::EnergyHydrology)
    function update_aux!(p, Y, t)
        (;
            ν,
            hydrology_cm,
            K_sat,
            S_s,
            θ_r,
            Ω,
            γ,
            α,
            β,
            ν_ss_om,
            ν_ss_gravel,
            ν_ss_quartz,
            γT_ref,
            κ_sat_frozen,
            κ_sat_unfrozen,
            ρc_ds,
            earth_param_set,
        ) = model.parameters

        @. p.soil.θ_l =
            volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i, θ_r)

        @. p.soil.κ = thermal_conductivity(
            model.parameters.κ_dry,
            kersten_number(
                Y.soil.θ_i,
                relative_saturation(p.soil.θ_l, Y.soil.θ_i, ν),
                α,
                β,
                ν_ss_om,
                ν_ss_quartz,
                ν_ss_gravel,
            ),
            κ_sat(p.soil.θ_l, Y.soil.θ_i, κ_sat_unfrozen, κ_sat_frozen),
        )

        @. p.soil.T = temperature_from_ρe_int(
            Y.soil.ρe_int,
            Y.soil.θ_i,
            volumetric_heat_capacity(
                p.soil.θ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            ),
            earth_param_set,
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
struct PhaseChange{FT} <: AbstractSoilSource{FT} end


"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::PhaseChange{FT},
             Y::ClimaCore.Fields.FieldVector,
             p::NamedTuple,
             model
             )

Computes the source terms for phase change.

"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::PhaseChange{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
) where {FT}
    params = model.parameters
    (; ν, ρc_ds, θ_r, hydrology_cm, earth_param_set) = params
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
    Δz_top = model.domain.fields.Δz_top # center face distance
    @. dY.soil.ϑ_l +=
        -phase_change_source(
            p.soil.θ_l,
            Y.soil.θ_i,
            p.soil.T,
            thermal_time(
                volumetric_heat_capacity(
                    p.soil.θ_l,
                    Y.soil.θ_i,
                    ρc_ds,
                    earth_param_set,
                ),
                2 * Δz_top, # the factor of 2 appears to get the face-face/layer thickness, Δz_top is center-face distance
                p.soil.κ,
            ),
            ν,
            θ_r,
            hydrology_cm,
            earth_param_set,
        )
    @. dY.soil.θ_i +=
        (_ρ_l / _ρ_i) * phase_change_source(
            p.soil.θ_l,
            Y.soil.θ_i,
            p.soil.T,
            thermal_time(
                volumetric_heat_capacity(
                    p.soil.θ_l,
                    Y.soil.θ_i,
                    ρc_ds,
                    earth_param_set,
                ),
                2 * Δz_top, #the factor of 2 appears to get the face-face/layer thickness, Δz_top is center-face distance
                p.soil.κ,
            ),
            ν,
            θ_r,
            hydrology_cm,
            earth_param_set,
        )
end


"""
    SoilSublimation{FT} <: AbstractSoilSource{FT}

Soil Sublimation source type. Used to defined a method
of `ClimaLand.source!` for soil sublimation.
"""
struct SoilSublimation{FT} <: AbstractSoilSource{FT} end

"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::SoilSublimation{FT},
             Y::ClimaCore.Fields.FieldVector,
             p::NamedTuple,
             model
             )

Updates dY.soil.θ_i in place with a term due to sublimation; this only affects
the surface layer of soil.

"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::SoilSublimation{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
) where {FT}
    _ρ_i = FT(LP.ρ_cloud_ice(model.parameters.earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(model.parameters.earth_param_set))
    z = model.domain.fields.z
    Δz_top = model.domain.fields.Δz_top # this returns the center-face distance, not layer thickness
    @. dY.soil.θ_i +=
        -p.soil.turbulent_fluxes.vapor_flux * p.soil.ice_frac * _ρ_l / _ρ_i *
        heaviside(z + 2 * Δz_top) # only apply to top layer, recall that z is negative
end

## The functions below are required to be defined
## for use with the `AtmosDrivenFluxBC` upper
## boundary conditions. They return the
## surface conditions necessary for computing
## radiative, sensible and latent heat fluxes
## as well as evaporation.

"""
    ClimaLand.surface_temperature(
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
function ClimaLand.surface_temperature(
    model::EnergyHydrology{FT},
    Y,
    p,
    t,
) where {FT}
    return ClimaLand.Domains.top_center_to_surface(p.soil.T)
end

"""
    ClimaLand.surface_resistance(
        model::EnergyHydrology{FT},
        Y,
        p,
        t,
    ) where {FT}

Returns the surface resistance field of the
`EnergyHydrology` soil model.
"""
function ClimaLand.surface_resistance(
    model::EnergyHydrology{FT},
    Y,
    p,
    t,
) where {FT}
    (; ν, θ_r, d_ds, earth_param_set, hydrology_cm) = model.parameters
    θ_l_sfc = p.soil.sfc_scratch
    ClimaLand.Domains.linear_interpolation_to_surface!(
        θ_l_sfc,
        p.soil.θ_l,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    # These are non-allocating
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    hydrology_cm_sfc = ClimaLand.Domains.top_center_to_surface(hydrology_cm)
    ν_sfc = ClimaLand.Domains.top_center_to_surface(ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(θ_r)
    return ClimaLand.Soil.soil_resistance.(
        θ_l_sfc,
        θ_i_sfc,
        hydrology_cm_sfc,
        ν_sfc,
        θ_r_sfc,
        d_ds,
        earth_param_set,
    )
end

"""
    ClimaLand.surface_emissivity(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface emissivity field of the
`EnergyHydrology` soil model.
"""
function ClimaLand.surface_emissivity(
    model::EnergyHydrology{FT},
    Y,
    p,
) where {FT}
    return model.parameters.emissivity
end

"""
    ClimaLand.surface_albedo(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface albedo field of the
`EnergyHydrology` soil model.
"""
function ClimaLand.surface_albedo(model::EnergyHydrology{FT}, Y, p) where {FT}
    return @. FT(
        0.5 * model.parameters.PAR_albedo + 0.5 * model.parameters.NIR_albedo,
    )
end

"""
    ClimaLand.surface_specific_humidity(
        model::EnergyHydrology{FT},
        Y,
        p,
        T_sfc,
        ρ_sfc
    ) where {FT}

Returns the surface specific humidity field of the
`EnergyHydrology` soil model.

This models the specific humidity over the soil liquid water as
the saturated value multiplied by
the factor `exp(ψ_sfc g M_w/(RT_sfc))` in accordance
with the Clausius-Clapeyron equation,
where `ψ_sfc` is the matric potential at the surface,
`T_sfc` the surface temperature, `g` the gravitational
acceleration on the surface of the Earth, `M_w` the molar
mass of water, and `R` the universal gas constant.

Over the soil ice, the specific humidity is the saturated value.

The total surface specific humidity of the soil is approximated by
q = q_over_ice * f + q_over_water * (1-f), where `f` is given by
the function `ice_fraction`.
"""
function ClimaLand.surface_specific_humidity(
    model::EnergyHydrology{FT},
    Y,
    p,
    T_sfc,
    ρ_sfc,
) where {FT}
    g = LP.grav(model.parameters.earth_param_set)
    R = LP.gas_constant(model.parameters.earth_param_set)
    M_w = LP.molar_mass_water(model.parameters.earth_param_set)
    thermo_params =
        LP.thermodynamic_parameters(model.parameters.earth_param_set)
    ψ_sfc = p.soil.sfc_scratch
    ClimaLand.Domains.linear_interpolation_to_surface!(
        ψ_sfc,
        p.soil.ψ,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    q_over_ice =
        Thermodynamics.q_vap_saturation_generic.(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Ice(),
        )
    q_over_liq =
        Thermodynamics.q_vap_saturation_generic.(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Liquid(),
        ) .* (@. exp(g * ψ_sfc * M_w / (R * T_sfc)))
    ice_frac = p.soil.ice_frac
    ν_sfc = ClimaLand.Domains.top_center_to_surface(model.parameters.ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(model.parameters.θ_r)
    S_c_sfc = ClimaLand.Domains.top_center_to_surface(
        model.parameters.hydrology_cm.S_c,
    )
    θ_l_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    @. ice_frac = ice_fraction(θ_l_sfc, θ_i_sfc, ν_sfc, θ_r_sfc)
    return @. (1 - ice_frac) * q_over_liq + ice_frac * q_over_ice
end


"""
    ClimaLand.surface_height(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface height of the `EnergyHydrology` model.
"""
function ClimaLand.surface_height(model::EnergyHydrology{FT}, Y, p) where {FT}
    return model.domain.fields.z_sfc
end

function ClimaLand.get_drivers(model::EnergyHydrology)
    bc = model.boundary_conditions.top
    if typeof(bc) <: AtmosDrivenFluxBC{
        <:PrescribedAtmosphere,
        <:AbstractRadiativeDrivers,
        <:AbstractRunoffModel,
    }
        return (bc.atmos, bc.radiation)
    else
        return ()
    end
end

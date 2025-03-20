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
    F <: Union{FT, ClimaCore.Fields.Field},
    SF <: Union{FT, ClimaCore.Fields.Field},
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
    "Soil PAR Albedo dry"
    PAR_albedo_dry::SF
    "Soil NIR Albedo dry"
    NIR_albedo_dry::SF
    "Soil PAR Albedo wet"
    PAR_albedo_wet::SF
    "Soil NIR Albedo wet"
    NIR_albedo_wet::SF
    "Thickness of top of soil used in albedo calculations (m)"
    albedo_calc_top_thickness::FT
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
        lateral_flow::Bool = false,
    ) where {FT, D, PS}

A constructor for a `EnergyHydrology` model, which sets the default value
of the `lateral_flow` flag to false.
"""
function EnergyHydrology{FT}(;
    parameters::EnergyHydrologyParameters{FT, PSE},
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
    lateral_flow::Bool = false,
) where {FT, D, PSE}
    @assert !lateral_flow
    top_bc = boundary_conditions.top
    if typeof(top_bc) <: AtmosDrivenFluxBC
        # If the top BC indicates atmospheric conditions are driving the model
        # add baseflow as a sink term, add sublimation as a sink term
        subl_source =
            sublimation_source(Val(top_bc.prognostic_land_components), FT)
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
        @. p.soil.top_bc_wvec = Geometry.WVector(rre_top_flux_bc)
        @. p.soil.bottom_bc_wvec = Geometry.WVector(rre_bottom_flux_bc)
        divf2c_rre = Operators.DivergenceF2C(
            top = Operators.SetValue(p.soil.top_bc_wvec),
            bottom = Operators.SetValue(p.soil.bottom_bc_wvec),
        )
        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ϑ_l =
            -(divf2c_rre(-interpc2f(p.soil.K) * gradc2f(p.soil.ψ + z)))

        # Heat equation RHS
        # Reuse the same scratch space:
        @. p.soil.top_bc_wvec = Geometry.WVector(heat_top_flux_bc)
        @. p.soil.bottom_bc_wvec = Geometry.WVector(heat_bottom_flux_bc)
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(p.soil.top_bc_wvec),
            bottom = Operators.SetValue(p.soil.bottom_bc_wvec),
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
    function compute_jacobian!(
        jacobian::MatrixFields.FieldMatrixWithSolver,
        Y,
        p,
        dtγ,
        t,
    )
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
        ∂ρeres∂ϑ = matrix[@name(soil.ρe_int), @name(soil.ϑ_l)]
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
                                ν - Y.soil.θ_i, #ν_eff
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
                            ν - Y.soil.θ_i, #ν_eff
                            θ_r,
                            S_s,
                        ),
                    )
                ) - (I,)
        end
        @. ∂ρeres∂ϑ =
            -dtγ * (
                divf2c_matrix() ⋅ MatrixFields.DiagonalMatrixRow(
                    -interpc2f_op(
                        volumetric_internal_energy_liq(
                            p.soil.T,
                            model.parameters.earth_param_set,
                        ) * p.soil.K,
                    ),
                ) ⋅ gradc2f_matrix() ⋅ MatrixFields.DiagonalMatrixRow(
                    ClimaLand.Soil.dψdϑ(
                        hydrology_cm,
                        Y.soil.ϑ_l,
                        ν - Y.soil.θ_i, #ν_eff
                        θ_r,
                        S_s,
                    ),
                )
            ) - (I,)

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
            PAR_albedo_dry,
            NIR_albedo_dry,
            PAR_albedo_wet,
            NIR_albedo_wet,
        ) = model.parameters

        @. p.soil.θ_l =
            volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i, θ_r)

        update_albedo!(
            model.boundary_conditions.top,
            p,
            model.domain,
            model.parameters,
        )

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
    src::PhaseChange,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
)
    params = model.parameters
    (; ν, ρc_ds, θ_r, hydrology_cm, earth_param_set) = params
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _ρ_i = LP.ρ_cloud_ice(earth_param_set)
    Δz = model.domain.fields.Δz # center face distance
    _LH_f0 = LP.LH_f0(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    _grav = LP.grav(earth_param_set)

    # Pass in physical constants as floats, to avoid issue broadcasting over fields and structs
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
                Δz,
                p.soil.κ,
            ),
            ν,
            θ_r,
            hydrology_cm,
            _ρ_i,
            _ρ_l,
            _LH_f0,
            _T_freeze,
            _grav,
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
                Δz,
                p.soil.κ,
            ),
            ν,
            θ_r,
            hydrology_cm,
            _ρ_i,
            _ρ_l,
            _LH_f0,
            _T_freeze,
            _grav,
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
        -p.soil.turbulent_fluxes.vapor_flux_ice * _ρ_l / _ρ_i *
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
    return @. (p.soil.PAR_albedo + p.soil.NIR_albedo) / 2
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



function turbulent_fluxes!(
    dest,
    atmos::PrescribedAtmosphere,
    model::EnergyHydrology{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
) where {FT}
    # Obtain surface quantities needed for computation; these should not allocate
    T_sfc = ClimaLand.surface_temperature(model, Y, p, t)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    d_sfc = ClimaLand.displacement_height(model, Y, p)
    u_air = p.drivers.u
    h_air = atmos.h
    (; K_sat, ν, θ_r, hydrology_cm, z_0m, z_0b, Ω, γ, γT_ref, earth_param_set) =
        model.parameters
    hydrology_cm_sfc = ClimaLand.Domains.top_center_to_surface(hydrology_cm)
    K_sat_sfc = ClimaLand.Domains.top_center_to_surface(K_sat)
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    ν_sfc = ClimaLand.Domains.top_center_to_surface(ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(θ_r)
    θ_l_sfc = p.soil.sfc_scratch
    ClimaLand.Domains.linear_interpolation_to_surface!(
        θ_l_sfc,
        p.soil.θ_l,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    dest .=
        soil_turbulent_fluxes_at_a_point.(
            T_sfc,
            θ_l_sfc,
            θ_i_sfc,
            h_sfc,
            d_sfc,
            hydrology_cm_sfc,
            ν_sfc,
            θ_r_sfc,
            K_sat_sfc,
            p.drivers.thermal_state,
            u_air,
            h_air,
            atmos.gustiness,
            z_0m,
            z_0b,
            Ω,
            γ,
            γT_ref,
            Ref(earth_param_set),
        )
    return nothing
end

"""
    soil_turbulent_fluxes_at_a_point(
                                T_sfc::FT,
                                θ_l_sfc::FT,
                                θ_i_sfc::FT,
                                h_sfc::FT,
                                d_sfc::FT,
                                hydrology_cm_sfc::C,
                                ν_sfc::FT,
                                θ_r_sfc::FT,
                                K_sat_sfc::FT,
                                thermal_state_air::Thermodynamics.PhaseEquil{FT},
                                u_air::FT,
                                h_air::FT,
                                gustiness::FT,
                                z_0m::FT,
                                z_0b::FT,
                                Ω::FT,
                                γ::FT,
                                γT_ref,::FT
                                earth_param_set::EP
                               ) where {FT <: AbstractFloat, C, EP}

Computes turbulent surface fluxes for soil at a point on a surface given
(1) Surface state conditions (`T_sfc`, `θ_l_sfc`, `θ_i_sfc`)
(2) Surface properties, such as the topographical height of the surface (`h_sfc`),
    the displacement height (`d_sfc`), hydraulic parameters (`hydrology_cm_sfc`,
    `ν_sfc, `θ_r_sfc`, `K_sat_sfc`)
(4) Atmospheric state conditions (`thermal_state_air`, `u_air`)
(5) Height corresponding to where atmospheric state is measured (`h_air`)
(6) Parameters: `gustiness`, roughness lengths `z_0m`, `z_0b`, several
    required to compute the soil conductivity `Ω`, γ`, γT_ref`, and the
    `earth_param_set`

This returns an energy flux and a  water volume flux, stored in
a tuple with self explanatory keys. If the temperature is above the freezing point,
the vapor flux comes from liquid water; if the temperature is below the freezing
point, it comes from the soil ice.
"""
function soil_turbulent_fluxes_at_a_point(
    T_sfc::FT,
    θ_l_sfc::FT,
    θ_i_sfc::FT,
    h_sfc::FT,
    d_sfc::FT,
    hydrology_cm_sfc,
    ν_sfc::FT,
    θ_r_sfc::FT,
    K_sat_sfc::FT,
    thermal_state_air::Thermodynamics.PhaseEquil{FT},
    u_air::FT,
    h_air::FT,
    gustiness::FT,
    z_0m::FT,
    z_0b::FT,
    Ω::FT,
    γ::FT,
    γT_ref::FT,
    earth_param_set::P,
) where {FT <: AbstractFloat, P}
    # Parameters
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    _LH_v0::FT = LP.LH_v0(earth_param_set)
    _LH_f0::FT = LP.LH_f0(earth_param_set)
    _ρ_liq::FT = LP.ρ_cloud_liq(earth_param_set)
    _ρ_ice::FT = LP.ρ_cloud_ice(earth_param_set)
    _T_freeze::FT = LP.T_freeze(earth_param_set)
    _grav::FT = LP.grav(earth_param_set)

    # Atmos air state
    state_air = SurfaceFluxes.StateValues(
        h_air - d_sfc - h_sfc,
        SVector{2, FT}(u_air, 0),
        thermal_state_air,
    )
    q_air::FT =
        Thermodynamics.total_specific_humidity(thermo_params, thermal_state_air)

    # Estimate surface air density from the atmos state using T_sfc
    ρ_sfc::FT = ClimaLand.compute_ρ_sfc(thermo_params, thermal_state_air, T_sfc)

    # Get the freezing point of the soil to determine if the water is sublimating or evaporating
    θtot_sfc = min(_ρ_ice / _ρ_liq * θ_i_sfc + θ_l_sfc, ν_sfc)
    # This is consistent with Equation (22) of Dall'Amico
    ψw0_sfc = matric_potential(
        hydrology_cm_sfc,
        effective_saturation(ν_sfc, θtot_sfc, θ_r_sfc),
    )
    Tf_depressed = _T_freeze * exp(_grav * ψw0_sfc / _LH_f0)

    # The following will be reset below
    β::FT = 1
    Ẽ_i::FT = 0 # vapor flux of liquid water, in units of volume flux of liquid water
    Ẽ_l::FT = 0 # vapor flux of frozen water, in units of volume flux of liquid water

    # Compute q_soil using ice or liquid as appropriate, and create the thermal state of the soil
    # For liquid water evap, β = 1, and for ice, β is a numerical factor which damps sublimation to zero as ice goes to zero,
    if T_sfc > Tf_depressed # liquid water evaporation
        liquid_evaporation = true
        q_sat_liq::FT = Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Liquid(),
        )
        thermal_state_sfc::Thermodynamics.PhaseEquil{FT} =
            Thermodynamics.PhaseEquil_ρTq(
                thermo_params,
                ρ_sfc,
                T_sfc,
                q_sat_liq,
            )
    else
        liquid_evaporation = false
        q_sat_ice::FT = Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Ice(),
        )
        thermal_state_sfc = Thermodynamics.PhaseEquil_ρTq(
            thermo_params,
            ρ_sfc,
            T_sfc,
            q_sat_ice,
        )
        if q_air < q_sat_ice
            β *= (θ_i_sfc / ν_sfc)^4
        end
    end

    # Now we compute the fluxes (E, LH, SH)
    # Thermal state and β encode ice vs liquid water vapor flux differences

    # SurfaceFluxes.jl expects a relative difference between where u_air = 0
    # and the atmosphere height. Here, we assume h and h_sfc are measured
    # relative to a common reference. Then d_sfc + h_sfc + z_0m is the apparent
    # source of momentum, and
    # Δh ≈ h_air - d_sfc - h_sfc is the relative height difference between the
    # apparent source of momentum and the atmosphere height.

    # In this we have neglected z_0m and z_0b (i.e. assumed they are small
    # compared to Δh).
    state_sfc = SurfaceFluxes.StateValues(
        FT(0),
        SVector{2, FT}(0, 0),
        thermal_state_sfc,
    )
    # State containers
    states = SurfaceFluxes.ValuesOnly(
        state_air,
        state_sfc,
        z_0m,
        z_0b,
        gustiness = gustiness,
        beta = β,
    )
    scheme = SurfaceFluxes.PointValueScheme()
    conditions =
        SurfaceFluxes.surface_conditions(surface_flux_params, states, scheme)

    SH::FT = SurfaceFluxes.sensible_heat_flux(
        surface_flux_params,
        conditions.Ch,
        states,
        scheme,
    )
    r_ae::FT = 1 / (conditions.Ch * SurfaceFluxes.windspeed(states))
    E_pot::FT =
        SurfaceFluxes.evaporation(surface_flux_params, states, conditions.Ch)# potential evaporation rate, mass flux
    Ẽ_pot::FT = E_pot / _ρ_liq # volume flux of liquid water

    # Adjust fluxes as needed for soil resistance; split between ice and liquid water loss
    if liquid_evaporation     # adjust the vapor loss to account for soil resistance, set sublimation to zero
        Ẽ_i = 0
        Ẽ_l = Ẽ_pot
        if q_air < q_sat_liq # adjust potential evaporation rate to account for soil resistance
            K_sfc::FT =
                impedance_factor(θ_i_sfc / (θ_l_sfc + θ_i_sfc - θ_r_sfc), Ω) *
                viscosity_factor(T_sfc, γ, γT_ref) *
                hydraulic_conductivity(
                    hydrology_cm_sfc,
                    K_sat_sfc,
                    effective_saturation(ν_sfc, θ_l_sfc, θ_r_sfc),
                )
            K_c::FT = hydraulic_conductivity(
                hydrology_cm_sfc,
                K_sat_sfc,
                hydrology_cm_sfc.S_c,
            )
            x::FT = 4 * K_sfc * (1 + Ẽ_pot / (4 * K_c))
            Ẽ_l *= x / (Ẽ_pot + x)
        end
    else # sublimation, set evaporation to zero
        Ẽ_i = Ẽ_pot
        Ẽ_l = 0
    end

    # Heat fluxes for soil
    LH::FT = _LH_v0 * (Ẽ_l + Ẽ_i) * _ρ_liq
    # Derivatives
    dshfdT::FT = 0
    dlhfdT::FT = 0
    return (
        lhf = LH,
        shf = SH,
        vapor_flux_liq = Ẽ_l,
        r_ae = r_ae,
        vapor_flux_ice = Ẽ_i,
    )
end

"""
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model::EnergyHydrology,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value for
the total liquid water volume per unit ground area for the `EnergyHydrology`.

The water in all phases is accounted for by converting ice volume to liquid water
volume using the ratio of the density of ice to the density of water.
"""
function ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::EnergyHydrology,
    Y,
    p,
    t,
)
    earth_param_set = model.parameters.earth_param_set
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    _ρ_ice = LP.ρ_cloud_ice(earth_param_set)
    ClimaCore.Operators.column_integral_definite!(
        surface_field,
        Y.soil.ϑ_l .+ Y.soil.θ_i .* _ρ_ice ./ _ρ_liq, # this line allocates
    )
    return nothing
end

"""
    ClimaLand.total_energy_per_area!(
        surface_field,
        model::EnergyHydrology,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value for
the total energy per unit ground area for the `EnergyHydrology`.
"""
function ClimaLand.total_energy_per_area!(
    surface_field,
    model::EnergyHydrology,
    Y,
    p,
    t,
)
    ClimaCore.Operators.column_integral_definite!(surface_field, Y.soil.ρe_int)
    return nothing
end

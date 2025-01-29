"""
    RichardsParameters{F <: Union{<: AbstractFloat, ClimaCore.Fields.Field}, C <: AbstractSoilHydrologyClosure}

A struct for storing parameters of the `RichardsModel`.
$(DocStringExtensions.FIELDS)
"""
struct RichardsParameters{
    F <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
    C,
}
    "The porosity of the soil (m^3/m^3)"
    ν::F
    "The hydrology closure model: vanGenuchten or BrooksCorey"
    hydrology_cm::C
    "The saturated hydraulic conductivity (m/s)"
    K_sat::F
    "The specific storativity (1/m)"
    S_s::F
    "The residual water fraction (m^3/m^3"
    θ_r::F
end

function RichardsParameters(;
    hydrology_cm::C,
    ν::F,
    K_sat::F,
    S_s::F,
    θ_r::F,
) where {F <: Union{<:AbstractFloat, ClimaCore.Fields.Field}, C}
    return RichardsParameters{F, typeof(hydrology_cm)}(
        ν,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
    )
end

"""
    RichardsModel

A model for simulating the flow of water in a porous medium
by solving the Richardson-Richards Equation.

A variety of boundary condition types are supported, including
FluxBC, RichardsAtmosDrivenFluxBC, MoistureStateBC, and FreeDrainage
(only for the bottom of the domain).

If you wish to
simulate soil hydrology under the context of a prescribed precipitation
volume flux (m/s) as a function of time, the RichardsAtmosDrivenFluxBC
type should be chosen. Please see the documentation for more details.

$(DocStringExtensions.FIELDS)
"""
struct RichardsModel{FT, PS, D, BC, S} <: AbstractSoilModel{FT}
    "the parameter set"
    parameters::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BC
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
    "A boolean flag which, when false, turns off the horizontal flow of water"
    lateral_flow::Bool
end

"""
    RichardsModel{FT}(;
        parameters::RichardsParameters,
        domain::D,
        boundary_conditions::NamedTuple,
        sources::Tuple,
        lateral_flow::Bool = true
    ) where {FT, D}

A constructor for a `RichardsModel`, which sets the
default value of `lateral_flow` to be true.
"""
function RichardsModel{FT}(;
    parameters::RichardsParameters,
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
    lateral_flow::Bool = false,
) where {FT, D}
    top_bc = boundary_conditions.top
    if typeof(top_bc) <: RichardsAtmosDrivenFluxBC
        # If the top BC indicates precipitation is driving the model,
        # add baseflow as a sink/source term
        subsurface_source = subsurface_runoff_source(top_bc.runoff)
        sources = append_source(subsurface_source, sources)
    end
    @assert !lateral_flow
    args = (parameters, domain, boundary_conditions, sources)
    RichardsModel{FT, typeof.(args)...}(args..., lateral_flow)
end

function make_update_boundary_fluxes(model::RichardsModel)
    function update_boundary_fluxes!(p, Y, t)
        z = model.domain.fields.z
        Δz_top = model.domain.fields.Δz_top
        Δz_bottom = model.domain.fields.Δz_bottom
        boundary_flux!(
            p.soil.top_bc,
            model.boundary_conditions.top,
            TopBoundary(),
            model,
            Δz_top,
            Y,
            p,
            t,
        )
        boundary_flux!(
            p.soil.bottom_bc,
            model.boundary_conditions.bottom,
            BottomBoundary(),
            model,
            Δz_bottom,
            Y,
            p,
            t,
        )

        # Update `p.soil.dfluxBCdY`
        set_dfluxBCdY!(
            model,
            model.boundary_conditions.top,
            TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
    end
    return update_boundary_fluxes!
end

"""
    make_compute_imp_tendency(model::RichardsModel)

An extension of the function `make_compute_imp_tendency`, for the Richardson-
Richards equation.

This function creates and returns a function which computes the entire
right hand side of the PDE for `ϑ_l`, and updates `dY.soil.ϑ_l` in place
with that value.
"""
function ClimaLand.make_compute_imp_tendency(model::RichardsModel)
    function compute_imp_tendency!(dY, Y, p, t)
        z = model.domain.fields.z
        top_flux_bc = p.soil.top_bc
        bottom_flux_bc = p.soil.bottom_bc
        @. p.soil.top_bc_wvec = Geometry.WVector(top_flux_bc)
        @. p.soil.bottom_bc_wvec = Geometry.WVector(bottom_flux_bc)
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()

        # We are setting a boundary value on a flux, which is a gradient of a scalar
        # Therefore, we should set boundary conditions in terms of a covariant vector
        # We set the third component first - supply a Covariant3Vector

        # Without topography only
        # In Cartesian coordinates, W (z^) = Cov3 (z^) = Contra3 (n^ = z^)
        # In spherical coordinates, W (r^) = Cov3 (r^) = Contra3 (n^ = r^)

        # It appears that the WVector is converted internally to a Covariant3Vector for the gradient value
        # at the boundary. Offline tests indicate that you get the same thing if
        # the bc is WVector(F) or Covariant3Vector(F*Δr) or Contravariant3Vector(F/Δr)

        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(p.soil.top_bc_wvec),
            bottom = Operators.SetValue(p.soil.bottom_bc_wvec),
        )

        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ϑ_l =
            -(divf2c_water(-interpc2f(p.soil.K) * gradc2f_water(p.soil.ψ + z)))
    end
    return compute_imp_tendency!
end

"""
    make_explicit_tendency(model::Soil.RichardsModel)

An extension of the function `make_compute_imp_tendency`, for the Richardson-
Richards equation.

Construct the tendency computation function for the explicit terms of the RHS,
which are horizontal components and source/sink terms.
"""
function ClimaLand.make_compute_exp_tendency(model::Soil.RichardsModel)
    # Currently, boundary conditions in the horizontal conditions
    # are restricted to be periodic. In this case, the explicit tendency
    # does not depend on boundary fluxes, and we do not need to update
    # the boundary_var variables prior to evaluation.
    function compute_exp_tendency!(dY, Y, p, t)
        # set dY before updating it
        dY.soil.ϑ_l .= eltype(dY.soil.ϑ_l)(0)
        z = model.domain.fields.z

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
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Union{HybridBox, SphericalShell},
                          lateral_flow::Val{true},
                          model::RichardsModel,
                          p::NamedTuple)

Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators for the RichardsModel,
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
    model::RichardsModel,
    p::NamedTuple,
    z::ClimaCore.Fields.Field,
)
    hdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    # The flux is already covariant, from hgrad, so no need to convert.
    @. dY.soil.ϑ_l += -hdiv(-p.soil.K * hgrad(p.soil.ψ + z))
end

"""
    prognostic_vars(soil::RichardsModel)

A function which returns the names of the prognostic variables
of `RichardsModel`.
"""
ClimaLand.prognostic_vars(soil::RichardsModel) = (:ϑ_l,)
ClimaLand.prognostic_types(soil::RichardsModel{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(soil::RichardsModel) = (:subsurface,)

"""
    auxiliary_vars(soil::RichardsModel)

A function which returns the names of the auxiliary variables
of `RichardsModel`.
"""
function ClimaLand.auxiliary_vars(soil::RichardsModel)
    return (
        :K,
        :ψ,
        boundary_vars(soil.boundary_conditions.top, ClimaLand.TopBoundary())...,
        boundary_vars(
            soil.boundary_conditions.bottom,
            ClimaLand.BottomBoundary(),
        )...,
    )
end

"""
    auxiliary_domain_names(soil::RichardsModel)

A function which returns the names of the auxiliary domain names
of `RichardsModel`.
"""
function ClimaLand.auxiliary_domain_names(soil::RichardsModel)
    return (
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
end

"""
    auxiliary_types(soil::RichardsModel)

A function which returns the names of the auxiliary types
of `RichardsModel`.
"""
function ClimaLand.auxiliary_types(soil::RichardsModel{FT}) where {FT}
    return (
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
end

"""
    make_update_aux(model::RichardsModel)

An extension of the function `make_update_aux`, for the Richardson-
Richards equation.

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_update_aux(model::RichardsModel)
    function update_aux!(p, Y, t)
        (; ν, hydrology_cm, K_sat, S_s, θ_r) = model.parameters
        @. p.soil.K = hydraulic_conductivity(
            hydrology_cm,
            K_sat,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(hydrology_cm, θ_r, Y.soil.ϑ_l, ν, S_s)
    end
    return update_aux!
end

"""
    ClimaLand.make_compute_jacobian(model::RichardsModel{FT}) where {FT}

Creates and returns the compute_jacobian! function for RichardsModel.
This updates the contribution for the soil liquid water content.

Using this Jacobian with a backwards Euler timestepper is equivalent
to using the modified Picard scheme of Celia et al. (1990).
"""
function ClimaLand.make_compute_jacobian(model::RichardsModel{FT}) where {FT}
    function compute_jacobian!(
        jacobian::MatrixFields.FieldMatrixWithSolver,
        Y,
        p,
        dtγ,
        t,
    )
        (; matrix) = jacobian
        (; ν, hydrology_cm, S_s, θ_r) = model.parameters

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
    end
    return compute_jacobian!
end

"""
    ClimaLand.get_drivers(model::RichardsModel)

Returns the driver variable symbols for the RichardsModel; these
depend on the boundary condition type and currently only are required
for the RichardsAtmosDrivenFluxBC, which is driven by
a prescribed time and space varying precipitation.
"""
function ClimaLand.get_drivers(model::RichardsModel)
    bc = model.boundary_conditions.top
    if typeof(bc) <: RichardsAtmosDrivenFluxBC{
        <:PrescribedPrecipitation,
        <:AbstractRunoffModel,
    }
        return (bc.precip,)
    else
        return ()
    end
end

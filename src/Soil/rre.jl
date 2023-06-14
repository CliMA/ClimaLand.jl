"""
    RichardsParameters{FT <: AbstractFloat}

A struct for storing parameters of the `RichardModel`.
$(DocStringExtensions.FIELDS)
"""
struct RichardsParameters{
    FT <: AbstractFloat,
    C <: AbstractSoilHydrologyClosure,
}
    "The porosity of the soil (m^3/m^3)"
    ν::FT
    "The hydrology closure model: vanGenuchten or BrooksCorey"
    hydrology_cm::C
    "The saturated hydraulic conductivity (m/s)"
    K_sat::FT
    "The specific storativity (1/m)"
    S_s::FT
    "The residual water fraction (m^3/m^3"
    θ_r::FT
end

function RichardsParameters(;
    hydrology_cm::C,
    ν::FT,
    K_sat::FT,
    S_s::FT,
    θ_r::FT,
) where {FT <: AbstractFloat, C <: AbstractSoilHydrologyClosure}
    return RichardsParameters{FT, typeof(hydrology_cm)}(
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
end

"""
    RichardsModel{FT}(;
        parameters::RichardsParameters{FT},
        domain::D,
        boundary_conditions::NamedTuple,
        sources::Tuple,
    ) where {FT, D}

A constructor for a `RichardsModel`.
"""
function RichardsModel{FT}(;
    parameters::RichardsParameters{FT},
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
) where {FT, D}
    args = (parameters, domain, boundary_conditions, sources)
    RichardsModel{FT, typeof.(args)...}(args...)
end


"""
    make_compute_exp_tendency(model::RichardsModel)

An extension of the function `make_compute_exp_tendency`, for the Richardson-
Richards equation.

This function creates and returns a function which computes the entire
right hand side of the PDE for `ϑ_l`, and updates `dY.soil.ϑ_l` in place
with that value.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_compute_exp_tendency(model::RichardsModel)
    function compute_exp_tendency!(dY, Y, p, t)
        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        Δz_top, Δz_bottom = get_Δz(z)

        top_flux_bc = boundary_flux(
            model.boundary_conditions.top.water,
            TopBoundary(),
            Δz_top,
            p,
            t,
            model.parameters,
        )
        bot_flux_bc = boundary_flux(
            model.boundary_conditions.bottom.water,
            BottomBoundary(),
            Δz_bottom,
            p,
            t,
            model.parameters,
        )

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
            top = Operators.SetValue(Geometry.WVector.(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(bot_flux_bc)),
        )

        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ϑ_l =
            -(divf2c_water(-interpc2f(p.soil.K) * gradc2f_water(p.soil.ψ + z)))
        # Horizontal contributions
        horizontal_components!(dY, model.domain, model, p, z)

        # Source terms
        for src in model.sources
            ClimaLSM.source!(dY, src, Y, p, model.parameters)
        end
    end
    return compute_exp_tendency!
end


"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::HybridBox,
                          model::RichardsModel,
                          p::ClimaCore.Fields.FieldVector)
Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a hybrid box domain, the horizontal contributions are
computed using the WeakDivergence and Gradient operators.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
    model::RichardsModel,
    p::ClimaCore.Fields.FieldVector,
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
ClimaLSM.prognostic_vars(soil::RichardsModel) = (:ϑ_l,)
ClimaLSM.prognostic_types(soil::RichardsModel{FT}) where {FT} = (FT,)
"""
    auxiliary_vars(soil::RichardsModel)

A function which returns the names of the auxiliary variables
of `RichardsModel`.

Note that auxiliary variables are not needed for such a simple model.
We could instead compute the conductivity and matric potential within the
tendency function explicitly, rather than compute and store them in the
auxiliary vector `p`. We did so in this case as a demonstration.
"""
ClimaLSM.auxiliary_vars(soil::RichardsModel) = (:K, :ψ)
ClimaLSM.auxiliary_types(soil::RichardsModel{FT}) where {FT} = (FT, FT)
"""
    make_update_aux(model::RichardsModel)

An extension of the function `make_update_aux`, for the Richardson-
Richards equation.

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::RichardsModel)
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
    RichardsTridiagonalW{R, J, W, T} <: ClimaLSM.AbstractTridiagonalW

A struct containing the necessary information for constructing a tridiagonal
Jacobian matrix (`W`) for solving Richards equation, treating only the vertical
diffusion term implicitly.

Note that the diagonal, upper diagonal, and lower diagonal entry values
are stored in this struct and updated in place.
$(DocStringExtensions.FIELDS)
"""
struct RichardsTridiagonalW{R, J, JA, T, A} <: ClimaLSM.AbstractTridiagonalW
    "Reference to dtγ, which is specified by the ODE solver"
    dtγ_ref::R
    "Diagonal entries of the Jacobian stored as a ClimaCore.Fields.Field"
    ∂ϑₜ∂ϑ::J
    "Array of tridiagonal matrices containing W for each column"
    W_column_arrays::JA
    "An allocated cache used to evaluate ldiv!"
    temp1::T
    "An allocated cache used to evaluate ldiv!"
    temp2::T
    "A flag indicating whether this struct is used to compute Wfact_t or Wfact"
    transform::Bool
    "A pre-allocated cache storing ones on the face space"
    ones_face_space::A
end

"""
    RichardsTridiagonalW(
        Y::ClimaCore.Fields.FieldVector;
        transform::Bool = false
)

Outer constructor for the RichardsTridiagonalW Jacobian
matrix struct.

Initializes all variables to zeros.
"""
function RichardsTridiagonalW(
    Y::ClimaCore.Fields.FieldVector;
    transform::Bool = false,
)
    FT = eltype(Y.soil.ϑ_l)
    center_space = axes(Y.soil.ϑ_l)
    N = Spaces.nlevels(center_space)

    tridiag_type = Operators.StencilCoefs{-1, 1, NTuple{3, FT}}
    ∂ϑₜ∂ϑ = Fields.Field(tridiag_type, center_space)

    ∂ϑₜ∂ϑ.coefs[1] .= FT(0)
    ∂ϑₜ∂ϑ.coefs[2] .= FT(0)
    ∂ϑₜ∂ϑ.coefs[3] .= FT(0)

    W_column_arrays = [
        LinearAlgebra.Tridiagonal(
            zeros(FT, N - 1) .+ FT(0),
            zeros(FT, N) .+ FT(0),
            zeros(FT, N - 1) .+ FT(0),
        ) for _ in 1:Threads.nthreads()
    ]
    dtγ_ref = Ref(FT(0))
    temp1 = similar(Y.soil.ϑ_l)
    temp1 .= FT(0)
    temp2 = similar(Y.soil.ϑ_l)
    temp2 .= FT(0)

    face_space = ClimaLSM.Domains.obtain_face_space(center_space)
    ones_face_space = ones(face_space)

    return RichardsTridiagonalW(
        dtγ_ref,
        ∂ϑₜ∂ϑ,
        W_column_arrays,
        temp1,
        temp2,
        transform,
        ones_face_space,
    )
end


"""
    ClimaLSM.make_update_jacobian(model::RichardsModel)

Creates and returns the update_jacobian! function for RichardsModel.

Using this Jacobian with a backwards Euler timestepper is equivalent
to using the modified Picard scheme of Celia et al. (1990).
"""
function ClimaLSM.make_update_jacobian(model::RichardsModel)
    function update_jacobian!(W::RichardsTridiagonalW, Y, p, dtγ, t)
        FT = eltype(Y.soil.ϑ_l)
        (; dtγ_ref, ∂ϑₜ∂ϑ) = W
        (; ν, hydrology_cm, S_s, θ_r) = model.parameters
        divf2c_op = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(FT(0))),
            bottom = Operators.SetValue(Geometry.WVector.(FT(0))),
        )
        divf2c_stencil = Operators.Operator2Stencil(divf2c_op)
        gradc2f_op = Operators.GradientC2F(
            top = Operators.SetGradient(Geometry.WVector.(FT(0))),
            bottom = Operators.SetGradient(Geometry.WVector.(FT(0))),
        )
        gradc2f_stencil = Operators.Operator2Stencil(gradc2f_op)
        interpc2f_op = Operators.InterpolateC2F(
            bottom = Operators.Extrapolate(),
            top = Operators.Extrapolate(),
        )
        compose = Operators.ComposeStencils()

        @. ∂ϑₜ∂ϑ = compose(
            divf2c_stencil(Geometry.Covariant3Vector(W.ones_face_space)),
            (
                interpc2f_op(p.soil.K) * ClimaLSM.to_scalar_coefs(
                    gradc2f_stencil(
                        ClimaLSM.Soil.dψdϑ(
                            hydrology_cm,
                            Y.soil.ϑ_l,
                            ν,
                            θ_r,
                            S_s,
                        ),
                    ),
                )
            ),
        )
        # Hardcoded for single column: FIX!
        z = Fields.coordinate_field(axes(Y.soil.ϑ_l)).z
        Δz_top, Δz_bottom = get_Δz(z)
        ∂T_bc∂YN = ClimaLSM.∂tendencyBC∂Y(
            model,
            model.boundary_conditions.top.water,
            ClimaLSM.TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
        #TODO: allocate space in W? See how final implementation of stencils with boundaries works out
        N = ClimaCore.Spaces.nlevels(axes(Y.soil.ϑ_l))
        parent(ClimaCore.Fields.level(∂ϑₜ∂ϑ.coefs.:2, N)) .=
            parent(ClimaCore.Fields.level(∂ϑₜ∂ϑ.coefs.:2, N)) .+
            parent(∂T_bc∂YN.soil.ϑ_l)

    end
    return update_jacobian!
end

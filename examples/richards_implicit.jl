using LinearAlgebra

using ClimaCore: Operators, Geometry, Spaces, Fields
using ClimaCore.Utilities: half
using ClimaLSM
using ClimaLSM.Soil: RichardsModel, RichardsParameters, FluxBC
using ClimaLSM.Domains: HybridBox, Column
using ClimaCore.Utilities: half
using UnPack
const compose = Operators.ComposeStencils()
const apply = Operators.ApplyStencil()
# Allow one() to be called on vectors.
Base.one(::T) where {T <: Geometry.AxisTensor} = one(T)
Base.one(::Type{T}) where {T′, A, S, T <: Geometry.AxisTensor{T′, 1, A, S}} =
    T(axes(T), S(one(T′)))

# As long as the boundary conditions are independent of state, we
# don't care about their values in computing Operator2Stencil
FT = Float64
gradc2f_op = Operators.GradientC2F(
    bottom = Operators.SetValue(FT(0.0)),
    top = Operators.SetValue(FT(0.0))
)
gradc2f_stencil = Operators.Operator2Stencil(gradc2f_op)
divf2c_op = Operators.DivergenceF2C(
    top = Operators.SetValue(Geometry.Covariant3Vector(FT(0))),
    bottom = Operators.SetValue(Geometry.Covariant3Vector(FT(0))),
)
divf2c_stencil = Operators.Operator2Stencil(divf2c_op)
interpc2f_op = Operators.InterpolateC2F(
    bottom = Operators.Extrapolate(),
    top = Operators.Extrapolate(),
)
interpc2f_stencil = Operators.Operator2Stencil(interpc2f_op)

struct TridiagonalJacobian{R, J1}
    # reference to dtγ, which is specified by the ODE solver
    dtγ_ref::R
    ∂ϑt∂ϑ::J1
    flag::Symbol
end

function TridiagonalJacobian(Y, flag::Symbol)
    FT = eltype(Y.soil.ϑ_l)
    dtγ_ref = Ref(zero(FT))
    space = axes(Y.soil.ϑ_l)
    Jtype = Operators.StencilCoefs{-1, 1, NTuple{3, FT}}
    J = Fields.Field(Jtype, space)
    args = (dtγ_ref,J)
    return TridiagonalJacobian{typeof.(args)...}(args..., flag)
end
Base.similar(w::TridiagonalJacobian) = w
to_scalar_coefs(vector_coefs) =
    map(vector_coef-> vector_coef.u₃, vector_coefs)
function make_Wfact(model::RichardsModel)
    function Wfact!(W, Y, p, dtγ, t)
        (; dtγ_ref, ∂ϑt∂ϑ, flag) = W
        dtγ_ref[] = dtγ
        # ∂ϑ∂t = ϑt = div(I(K(θ)) grad(h(θ)))
        if flag == :exact
            @. ∂ϑt∂ϑ = compose(
                divf2c_stencil(one.(gradc2f_op.(p.soil.ψ + model.coordinates.z))),
                to_scalar_coefs(
                    interpc2f_op(dKdθ((Y.soil.ϑ_l - θ_r)/(ν-θ_r),
                                      ν,
                                      θ_r,
                                      vg_m,
                                      K_sat
                                      )
                                 )*interpc2f_stencil(one(Y.soil.ϑ_l))*
                    gradc2f_op(p.soil.ψ + model.coordinates.z) +
                    interpc2f_op(p.soil.K*dψdθ((Y.soil.ϑ_l - θ_r)/(ν-θ_r),
                                               ν,
                                               θ_r,
                                               vg_α,
                                               vg_n,
                                               vg_m,
                                               S_s
                                               )
                                 )* gradc2f_stencil(one.(Y.soil.ϑ_l))
                )
            )
        elseif flag == :modifiedpicard
            @. ∂ϑt∂ϑ = compose(
                divf2c_stencil(one.(gradc2f_op.(p.soil.ψ + model.coordinates.z))),
                to_scalar_coefs(
                    interpc2f_op(p.soil.K*dψdθ((Y.soil.ϑ_l-θ_r)/(ν-θ_r),
                                               ν,
                                               θ_r,
                                               vg_α,
                                               vg_n,
                                               vg_m,
                                               S_s)
                                 )*gradc2f_stencil(one.(Y.soil.ϑ_l))
                )
            )
        else
            println("Unsupported Jacobian approximation")
        end
        

        # CLM approach
        # Newton's method with one iteration, and use the previous time step as your first guess.
        
    end
    return Wfact!
end
#S = (θ - θ_r) / (ν - θ_r), but the deriv is with respect to θ
function dKdθ(S::FT, ν::FT, θ_r::FT, vg_m::FT, Ksat::FT)::FT where {FT}
    if S < 1
        return (((1.0 - (1.0 - S^(1.0 / vg_m))^vg_m))^2.0 / 2.0 / S^0.5
                + 2 * S^(1 / vg_m - 1 / 2) * (1.0 - (1.0 - S^(1.0 / vg_m))^vg_m) / (1 - S^(1 / vg_m))^(1.0 - vg_m)
                )* K_sat / (ν - θ_r)
    else
        return 0.0
    end
    
end
#S = (θ - θ_r) / (ν - θ_r), but the deriv is with respect to θ
function dψdθ(S::FT,  ν::FT, θ_r::FT, vg_α::FT, vg_n::FT, vg_m::FT, S_s::FT)::FT where {FT}
    if S < 1.0
        return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
            (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
            S^(-1 / vg_m - 1)
    else
        return 1.0 / S_s
    end
end


function make_implicit_tendency(model::RichardsModel)
    function implicit_tendency!(dY,Y,p,t)
        @unpack ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r = model.parameters
        # when does aux get updated?
        @. p.soil.K = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
        top_flux_bc, bot_flux_bc =
            boundary_fluxes(model.boundary_conditions, p, t)
        z = model.coordinates.z
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(bot_flux_bc)),
        )
        @. dY.soil.ϑ_l =
            -(divf2c_water(-interpc2f(p.soil.K) * gradc2f_water(p.soil.ψ + z)))
    end
    return implicit_tendency!
end


function make_explicit_tendency(model::RichardsModel)
    function explicit_tendency!(dY,Y,p,t)
        #when does aux get updated?
        @. p.soil.K = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
        hdiv = Operators.WeakDivergence()
        hgrad = Operators.Gradient()
        @. dY.soil.ϑ_l += -hdiv(-p.soil.K * hgrad(p.soil.ψ + model.coordinates.z))
        Spaces.weighted_dss!(dY.soil.ϑ_l)
    end
    return explicit_tendency!
end
#domain = HybridBox(;
#                   zlim = (0.0,1.0),
#                   xlim = (0.0,100.0),
#                   ylim = (0.0,100.0),
#                   nelements = (10,10,10),
#                   npolynomial = 1,
#                   periodic = (true,true)
#                   )
ν = FT(0.495)
K_sat = FT(0.0443 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(2.0)
vg_α = FT(2.6) # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0)
zmax = FT(0)
zmin = FT(-10)
nelems = 10

soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
top_flux_bc = FT(0.0)
bot_flux_bc = FT(0.0)
sources = ()
boundary_conditions = FluxBC{FT}(top_flux_bc, bot_flux_bc)
params = RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

soil = RichardsModel{FT}(;
                              parameters = params,
                              domain = soil_domain,
                              boundary_conditions = boundary_conditions,
                              sources = sources,
                              )

Y, p, coords = initialize(soil)
function init_soil!(Y, z, params)
    function hydrostatic_profile(
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack ν, vg_α, vg_n, vg_m, θ_r = params
        #unsaturated zone only, assumes water table starts at z_∇
        z_∇ = FT(-10)# matches zmin
        S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
        ϑ_l = S * (ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Y.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
end

init_soil!(Y, coords.z, soil.parameters)
update_aux! = make_update_aux(soil)
update_aux!(p, Y, 0.0)

W = TridiagonalJacobian(Y, :exact)
Wfact! = make_Wfact(soil)


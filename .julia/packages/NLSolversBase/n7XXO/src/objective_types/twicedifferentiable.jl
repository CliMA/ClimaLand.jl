# Used for objectives and solvers where the gradient and Hessian is available/exists
mutable struct TwiceDifferentiable{T,TDF,TH,TX} <: AbstractObjective
    f
    df
    fdf
    dfh
    fdfh
    h
    F::T
    DF::TDF
    H::TH
    x_f::TX
    x_df::TX
    x_h::TX
    f_calls::Vector{Int}
    df_calls::Vector{Int}
    h_calls::Vector{Int}
end
# compatibility with old constructor
function TwiceDifferentiable(f, g, fg, h, x::TX, F::T = real(zero(eltype(x))), G::TG = alloc_DF(x, F), H::TH = alloc_H(x, F); inplace = true) where {T, TG, TH, TX}
    x_f, x_df, x_h = x_of_nans(x), x_of_nans(x), x_of_nans(x)

    g! = df!_from_df(g, F, inplace)
    fg! = fdf!_from_fdf(fg, F, inplace)
    h! = h!_from_h(h, F, inplace)

    TwiceDifferentiable{T,TG,TH,TX}(f, g!, fg!, nothing, nothing, h!,
                                        copy(F), copy(G), copy(H),
                                        x_f, x_df, x_h,
                                        [0,], [0,], [0,])
end

function TwiceDifferentiable(f, g, h,
                             x::AbstractVector{TX},
                             F::Real = real(zero(eltype(x))),
                             G = alloc_DF(x, F),
                             H = alloc_H(x, F); inplace = true) where {TX}
    g! = df!_from_df(g, F, inplace)
    h! = h!_from_h(h, F, inplace)

    fg! = make_fdf(x, F, f, g!)
    x_f, x_df, x_h = x_of_nans(x), x_of_nans(x), x_of_nans(x)

    return TwiceDifferentiable(f, g!, fg!, nothing, nothing, h!, F, G, H, x_f, x_df, x_h, [0,], [0,], [0,])
end



function TwiceDifferentiable(f, g,
                             x_seed::AbstractVector{T},
                             F::Real = real(zero(T)); autodiff = :finite, inplace = true) where T
    n_x = length(x_seed)

    g! = df!_from_df(g, F, inplace)
    fg! = make_fdf(x_seed, F, f, g!)

    backend = get_adtype(autodiff)
    hess_prep = DI.prepare_hessian(f, backend, x_seed)
    function h!(_h, _x)
        DI.hessian!(f, _h, hess_prep, backend, _x)
        return _h
    end
    TwiceDifferentiable(f, g!, fg!, h!, x_seed, F)
end

TwiceDifferentiable(d::NonDifferentiable, x_seed::AbstractVector{T} = d.x_f, F::Real = real(zero(T)); autodiff = :finite) where {T<:Real} =
    TwiceDifferentiable(d.f, x_seed, F; autodiff = autodiff)

function TwiceDifferentiable(d::OnceDifferentiable, x_seed::AbstractVector{T} = d.x_f,
                             F::Real = real(zero(T)); autodiff = :finite) where T<:Real
    backend = get_adtype(autodiff)
    hess_prep = DI.prepare_hessian(d.f, backend, x_seed)
    function h!(_h, _x)
        DI.hessian!(d.f, _h, hess_prep, backend, _x)
        return _h
    end
    return TwiceDifferentiable(d.f, d.df, d.fdf, h!, x_seed, F, gradient(d))
end

function TwiceDifferentiable(f, x::AbstractArray, F::Real = real(zero(eltype(x)));
                             autodiff = :finite, inplace = true)
    backend = get_adtype(autodiff)
    grad_prep = DI.prepare_gradient(f, backend, x)
    hess_prep = DI.prepare_hessian(f, backend, x)
    function g!(_g, _x)
        DI.gradient!(f, _g, grad_prep, backend, _x)
        return nothing
    end
    function fg!(_g, _x)
        y, _ = DI.value_and_gradient!(f, _g, grad_prep, backend, _x)
        return y
    end
    function h!(_h, _x)
        DI.hessian!(f, _h, hess_prep, backend, _x)
        return _h
    end
    TwiceDifferentiable(f, g!, fg!, h!, x, F)
end

function hv_product!(obj::TwiceDifferentiable, x, v)
    H = hessian!(obj, x)
    return H*v
end

# Used for objectives and solvers where the gradient and Hessian is available/exists
mutable struct TwiceDifferentiableHV{T,TDF,THv,TX} <: AbstractObjective
    f
    fdf
    hv
    F::T
    DF::TDF
    Hv::THv
    x_f::TX
    x_df::TX
    x_hv::TX
    v_hv::TX
    f_calls::Vector{Int}
    df_calls::Vector{Int}
    hv_calls::Vector{Int}
end

# compatibility with old constructor
function TwiceDifferentiableHV(f, fdf, h, x::TX, F::T, G::TG = copy(x), H::THv = copy(x)) where {T, TG, THv, TX}
    x_f, x_df, x_hv, v_hv = x_of_nans(x), x_of_nans(x), x_of_nans(x), x_of_nans(x)
    TwiceDifferentiableHV{T,TG, THv, TX}(f, fdf, h,
                                        copy(F), copy(G), copy(H),
                                        x_f, x_df, x_hv, v_hv,
                                        [0,], [0,], [0,])
end

function TwiceDifferentiableHV(f, fdf, h, x::AbstractVector{T}) where T
    return TwiceDifferentiableHV(f, fdf, h, x, real(zero(T)))
end

# assumption: fg takes (f, g, x)
function TwiceDifferentiableHV(::Nothing, fg, hv, x::AbstractVector{T}) where T
    f   =     x  -> fg(0, nothing, x)
    fg! = (g, x) -> fg(0, g, x)
    return TwiceDifferentiableHV(f, fg!, hv, x)
end

function gradient!!(obj::TwiceDifferentiableHV, x)
    obj.df_calls .+= 1
    copyto!(obj.x_df, x)
    obj.fdf(obj.DF, x)
end

function hv_product!(obj::TwiceDifferentiableHV, x, v)
    if x != obj.x_hv ||  v != obj.v_hv
        hv_product!!(obj, x, v)
    end
    obj.Hv
end
function hv_product!!(obj::TwiceDifferentiableHV, x, v)
    obj.hv_calls .+= 1
    copyto!(obj.x_hv, x)
    copyto!(obj.v_hv, v)
    obj.hv(obj.Hv, x, v)
end
# Deprecate the following?
# Using it makes code non-generic (it requires storage for the result)
hv_product(obj::TwiceDifferentiableHV) = obj.Hv

"""
    f!_from_f(f, F::Abstractarray)

Return an inplace version of f
"""
function f!_from_f(f, F::AbstractArray, inplace)
    if inplace
        return f
    else
        return function ff!(F, x)
            copyto!(F, f(x))
            F
        end
    end
end
function df!_from_df(g, F::Real, inplace)
    if inplace
        return g
    else
        return function gg!(G, x)
            copyto!(G, g(x))
            G
        end
    end
end
function df!_from_df(j, F::AbstractArray, inplace)
    if inplace
        return j
    else
        return function jj!(J, x)
            copyto!(J, j(x))
            J
        end
    end
end
function fdf!_from_fdf(fg, F::Real, inplace)
    if inplace
        return fg
    else
        return function ffgg!(G, x)
            fx, gx = fg(x)
            copyto!(G, gx)
            fx
        end
    end
end
function fdf!_from_fdf(fj, F::AbstractArray, inplace)
    if inplace
        return fj
    else
        return function ffjj!(F, J, x)
            fx, jx = fj(x)
            copyto!(J, jx)
            copyto!(F, fx)
        end
    end
end
function h!_from_h(h, F::Real, inplace)
    if inplace
        return h
    else
        return function hh!(H, x)
            copyto!(H, h(x))
            H
        end
    end
end

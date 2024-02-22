export searchsortednearest
"""
    searchsortednearest(a, x)

Find the index corresponding to the nearest value to `x` in `a`.
"""
function searchsortednearest(a, x)
    i = searchsortedfirst(a, x)
    if i == 1            # x <= a[1]
        return i
    elseif i > length(a) # x > a[end]
        return length(a)
    elseif a[i] == x     # x is one of the elements
        return i
    else                 # general case
        return abs(a[i] - x) < abs(a[i - 1] - x) ? i : i - 1
    end
end


"""
    linear_interpolation(indep_vars, dep_vars, indep_value)

Carries out linear interpolation to obtain a value at
location `indep_value`, using a independent variable
1-d vector `indep_vars` and a dependent variable 
1-d vector `dep_vars`.

If the `indep_value` is outside the range of `indep_vars`, this
returns the endpoint value closest.
"""
function linear_interpolation(indep_vars, dep_vars, indep_value)
    N = length(indep_vars)
    id = searchsortedfirst(indep_vars, indep_value)
    if id == 1
        dep_vars[begin]
    elseif id == N + 1
        dep_vars[end]
    else
        id_prev = id - 1
        x0, x1 = indep_vars[id_prev], indep_vars[id]
        y0, y1 = dep_vars[id_prev], dep_vars[id]
        y0 + (y1 - y0) / (x1 - x0) * (indep_value - x0)
    end
end

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

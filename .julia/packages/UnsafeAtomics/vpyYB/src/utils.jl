if !@isdefined(⊼)
    ⊼(a, b) = ~(a & b)
end

if !@isdefined(LazyString)
    const LazyString = string
end

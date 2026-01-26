function sum(x::Vector)
    acc = zero(eltype(x))
    for el in x
        acc = acc + el
    end
    return acc
end

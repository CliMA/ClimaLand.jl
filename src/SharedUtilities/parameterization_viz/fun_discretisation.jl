function mat(x, y, r, fun, params)    
    x = collect(range(x[1], length=r, stop=x[2])) 
    y = collect(range(y[1], length=r, stop=y[2])) 
    X = repeat(1:r, inner=r)  
    Y = repeat(1:r, outer=r) 
    X2 = repeat(x, inner=r)    
    Y2 = repeat(y, outer=r) 
    FMatrix = Matrix(sparse(X, Y, fun.(FT.(X2), FT.(Y2), repeat([params], r*r))))
    return x, y, FMatrix
end

function d1_vec(x, y, fun, params)
    vec = fun.(FT.(x), FT.(repeat([y], 31)), repeat([params], 31))
    return vec
end

function d2_vec(x, y, fun, params)
    vecM = fun.(FT.(repeat([x], 31)), FT.(y), repeat([params], 31))
    return vecM
end

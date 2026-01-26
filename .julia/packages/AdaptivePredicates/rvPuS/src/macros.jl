@inline function Absolute(a)
    return abs(a)
end

@inline function Fast_Two_Sum_Tail(a, b, x)
    bvirt = x - a
    y = b - bvirt
    return y
end

@inline function Fast_Two_Sum(a, b)
    x = a + b
    y = Fast_Two_Sum_Tail(a, b, x)
    return x, y
end

@inline function Fast_Two_Diff_Tail(a, b, x)
    bvirt = a - x
    y = bvirt - b
    return y
end

@inline function Fast_Two_Diff(a, b)
    x = a - b
    y = Fast_Two_Diff_Tail(a, b, x)
    return x, y
end

@inline function Two_Sum_Tail(a, b, x)
    bvirt = x - a
    avirt = x - bvirt
    bround = b - bvirt
    around = a - avirt
    y = around + bround
    return y
end

@inline function Two_Sum(a, b)
    x = a + b
    y = Two_Sum_Tail(a, b, x)
    return x, y
end

@inline function Two_Diff_Tail(a, b, x)
    bvirt = a - x
    avirt = x + bvirt
    bround = bvirt - b
    around = a - avirt
    y = around + bround
    return y
end

@inline function Two_Diff(a, b)
    x = a - b
    y = Two_Diff_Tail(a, b, x)
    return x, y
end

@inline function Split(a)
    c = splitter(eltype(a)) * a
    abig = c - a
    ahi = c - abig
    alo = a - ahi
    return ahi, alo
end

@inline function Two_Product_Tail(a, b, x)
    ahi, alo = Split(a)
    bhi, blo = Split(b)
    err1 = x - (ahi * bhi)
    err2 = err1 - (alo * bhi)
    err3 = err2 - (ahi * blo)
    y = (alo * blo) - err3
    return y
end

@inline function Two_Product(a, b)
    x = a * b
    y = Two_Product_Tail(a, b, x)
    return x, y
end

@inline function Two_Product_Presplit(a, b, bhi, blo)
    x = a * b
    ahi, alo = Split(a)
    err1 = x - (ahi * bhi)
    err2 = err1 - (alo * bhi)
    err3 = err2 - (ahi * blo)
    y = (alo * blo) - err3
    return x, y
end

@inline function Two_Product_2Presplit(a, ahi, alo, b, bhi, blo)
    x = a * b
    err1 = x - (ahi * bhi)
    err2 = err1 - (alo * bhi)
    err3 = err2 - (ahi * blo)
    y = (alo * blo) - err3
    return x, y
end

@inline function Square_Tail(a, x)
    ahi, alo = Split(a)
    err1 = x - (ahi * ahi)
    err3 = err1 - ((ahi + ahi) * alo)
    y = (alo * alo) - err3
    return y
end

@inline function Square(a)
    x = a * a
    y = Square_Tail(a, x)
    return x, y
end

@inline function Two_One_Sum(a1, a0, b)
    _i, x0 = Two_Sum(a0, b)
    x2, x1 = Two_Sum(a1, _i)
    return x2, x1, x0
end

@inline function Two_One_Diff(a1, a0, b)
    _i, x0 = Two_Diff(a0, b)
    x2, x1 = Two_Sum(a1, _i)
    return x2, x1, x0
end

@inline function Two_Two_Sum(a1, a0, b1, b0)
    _j, _0, x0 = Two_One_Sum(a1, a0, b0)
    x3, x2, x1 = Two_One_Sum(_j, _0, b1)
    return x3, x2, x1, x0
end

@inline function Two_Two_Diff(a1, a0, b1, b0)
    _j, _0, x0 = Two_One_Diff(a1, a0, b0)
    x3, x2, x1 = Two_One_Diff(_j, _0, b1)
    return x3, x2, x1, x0
end

@inline function Four_One_Sum(a3, a2, a1, a0, b)
    _j, x1, x0 = Two_One_Sum(a1, a0, b)
    x4, x3, x2 = Two_One_Sum(a3, a2, _j)
    return x4, x3, x2, x1, x0
end

@inline function Four_Two_Sum(a3, a2, a1, a0, b1, b0)
    _k, _2, _1, _0, x0 = Four_One_Sum(a3, a2, a1, a0, b0)
    x5, x4, x3, x2, x1 = Four_One_Sum(_k, _2, _1, _0, b1)
    return x5, x4, x3, x2, x1, x0
end

@inline function Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0)
    _l, _2, _1, _0, x1, x0 = Four_Two_Sum(a3, a2, a1, a0, b1, b0)
    x7, x6, x5, x4, x3, x2 = Four_Two_Sum(_l, _2, _1, _0, b4, b3)
    return x7, x6, x5, x4, x3, x2, x1, x0
end

@inline function Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b)
    _j, x3, x2, x1, x0 = Four_One_Sum(a3, a2, a1, a0, b)
    x8, x7, x6, x5, x4 = Four_One_Sum(a7, a6, a5, a4, _j)
    return x8, x7, x6, x5, x4, x3, x2, x1, x0
end

@inline function Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0)
    _k, _6, _5, _4, _3, _2, _1, _0, x0 = Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0)
    x9, x8, x7, x6, x5, x4, x3, x2, x1 = Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1)
    return x9, x8, x7, x6, x5, x4, x3, x2, x1, x0
end

@inline function Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0)
    _l, _6, _5, _4, _3, _2, _1, _0, x1, x0 = Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0)
    x11, x10, x9, x8, x7, x6, x5, x4, x3, x2 = Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3)
    return x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0
end

@inline function Two_One_Product(a1, a0, b)
    bhi, blo = Split(b)
    _i, x0 = Two_Product_Presplit(a0, b, bhi, blo)
    _j, _0 = Two_Product_Presplit(a1, b, bhi, blo)
    _k, x1 = Two_Sum(_i, _0)
    x3, x2 = Fast_Two_Sum(_j, _k)
    return x3, x2, x1, x0
end

@inline function Four_One_Product(a3, a2, a1, a0, b)
    bhi, blo = Split(b)
    _i, x0 = Two_Product_Presplit(a0, b, bhi, blo)
    _j, _0 = Two_Product_Presplit(a1, b, bhi, blo)
    _k, x1 = Two_Sum(_i, _0)
    _i, x2 = Fast_Two_Sum(_j, _k)
    _j, _0 = Two_Product_Presplit(a2, b, bhi, blo)
    _k, x3 = Two_Sum(_i, _0)
    _i, x4 = Fast_Two_Sum(_j, _k)
    _j, _0 = Two_Product_Presplit(a3, b, bhi, blo)
    _k, x5 = Two_Sum(_i, _0)
    x7, x6 = Fast_Two_Sum(_j, _k)
    return x7, x6, x5, x4, x3, x2, x1, x0
end

@inline function Two_Two_Product(a1, a0, b1, b0)
    a0hi, a0lo = Split(a0)
    bhi, blo = Split(b0)
    _i, x0 = Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo)
    a1hi, a1lo = Split(a1)
    _j, _0 = Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo)
    _k, _1 = Two_Sum(_i, _0)
    _l, _2 = Fast_Two_Sum(_j, _k)
    bhi, blo = Split(b1)
    _i, _0 = Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo)
    _k, x1 = Two_Sum(_1, _0)
    _j, _1 = Two_Sum(_2, _k)
    _m, _2 = Two_Sum(_l, _j)
    _j, _0 = Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo)
    _n, _0 = Two_Sum(_i, _0)
    _i, x2 = Two_Sum(_1, _0)
    _k, _1 = Two_Sum(_2, _i)
    _l, _2 = Two_Sum(_m, _k)
    _k, _0 = Two_Sum(_j, _n)
    _j, x3 = Two_Sum(_1, _0)
    _i, _1 = Two_Sum(_2, _j)
    _m, _2 = Two_Sum(_l, _i)
    _i, x4 = Two_Sum(_1, _k)
    _k, x5 = Two_Sum(_2, _i)
    x7, x6 = Two_Sum(_m, _k)
    return x7, x6, x5, x4, x3, x2, x1, x0
end

@inline function Two_Square(a1, a0)
    _j, x0 = Square(a0)
    _0 = a0 + a0
    _k, _1 = Two_Product(a1, _0)
    _l, _2, x1 = Two_One_Sum(_k, _1, _j)
    _j, _1 = Square(a1)
    x5, x4, x3, x2 = Two_Two_Sum(_j, _1, _l, _2)
    return x5, x4, x3, x2, x1, x0
end
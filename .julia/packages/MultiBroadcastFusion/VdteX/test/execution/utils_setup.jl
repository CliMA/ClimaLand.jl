get_array(AType, FT, s) = AType(zeros(FT, s...))

function get_arrays(sym, AType, FT, s, n = 4)
    println("array_type = $AType")
    fn = ntuple(i -> Symbol(sym, i), n)
    return (; zip(fn, ntuple(_ -> get_array(AType, FT, s), n))...)
end

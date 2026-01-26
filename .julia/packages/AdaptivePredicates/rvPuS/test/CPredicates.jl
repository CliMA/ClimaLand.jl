module CPredicates
const C64 = joinpath(@__DIR__, "original", "libpredicates64.so")
const C32 = joinpath(@__DIR__, "original", "libpredicates32.so")
const F64 = Float64
const F32 = Float32

for (C, T, F) in ((C64, F64, Cdouble), (C32, F32, Cfloat))
    @eval begin
        function Absolute(a::$T)
            return @ccall $C._Absolute(a::$F)::$F
        end

        function Fast_Two_Sum_Tail(a::$T, b::$T, x::$T)
            return @ccall $C._Fast_Two_Sum_Tail(a::$F, b::$F, x::$F)::$F
        end

        function Fast_Two_Sum(a::$T, b::$T)
            return @ccall $C._Fast_Two_Sum(a::$F, b::$F)::NTuple{2,$F}
        end

        function Fast_Two_Diff_Tail(a::$T, b::$T, x::$T)
            return @ccall $C._Fast_Two_Diff_Tail(a::$F, b::$F, x::$F)::$F
        end

        function Fast_Two_Diff(a::$T, b::$T)
            return @ccall $C._Fast_Two_Diff(a::$F, b::$F)::NTuple{2,$F}
        end

        function Two_Sum_Tail(a::$T, b::$T, x::$T)
            return @ccall $C._Two_Sum_Tail(a::$F, b::$F, x::$F)::$F
        end

        function Two_Sum(a::$T, b::$T)
            return @ccall $C._Two_Sum(a::$F, b::$F)::NTuple{2,$F}
        end

        function Two_Diff_Tail(a::$T, b::$T, x::$T)
            return @ccall $C._Two_Diff_Tail(a::$F, b::$F, x::$F)::$F
        end

        function Two_Diff(a::$T, b::$T)
            return @ccall $C._Two_Diff(a::$F, b::$F)::NTuple{2,$F}
        end

        function Split(a::$T)
            return @ccall $C._Split(a::$F)::NTuple{2,$F}
        end

        function Two_Product_Tail(a::$T, b::$T, x::$T)
            return @ccall $C._Two_Product_Tail(a::$F, b::$F, x::$F)::$F
        end

        function Two_Product(a::$T, b::$T)
            return @ccall $C._Two_Product(a::$F, b::$F)::NTuple{2,$F}
        end

        function Two_Product_Presplit(a::$T, b::$T, bhi::$T, blo::$T)
            return @ccall $C._Two_Product_Presplit(a::$F, b::$F, bhi::$F, blo::$F)::NTuple{2,$F}
        end

        function Two_Product_2Presplit(a::$T, ahi::$T, alo::$T, b::$T, bhi::$T, blo::$T)
            return @ccall $C._Two_Product_2Presplit(a::$F, ahi::$F, alo::$F, b::$F, bhi::$F, blo::$F)::NTuple{2,$F}
        end

        function Square_Tail(a::$T, x::$T)
            return @ccall $C._Square_Tail(a::$F, x::$F)::$F
        end

        function Square(a::$T)
            return @ccall $C._Square(a::$F)::NTuple{2,$F}
        end

        function Two_One_Sum(a1::$T, a0::$T, b::$T)
            return @ccall $C._Two_One_Sum(a1::$F, a0::$F, b::$F)::NTuple{3,$F}
        end

        function Two_One_Diff(a1::$T, a0::$T, b::$T)
            return @ccall $C._Two_One_Diff(a1::$F, a0::$F, b::$F)::NTuple{3,$F}
        end

        function Two_Two_Sum(a1::$T, a0::$T, b1::$T, b0::$T)
            return @ccall $C._Two_Two_Sum(a1::$F, a0::$F, b1::$F, b0::$F)::NTuple{4,$F}
        end

        function Two_Two_Diff(a1::$T, a0::$T, b1::$T, b0::$T)
            return @ccall $C._Two_Two_Diff(a1::$F, a0::$F, b1::$F, b0::$F)::NTuple{4,$F}
        end

        function Four_One_Sum(a3::$T, a2::$T, a1::$T, a0::$T, b::$T)
            return @ccall $C._Four_One_Sum(a3::$F, a2::$F, a1::$F, a0::$F, b::$F)::NTuple{5,$F}
        end

        function Four_Two_Sum(a3::$T, a2::$T, a1::$T, a0::$T, b1::$T, b0::$T)
            return @ccall $C._Four_Two_Sum(a3::$F, a2::$F, a1::$F, a0::$F, b1::$F, b0::$F)::NTuple{6,$F}
        end

        function Four_Four_Sum(a3::$T, a2::$T, a1::$T, a0::$T, b3::$T, b2::$T, b1::$T, b0::$T)
            return @ccall $C._Four_Four_Sum(a3::$F, a2::$F, a1::$F, a0::$F, b3::$F, b2::$F, b1::$F, b0::$F)::NTuple{8,$F}
        end

        function Eight_One_Sum(a7::$T, a6::$T, a5::$T, a4::$T, a3::$T, a2::$T, a1::$T, a0::$T, b::$T)
            return @ccall $C._Eight_One_Sum(a7::$F, a6::$F, a5::$F, a4::$F, a3::$F, a2::$F, a1::$F, a0::$F, b::$F)::NTuple{9,$F}
        end

        function Eight_Two_Sum(a7::$F, a6::$F, a5::$F, a4::$F, a3::$T, a2::$T, a1::$T, a0::$T, b1::$T, b0::$T)
            return @ccall $C._Eight_Two_Sum(a7::$F, a6::$F, a5::$F, a4::$F, a3::$F, a2::$F, a1::$F, a0::$F, b1::$F, b0::$F)::NTuple{10,$F}
        end

        function Eight_Four_Sum(a7::$T, a6::$T, a5::$T, a4::$T, a3::$T, a2::$T, a1::$T, a0::$T, b3::$T, b2::$T, b1::$T, b0::$T)
            return @ccall $C._Eight_Four_Sum(a7::$F, a6::$F, a5::$F, a4::$F, a3::$F, a2::$F, a1::$F, a0::$F, b3::$F, b2::$F, b1::$F, b0::$F)::NTuple{12,$F}
        end

        function Two_One_Product(a1::$T, a0::$T, b::$T)
            return @ccall $C._Two_One_Product(a1::$F, a0::$F, b::$F)::NTuple{4,$F}
        end

        function Four_One_Product(a3::$T, a2::$T, a1::$T, a0::$T, b::$T)
            return @ccall $C._Four_One_Product(a3::$F, a2::$F, a1::$F, a0::$F, b::$F)::NTuple{8,$F}
        end

        function Two_Two_Product(a1::$T, a0::$T, b1::$T, b0::$T)
            return @ccall $C._Two_Two_Product(a1::$F, a0::$F, b1::$F, b0::$F)::NTuple{8,$F}
        end

        function Two_Square(a1::$T, a0::$T)
            return @ccall $C._Two_Square(a1::$F, a0::$F)::NTuple{6,$F}
        end

        function grow_expansion(elen::Int, e::Vector{$T}, b::$T, h::Vector{$T})
            hindex = @ccall $C.grow_expansion(elen::Int, e::Ptr{$F}, b::$F, h::Ptr{$F})::Int
            return h, hindex
        end

        function grow_expansion_zeroelim(elen::Int, e::Vector{$T}, b::$T, h::Vector{$T})
            hindex = @ccall $C.grow_expansion_zeroelim(elen::Int, e::Ptr{$F}, b::$F, h::Ptr{$F})::Int
            return h, hindex
        end

        function expansion_sum(elen::Int, e::Vector{$T}, flen::Int, f::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.expansion_sum(elen::Int, e::Ptr{$F}, flen::Int, f::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function expansion_sum_zeroelim1(elen::Int, e::Vector{$T}, flen::Int, f::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.expansion_sum_zeroelim1(elen::Int, e::Ptr{$F}, flen::Int, f::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function expansion_sum_zeroelim2(elen::Int, e::Vector{$T}, flen::Int, f::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.expansion_sum_zeroelim2(elen::Int, e::Ptr{$F}, flen::Int, f::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function fast_expansion_sum(elen::Int, e::Vector{$T}, flen::Int, f::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.fast_expansion_sum(elen::Int, e::Ptr{$F}, flen::Int, f::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function fast_expansion_sum_zeroelim(elen::Int, e::Vector{$T}, flen::Int, f::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.fast_expansion_sum_zeroelim(elen::Int, e::Ptr{$F}, flen::Int, f::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function linear_expansion_sum(elen::Int, e::Vector{$T}, flen::Int, f::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.linear_expansion_sum(elen::Int, e::Ptr{$F}, flen::Int, f::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function linear_expansion_sum_zeroelim(elen::Int, e::Vector{$T}, flen::Int, f::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.linear_expansion_sum_zeroelim(elen::Int, e::Ptr{$F}, flen::Int, f::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function scale_expansion(elen::Int, e::Vector{$T}, b::$T, h::Vector{$T})
            hindex = @ccall $C.scale_expansion(elen::Int, e::Ptr{$F}, b::$F, h::Ptr{$F})::Int
            return h, hindex
        end

        function scale_expansion_zeroelim(elen::Int, e::Vector{$T}, b::$T, h::Vector{$T})
            hindex = @ccall $C.scale_expansion_zeroelim(elen::Int, e::Ptr{$F}, b::$F, h::Ptr{$F})::Int
            return h, hindex
        end

        function compress(elen::Int, e::Vector{$T}, h::Vector{$T})
            hindex = @ccall $C.compress(elen::Int, e::Ptr{$F}, h::Ptr{$F})::Int
            return h, hindex
        end

        function estimate(elen::Int, e::Vector{$T})
            h = @ccall $C.estimate(elen::Int, e::Ptr{$F})::$F
            return h
        end
    end

    for post in (:fast, :exact, :slow, Symbol(""))
        @eval begin
            function $(Symbol("orient2d$post"))(pa::NTuple{2,$T}, pb::NTuple{2,$T}, pc::NTuple{2,$T})
                return @ccall $C.$(Symbol("orient2d$post"))(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F})::$F
            end

            function $(Symbol("orient3d$post"))(pa::NTuple{3,$T}, pb::NTuple{3,$T}, pc::NTuple{3,$T}, pd::NTuple{3,$T})
                return @ccall $C.$(Symbol("orient3d$post"))(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F}, Ref(pd)::Ptr{$F})::$F
            end

            function $(Symbol("incircle$post"))(pa::NTuple{2,$T}, pb::NTuple{2,$T}, pc::NTuple{2,$T}, pd::NTuple{2,$T})
                return @ccall $C.$(Symbol("incircle$post"))(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F}, Ref(pd)::Ptr{$F})::$F
            end

            function $(Symbol("insphere$post"))(pa::NTuple{3,$T}, pb::NTuple{3,$T}, pc::NTuple{3,$T}, pd::NTuple{3,$T}, pe::NTuple{3,$T})
                return @ccall $C.$(Symbol("insphere$post"))(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F}, Ref(pd)::Ptr{$F}, Ref(pe)::Ptr{$F})::$F
            end
        end
    end

    @eval begin
        function orient2dadapt(pa::NTuple{2,$T}, pb::NTuple{2,$T}, pc::NTuple{2,$T}, detsum::$T)
            return @ccall $C.orient2dadapt(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F}, detsum::$F)::$F
        end

        function orient3dadapt(pa::NTuple{3,$T}, pb::NTuple{3,$T}, pc::NTuple{3,$T}, pd::NTuple{3,$T}, permanent::$T)
            return @ccall $C.orient3dadapt(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F}, Ref(pd)::Ptr{$F}, permanent::$F)::$F
        end

        function incircleadapt(pa::NTuple{2,$T}, pb::NTuple{2,$T}, pc::NTuple{2,$T}, pd::NTuple{2,$T}, permanent::$T)
            return @ccall $C.incircleadapt(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F}, Ref(pd)::Ptr{$F}, permanent::$F)::$F
        end

        function insphereadapt(pa::NTuple{3,$T}, pb::NTuple{3,$T}, pc::NTuple{3,$T}, pd::NTuple{3,$T}, pe::NTuple{3,$T}, permanent::$T)
            return @ccall $C.insphereadapt(Ref(pa)::Ptr{$F}, Ref(pb)::Ptr{$F}, Ref(pc)::Ptr{$F}, Ref(pd)::Ptr{$F}, Ref(pe)::Ptr{$F}, permanent::$F)::$F
        end
    end
end

const C64_consts = @ccall C64.exactinit2()::NTuple{15,Float64}
const C32_consts = @ccall C32.exactinit2()::NTuple{15,Float32}

function __init__()
    @ccall C64.exactinit()::Cvoid
    @ccall C32.exactinit()::Cvoid
end

end # module CPredicates

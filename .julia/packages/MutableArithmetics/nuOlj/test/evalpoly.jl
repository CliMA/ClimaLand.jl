# Copyright (c) 2023 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

abstract type OpSignature end
struct RegularSignature <: OpSignature end
struct ToSignature <: OpSignature end

function signature_type_of(
    ::Union{typeof(MA.operate),typeof(MA.operate!),typeof(MA.operate!!)},
)
    return RegularSignature()
end

function signature_type_of(
    ::Union{typeof(MA.operate_to!),typeof(MA.operate_to!!)},
)
    return ToSignature()
end

function op_may_modify_first_argument(op::Function)
    return (op != MA.operate) & (signature_type_of(op) == RegularSignature())
end

function op_is_allowed_for_arithmetic(
    ::typeof(evalpoly),
    ::Union{typeof(MA.operate!),typeof(MA.operate_to!)},
    ::Type{F},
) where {F<:Any}
    return MA.mutability(F) == MA.IsMutable()
end

function op_is_allowed_for_arithmetic(
    ::typeof(evalpoly),
    ::Union{typeof(MA.operate),typeof(MA.operate!!),typeof(MA.operate_to!!)},
    ::Type{F},
) where {F<:Any}
    return true
end

function allowed_allocated_byte_count(
    ::typeof(evalpoly),
    ::Union{typeof(MA.operate),typeof(MA.operate!),typeof(MA.operate!!)},
    ::Type{F},
) where {F<:Any}
    return @allocated zero(F)
end

function allowed_allocated_byte_count(
    ::typeof(evalpoly),
    ::Union{typeof(MA.operate_to!),typeof(MA.operate_to!!)},
    ::Type{F},
) where {F<:Any}
    return @allocated nothing
end

const evalpoly_supported_arithmetics =
    (BigFloat, Rational{Int}, Float64, Float32, Int)

const evalpoly_operations =
    (MA.operate, MA.operate!, MA.operate_to!, MA.operate!!, MA.operate_to!!)

@testset "evalpoly with $op and $F" for F in evalpoly_supported_arithmetics,
    op in evalpoly_operations

    op_is_allowed_for_arithmetic(evalpoly, op, F) || continue

    sig = signature_type_of(op)

    if MA.mutability(F) == MA.IsMutable()
        @testset "empty coefficients $collection_type" for collection_type in
                                                           ("vector", "tuple")
            out = one(F)
            x = one(F)
            backup = MA.copy_if_mutable(x)
            coefs = F[]
            if collection_type == "tuple"
                coefs = ()
            end

            # Check that the result value is OK
            if sig == RegularSignature()
                @test iszero(@inferred op(evalpoly, x, coefs))
            elseif sig == ToSignature()
                @test iszero(@inferred op(out, evalpoly, x, coefs))
            else
                error("unexpected")
            end

            # Check that the arguments are unmodified
            op_may_modify_first_argument(op) || @test backup == x
        end
    end

    @testset "exact values: $degree, $x_int" for degree in 0:4, x_int in -5:5
        coefs_int = rand(-6:6, degree + 1)

        coefs = map(F, coefs_int)
        x = F(x_int)

        reference = evalpoly(x, coefs)

        @test reference == evalpoly(x_int, coefs_int)

        @testset "collection type $collection_type" for collection_type in
                                                        ("vector", "tuple")
            out = zero(F)
            if collection_type == "vector"
                coefs_arg = coefs
            else
                coefs_arg = (coefs...,)
            end
            if sig == RegularSignature()
                out = @inferred op(evalpoly, x, coefs_arg)
            elseif sig == ToSignature()
                out = @inferred op(out, evalpoly, x, coefs_arg)
            else
                error("unexpected")
            end

            # Check that the result value is OK
            @test reference == out

            # Check that the arguments are unmodified
            if op_may_modify_first_argument(op)
                x = F(x_int)
            else
                @test x == F(x_int)
            end
            @test coefs == map(F, coefs_int)
        end
    end

    @testset "allocation" begin
        byte_cnt = allowed_allocated_byte_count(evalpoly, op, F)
        coefs_tuple = map(F, (0, 1, 0, 1, 1))
        @testset "collection type $collection_type" for collection_type in
                                                        ("vector", "tuple")
            if collection_type == "vector"
                coefs = collect(coefs_tuple)
            else
                coefs = coefs_tuple
            end
            local tested_fun
            if sig == RegularSignature()
                tested_fun = let x = one(F), coefs = coefs
                    () -> op(evalpoly, x, coefs)
                end
            elseif sig == ToSignature()
                tested_fun = let o = one(F), x = one(F), coefs = coefs
                    () -> op(o, evalpoly, x, coefs)
                end
            else
                error("unexpected")
            end
            alloc_test(tested_fun, byte_cnt)
        end
    end
end

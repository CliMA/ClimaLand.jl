# This file is a part of InverseFunctions.jl, licensed under the MIT License (MIT).

using Test
using InverseFunctions
using Unitful
using Dates


foo(x) = inv(exp(-x) + 1)
inv_foo(y) = log(y / (1 - y))

InverseFunctions.inverse(::typeof(foo)) = inv_foo
InverseFunctions.inverse(::typeof(inv_foo)) = foo


struct Bar{MT<:AbstractMatrix}
    A::MT
end

(f::Bar)(x) = f.A * x
InverseFunctions.inverse(f::Bar) = Bar(inv(f.A))

@static if VERSION >= v"1.6"
    _bc_func(f) = Base.Broadcast.BroadcastFunction(f)
else
    _bc_func(f) = Base.Fix1(broadcast, f)
end

@testset "inverse" begin
    f_without_inverse(x) = 1
    @test inverse(f_without_inverse) isa NoInverse
    @test_throws ErrorException inverse(f_without_inverse)(42)
    @test inverse(inverse(f_without_inverse)) === f_without_inverse

    for f in (f_without_inverse ∘ exp, exp ∘ f_without_inverse, _bc_func(f_without_inverse), Base.Fix1(broadcast, f_without_inverse), Base.Fix1(map, f_without_inverse))
        @test inverse(f) == NoInverse(f)
        @test inverse(inverse(f)) == f
    end

    @test @inferred(inverse(Complex)) isa NoInverse{Type{Complex}}
    @test @inferred(NoInverse(Complex)) isa NoInverse{Type{Complex}}

    InverseFunctions.test_inverse(inverse, log, compare = ===)
end

@testset "maths" begin
    InverseFunctions.test_inverse(!, false)

    x = rand()
    for f in (
            foo, inv_foo, log, log2, log10, log1p, sqrt,
            Base.Fix2(^, 3*rand() - 0.5), Base.Fix2(^, rand(float.([-10:-1; 1:10]))), Base.Fix1(^, rand()), Base.Fix1(log, rand()), Base.Fix1(log, 1/rand()), Base.Fix2(log, rand()),
        )
        InverseFunctions.test_inverse(f, x)
    end
    for f in (
            +, -, exp, exp2, exp10, expm1, cbrt, deg2rad, rad2deg, conj,
            sinh, tanh, coth, csch, asinh, atanh, acsch,  # all invertible hyperbolic functions aside from acoth
            Base.Fix1(+, rand()), Base.Fix2(+, rand()), Base.Fix1(-, rand()), Base.Fix2(-, rand()),
            Base.Fix1(*, rand()), Base.Fix2(*, rand()), Base.Fix1(/, rand()), Base.Fix2(/, rand()), Base.Fix1(\, rand()), Base.Fix2(\, rand()),
            Base.Fix2(^, rand(-11:2:11)),
        )
        InverseFunctions.test_inverse(f, x)
        InverseFunctions.test_inverse(f, -x)
    end
    # acoth only defined for |x| > 1
    InverseFunctions.test_inverse(acoth, 1 + x)
    InverseFunctions.test_inverse(acoth, -1 - x)

    InverseFunctions.test_inverse(conj, 2 - 3im)
    InverseFunctions.test_inverse(reverse, [10, 20, 30])

    x = rand(0:10)
    for f in (Base.Fix2(divrem, rand([-5:-1; 1:5])), Base.Fix2(fldmod, rand([-5:-1; 1:5])), Base.Fix2(divrem, 0.123), Base.Fix2(fldmod, 0.123))
        compare = (a, b) -> all(isapprox.(a, b))
        InverseFunctions.test_inverse(f, x; compare=compare)
        InverseFunctions.test_inverse(f, -x; compare=compare)
        InverseFunctions.test_inverse(f, x/9; compare=compare)
        InverseFunctions.test_inverse(f, -x/9; compare=compare)
    end

    # ensure that inverses have domains compatible with original functions
    @test_throws DomainError inverse(sqrt)(-1.0)
    InverseFunctions.test_inverse(sqrt, complex(-1.0))
    InverseFunctions.test_inverse(sqrt, complex(1.0))
    @test_throws DomainError inverse(Base.Fix1(*, 0))
    @test_throws DomainError inverse(Base.Fix2(^, 0))
    @test_throws DomainError inverse(Base.Fix1(log, -2))(5)
    @test_throws DomainError inverse(Base.Fix1(log, 2))(-5)
    InverseFunctions.test_inverse(inverse(Base.Fix1(log, 2)), complex(-5))
    @test_throws DomainError inverse(Base.Fix2(^, 0.5))(-5)
    @test_throws DomainError inverse(Base.Fix2(^, 0.51))(complex(-5))
    @test_throws DomainError inverse(Base.Fix2(^, 2))(complex(-5))
    InverseFunctions.test_inverse(Base.Fix2(^, 0.5), complex(-5))
    @test_throws DomainError inverse(Base.Fix2(^, 2))
    @test_throws DomainError inverse(Base.Fix2(^, -4))
    InverseFunctions.test_inverse(Base.Fix2(^, 2.0), 4)
    @test_throws DomainError inverse(Base.Fix1(^, 2.0))(-4)
    @test_throws DomainError inverse(Base.Fix1(^, -2.0))(4)
    @test_throws DomainError inverse(Base.Fix1(^, 0))(4)
    @test_throws DomainError inverse(Base.Fix1(log, -2))(4)
    @test_throws DomainError inverse(Base.Fix1(log, 1))(4)
    @test_throws DomainError inverse(Base.Fix2(^, 0))(4)
    @test_throws DomainError inverse(Base.Fix2(log, -2))(4)
    @test_throws DomainError inverse(Base.Fix2(log, 1))(4)
    InverseFunctions.test_inverse(Base.Fix2(^, -1), complex(-5.))
    @test_throws DomainError inverse(Base.Fix2(^, 2))(-5)
    @test_throws DomainError inverse(Base.Fix1(^, 2))(-5)
    @test_throws DomainError inverse(Base.Fix1(^, -2))(3)
    @test_throws DomainError inverse(Base.Fix1(^, -2))(3)

    @test_throws DomainError inverse(Base.Fix2(divrem, 5))((-3, 2))
    @test_throws DomainError inverse(Base.Fix2(fldmod, 5))((-3, -2))
    InverseFunctions.test_inverse(inverse(Base.Fix2(divrem, 5)), (-3, -2); compare=(==))
    InverseFunctions.test_inverse(inverse(Base.Fix2(fldmod, 5)), (-3, 2); compare=(==))

    InverseFunctions.test_inverse(reim, -3; compare=(==))
    InverseFunctions.test_inverse(reim, -3+2im; compare=(==))
    InverseFunctions.test_inverse(Base.splat(complex), (-3, 2); compare=(==))

    A = rand(5, 5)
    for f in (
            identity, inv, adjoint, transpose,
            log, sqrt, +, -, exp,
            Base.Fix1(+, rand(5, 5)), Base.Fix2(+, rand(5, 5)), Base.Fix1(-, rand(5, 5)), Base.Fix2(-, rand(5, 5)),
            Base.Fix1(*, rand()), Base.Fix2(*, rand()), Base.Fix1(*, rand(5, 5)), Base.Fix2(*, rand(5, 5)),
            Base.Fix2(/, rand()), Base.Fix1(/, rand(5, 5)), Base.Fix2(/, rand(5, 5)),
            Base.Fix1(\, rand()), Base.Fix1(\, rand(5, 5)), Base.Fix2(\, rand(5, 5)),
        )
        if f != log || VERSION >= v"1.6"
            # exp(log(A::AbstractMatrix)) ≈ A is broken on at least Julia v1.0
            InverseFunctions.test_inverse(f, A)
        end
    end

    X = rand(5)
    for f in (_bc_func(foo), Base.Fix1(broadcast, foo), Base.Fix1(map, foo))
        for x in (x, fill(x, 3), X)
            InverseFunctions.test_inverse(f, x)
        end
    end

    InverseFunctions.test_inverse(Bar(rand(3,3)), rand(3))

    @static if VERSION >= v"1.6"
        InverseFunctions.test_inverse(log ∘ foo, x)
    end
end

@testset "unitful" begin
    # the majority of inverse just propagate to underlying mathematical functions and don't have any issues with unitful numbers
    # only those that behave treat real numbers differently have to be tested here
    x = rand()u"m"
    InverseFunctions.test_inverse(sqrt, x)
    @test_throws DomainError inverse(sqrt)(-x)

    InverseFunctions.test_inverse(Base.Fix2(^, 3), x)
    InverseFunctions.test_inverse(Base.Fix2(^, 3), -x)
    InverseFunctions.test_inverse(Base.Fix2(^, -3.5), x)
    @test_throws DomainError inverse(Base.Fix2(^, 2))(-x)
end

@testset "dates" begin
    InverseFunctions.test_inverse(Dates.date2epochdays, Date(2020, 1, 2); compare = ===)
    InverseFunctions.test_inverse(Dates.datetime2epochms, DateTime(2020, 1, 2, 12, 34, 56); compare = ===)
    InverseFunctions.test_inverse(Dates.epochdays2date, Int64(1234); compare = ===)
    InverseFunctions.test_inverse(Dates.epochms2datetime, Int64(1234567890); compare = ===)

    InverseFunctions.test_inverse(datetime2unix, DateTime(2020, 1, 2, 12, 34, 56); compare = ===)
    InverseFunctions.test_inverse(unix2datetime, 1234.56; compare = ===)
    InverseFunctions.test_inverse(datetime2julian, DateTime(2020, 1, 2, 12, 34, 56); compare = ===)
    InverseFunctions.test_inverse(julian2datetime, 1234.56; compare = ===)
end

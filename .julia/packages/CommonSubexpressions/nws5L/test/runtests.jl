using Test
using CommonSubexpressions
using CommonSubexpressions: binarize

@testset "basic usage" begin
    f_count = 0
    function f(x)
        f_count += 1
        1
    end

    g_count = 0
    function g(x)
        g_count += 1
        2
    end

    @test f_count == 0
    @test g_count == 0

    x = 1
    @cse [f(g(x)) for i in 1:10]
    @test f_count == 1
    @test g_count == 1

    f_count = 0
    g_count = 0
    [f(g(x)) for i in 1:10]
    @test f_count == 10
    @test g_count == 10

    h_count = 0
    function h(x, y)
        h_count += 1
        3
    end

    @cse function foo(x)
        [h(h(f(g(x)), g(x)), i) == i for i in 1:5]
    end
    f_count = 0
    g_count = 0
    h_count = 0
    @test foo(1) == [false, false, true, false, false]
    @test f_count == 1
    @test g_count == 1
    @test h_count == 6
end

@testset "function definition" begin
    f_count = 0
    function f(x)
        f_count += 1
        1
    end

    @cse function foo(x)
        [f(x) == i for i in 1:5]
    end

    @test foo(1) == [true, false, false, false, false]
    @test f_count == 1
end

@testset "algebra" begin
    H = [1 2; 3 4]
    W = [2 3; 4 5]
    G = [4 5; 6 7]

    @test (@cse inv(H) * G + W * (inv(H) + W)) == begin
        invH = inv(H)
        invH * G + W * (inv(H) + W)
    end
end

@testset "dict" begin
    @test (@cse begin
        x = Dict("foo" => sin(pi), "bar" => 2 + 2)
        x["bar"] = 100
        x["foo"], x["bar"]
    end) == begin
        x = Dict("foo" => sin(pi), "bar" => 2 + 2)
        x["bar"] = 100
        x["foo"], x["bar"]
    end
end

@testset "int" begin
    @test (@cse 1) == 1
end

@testset "nested usage" begin
    f_count = 0
    function f(x)
        f_count += 1
        x
    end

    g_count = 0
    function g(x)
        g_count += 1
        2x
    end

    @test f_count == 0
    @test g_count == 0

    result = @cse for i in 1:5
        x = f(2)
        @cse for j in 1:5
            y = g(i)
        end
    end
    @test f_count == 1
    @test g_count == 5

    @test result == for i in 1:5
        x = f(2)
        for j in 1:5
            y = g(i)
        end
    end
end

@testset "cse function" begin
    expr = cse(:((x + 1) * (x + 1)))
    @test expr.head == :block
    @test length(expr.args) == 3
    @test expr.args[1].head == :(=)
    @test expr.args[2].head == :(=)
    @test expr.args[1].args[2].head == :call
    @test expr.args[1].args[2].args == [:(+), :x, 1]
    @test expr.args[2].args[2].head == :call
    @test expr.args[2].args[2].args[1] == :(*)
end

@testset "warnings" begin
    @test_logs (:warn, "CommonSubexpressions can't yet handle expressions of this form: foo") cse(Expr(:foo, 1, 2, 3))

    @cse(1 + 2 + 3, false)
    @cse(1 + 2 + 3, warn=true)
end

@testset "inplace" begin
    x = 1
    @cse begin
        y = 2 + x
        x += y * x
    end
    @test y == 3
    @test x == 4
end

@testset "binarize" begin
    @test binarize(:a) == :a
    @test binarize(:(a)) == :(a)
    @test binarize(:(+a)) == :(+a)
    @test binarize(:(a + b)) == :(a + b)
    @test binarize(:(a + b + c)) == :((a + b) + c)
    @test binarize(:(a + b + c + d)) == :(((a + b) + c) + d)
    @test binarize(:(a + b + c + d + e)) == :((((a + b) + c) + d) + e)

    # Nested function calls
    @test binarize(:((a + b + c) + (d + e + f))) == :(((a + b) + c) + ((d + e) + f))

    # Arbirary binary functions should be left as-is
    @test binarize(:(f(a + b, c))) == :(f(a + b, c))
    @test binarize(:(f(a + b, g(c, d + e)))) == :(f(a + b, g(c, d + e)))

    # Make sure we can binarize functions whose arguments are expressions
    @test binarize(:(f(a + b, c + d, g(e, f)))) == :(f(f(a + b, c + d), g(e, f)))
end

# Test nesting `@cse` with a macro like `@fastmath`
module NestedMacroTest
    using CommonSubexpressions
    using MacroTools: postwalk, @capture
    using Test

    const special_plus_calls = Ref(0)

    function special_plus(args...)
        special_plus_calls[] += 1
        +(args...)
    end

    macro special_math(expr)
        esc(postwalk(expr) do ex
            if @capture(ex, +(x__))
                :(special_plus($(x...)))
            else
                ex
            end
        end)
    end

    @testset "Nested Macros" begin
        @test(@cse(@special_math(2 + 2)) == 4)
        @test special_plus_calls[] == 1

        special_plus_calls[] = 0
        @test(@cse(@binarize(@special_math(1 + 2 + 3 + 4))) == 10)
        # Three more calls to `special_plus`
        @test special_plus_calls[] == 3

        special_plus_calls[] = 0
        # Test that CSE is actually eliminated common subexpressions even
        # when nested with another macro.
        @test(@cse(@special_math((1 + 2) + (3 + 4))) == 10)
        # Three more calls: 1 + 2, 3 + 4, and 3 + 7
        @test special_plus_calls[] == 3

        special_plus_calls[] = 0
        @test(@cse(@special_math((2 + 2) + (2 + 2))) == 8)
        # Only two more calls because of CSE: 2 + 2, 4 + 4
        @test special_plus_calls[] == 2

        # These particular macros also commute:

        special_plus_calls[] = 0
        @test(@special_math(@cse((1 + 2) + (3 + 4))) == 10)
        @test special_plus_calls[] == 3

        special_plus_calls[] = 0
        @test(@special_math(@cse((2 + 2) + (2 + 2))) == 8)
        @test special_plus_calls[] == 2

        special_plus_calls[] = 0
        @test(@cse(@binarize(@special_math((1 + 2 + 3) + (1 + 2 + 4) + (1 + 2 + 5)))) == 21)
        # Test that the duplicate calls to `1 + 2` were eliminated
        @test special_plus_calls[] == 6
    end
end

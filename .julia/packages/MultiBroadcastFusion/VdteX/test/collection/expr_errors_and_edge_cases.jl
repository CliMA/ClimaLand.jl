#=
using Revise; include(joinpath("test", "collection" "expr_errors_and_edge_cases.jl"))
=#
using Test
import MultiBroadcastFusion as MBF
@testset "loops" begin
    # loops are not allowed, because
    # code transformation occurs at macro
    # expansion time, and we can't generally
    # know how many times the loop will be
    # executed at this time.

    # We could try to specialize on literal ranges, e.g.,
    # `for i in 1:10`, but that is likely an uncommon
    # edge case.

    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        for i in 1:10
            @. y2 = x2 + x3 + x4 + x5
        end
        @. y1 = x1 + x2 + x3 + x4
    end
    @test_throws ErrorException("Loops are not allowed inside fused blocks") MBF.fused_direct(
        expr_in,
    )
end

struct Foo end

@testset "If-statements" begin
    # If-statements are not allowed, because
    # code transformation occurs at macro
    # expansion time, and Bools, even types,
    # are not known at this time.

    # We could specialize on literals, e.g.,
    # `if true`, but that is likely an uncommon
    # edge case.
    foo = Foo()
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        if foo isa Foo
            @. y2 = x2 + x3 + x4 + x5
        end
        @. y1 = x1 + x2 + x3 + x4
    end
    @test_throws ErrorException(
        "If-statements are not allowed inside fused blocks",
    ) MBF.fused_direct(expr_in)
end

bar() = nothing
@testset "Function calls" begin
    # Function calls are not allowed, because
    # this could lead to subtle bugs (order of compute).
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        bar()
        @. y1 = x1 + x2 + x3 + x4
    end
    @test_throws ErrorException(
        "Function calls are not allowed inside fused blocks",
    ) MBF.fused_direct(expr_in)
end

@testset "Non-broadcast variable assignments" begin
    # Non-broadcast variable assignments is not allowed, because
    # this could lead to subtle bugs (order of compute).
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        x1 = x2
        @. y1 = x1 + x2 + x3 + x4
    end
    @test_throws ErrorException(
        "Non-broadcast assignments are not allowed inside fused blocks",
    ) MBF.fused_direct(expr_in)
end

@testset "No let-blocks" begin
    # Let-blocks could hide other non-allowed things
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        let z = 1
        end
        @. y1 = x1 + x2 + x3 + x4
    end
    @test_throws ErrorException(
        "Let-blocks are not allowed inside fused blocks",
    ) MBF.fused_direct(expr_in)
end

@testset "Dangling symbols" begin
    # While inaucuous, we prohibit dangling symbols
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        :a
        @. y1 = x1 + x2 + x3 + x4
    end
    @test_throws ErrorException(
        "Dangling symbols are not allowed inside fused blocks",
    ) MBF.fused_direct(expr_in)
end

@testset "quote" begin
    # Not sure why this would be needed, so
    # we don't allow quotes inside fused blocks.
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        quote end
        @. y1 = x1 + x2 + x3 + x4
    end
    @test_throws ErrorException("Quotes are not allowed inside fused blocks") MBF.fused_direct(
        expr_in,
    )
end

@testset "Comments" begin
    expr_in = quote
        @. y1 = x1 + x2 + x3 + x4
        # Foo bar baz
        # if i in 1:N
        @. y2 = x2 + x3 + x4 + x5
    end

    expr_out = :(tuple(
        Pair(y1, Base.broadcasted(+, x1, x2, x3, x4)),
        Pair(y2, Base.broadcasted(+, x2, x3, x4, x5)),
    ))
    @test MBF.fused_direct(expr_in) == expr_out
end

@testset "Empty" begin
    expr_in = quote end
    @test MBF.fused_direct(expr_in) == :(tuple())
end

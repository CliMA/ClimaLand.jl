using FastClosures, Test

# Test utility wrapping find_var_uses!
function find_var_uses(ex)
    vars = FastClosures.Var[]
    bound_vars = FastClosures.Var[]
    FastClosures.find_var_uses!(vars, bound_vars, ex, 0)
    Symbol[v.name for v in vars]
end

# Check that a particular error occurs during macro invocation.  In 0.7 / 1.0
# this is wrapped in a LoadError, so unwrap that before rethrowing.  See also
# https://discourse.julialang.org/t/exceptions-in-macros-in-julia-0-7-1-0/14145/2
macro test_macro_throws(typ, expr)
    quote
        @test_throws $typ begin
            try
                $(esc(expr))
            catch e
                while e isa LoadError
                    e = e.error
                end
                rethrow(e)
            end
        end
    end
end

# code_warntype issues
function f1()
    if true
    end
    r = 1
    cb = ()->r
    identity(cb)
end

# code_warntype clean
function f2()
    if true
    end
    r = 1
    cb = @closure ()->r
    identity(cb)
end

macro nested_closure(ex)
    quote
        @closure x->x+$(esc(ex))
    end
end

@testset "FastClosures" begin

@testset "Basic closure tests" begin
    @test (@closure ()->1.0)() === 1.0
    a = 0.5
    @test (@closure ()->a)() === a
    @test (@closure (x)->a)(1) === a
    @test (@closure (x,y)->a)(1,2) === a
    @test (@closure (x,y,z)->a)(1,2,3) === a
    @test (@closure (x,y,z)->x+a)(1,2,3) === 1.5
    #@test_broken (@closure (x::Int,y::Int,z::Int)->x+a)(1,2,3) === 1.5
    @test (@closure function blah(x,y,z); x+a; end)(1,2,3) === 1.5
    @test_broken begin
        @closure function blah(x,y,z)
            x+a
        end
        blah(1,2,3)
    end === 1.5
    @test_macro_throws ArgumentError @eval(@closure(1+2))
    # Test that when macroexpand() fails inside the @closure macro, the correct
    # error is generated
    @test_macro_throws UndefVarError @eval(@closure () -> @nonexistent_macro)
end

@testset "@closure use inside a macro" begin
    # Test workaround for
    # https://github.com/JuliaLang/julia/issues/23221
    y = 10
    z = @nested_closure(y)(2)
    @test z == 12
end

@testset "internal find_var_uses torture tests" begin

@testset "Basic usage" begin
    # Basic variable use
    @test find_var_uses(:(a+b*c)) == [:a, :b, :c]
    # Deduplication
    @test find_var_uses(:(a*a)) == [:a]
    # Assignment rebinds, doesn't count as a "use".
    @test find_var_uses(:(a=1)) == []
    # keyword args
    @test find_var_uses(:(f(a, b=1, c=2, d=e))) == [:a, :e]
end


@testset "New bindings" begin
    # Introduce new bindings and use them
    @test find_var_uses(:(
        begin
            a = 1
            b = 2
            c = a+b
        end)) == []
    # `local` qualifies a binding, it's not a variable use
    @test find_var_uses(:(
        begin
            local i=1
            local j
            i+1
        end
        )) == []
    # Declaring variables of given type
    @test find_var_uses(:(
        begin
            y::S = 0.0
            x::T
            z
        end
        )) == [:S, :T, :z]
end


@testset "Scopes" begin
    @test find_var_uses(:(
        begin
            b = 1
            for i=1:10
                c = b   # uses b as bound above
            end
        end)) == []
    @test find_var_uses(:(
        begin
            a = 1
            for i=1:10
                b = a
            end
            c = b # uses b bound from outside the expression
        end)) == [:b]
    @test find_var_uses(:(
        try
            a = 1+c
        catch err
            show(err)
        finally
            b = 1+d
        end)) == [:c,:d]
end


@testset "Assignment and tuple unpacking" begin
    @test find_var_uses(:(
        begin
            a[i] = 10
        end
        )) == [:a,:i]
    @test find_var_uses(:(
        begin
            a,b = 10,11
            (c,(d,e[f])) = (1,(2,3))
            a+b
        end
        )) == [:e,:f]
end

@testset "Excluded expr heads" begin
    # Expr(:inbounds) excluded
    @test find_var_uses(macroexpand(@__MODULE__, quote
        b = @inbounds a[i]
    end)) == [:a, :i]

    # Expr(:quote) excluded
    @test find_var_uses(quote
        b = :(x + y)
    end) == []
end

@testset "do syntax" begin
    function test_do(a)
        b = zeros(10)
        c = (rand(),rand())
        c = @closure ntuple(Val{2}()) do i
            c[i] + b[i] + a
        end
    end
    @inferred test_do(22)
end

end


end

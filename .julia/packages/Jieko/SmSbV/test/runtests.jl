using Test
using REPL # need to load this module to make docstring rendering behavior consistent after v1.11
using Jieko: Jieko

module Example
include("example/basic.jl")
include("example/readme.jl")
include("example/empty.jl")
include("example/multi.jl")
end # module

@testset "err" begin
    @test_throws Jieko.NotImplementedError Jieko.not_implemented_error()
    @test contains(sprint(showerror, Jieko.NotImplementedError()), "Not implemented yet")
end

@testset "names" begin
    @static if VERSION > v"1.11-"
        # interface is public
        @test names(Example.Basic) == [Symbol("@goo"), Symbol("@moo"), :Basic, :Foo, :Prelude, :X, :foo]
    else
        @test names(Example.Basic) == [:Basic]
    end

    @test names(Example.Basic.Prelude) == [Symbol("@goo"), Symbol("@moo"), :Basic, :Foo, :Prelude, :X, :foo]
    @test names(Example.EmptyPrelude.Prelude) == [:Prelude]
    @test names(Example.Empty) == [:Empty]
    stub = Jieko.stub(Example.Empty)
    @test isempty(stub.macros)
    @test isempty(stub.interface)
    stub = Jieko.stub(Example.Basic)
    @test !isempty(stub.macros)
    @test !isempty(stub.interface)
end

@testset "doc" begin
    @test contains(sprint(show, @doc(Example.Basic.@goo)), "@goo <x::Int> <y::String> <zs>...")
    @test contains(sprint(show, @doc(Example.Basic.@moo)), "@moo <x::Int> [<y::String> = \"aaa\"]")
    @test contains(sprint(show, @doc(Example.Basic.Foo)), "struct Foo <: Real")
    @test contains(sprint(show, @doc(Example.Basic.foo)), "foo(x::Float64) -> Int")
    @test contains(sprint(show, @doc(Example.Basic.X)), "X")

    moduledoc = sprint(show, @doc(Example.Basic))
    @test contains(moduledoc, "### Prelude")
    @test contains(moduledoc, "### Definitions")
    @test contains(moduledoc, "#### Constants")
    @test contains(moduledoc, "#### Macros")
    @test contains(moduledoc, "#### Structs")
    @test contains(moduledoc, "#### Interfaces")
    @test contains(moduledoc, """```julia
    X
    ```
    """)
    @test contains(moduledoc, """```julia
    @goo <x::Int> <y::String> <zs>...
    ```""")
    @test contains(moduledoc, """```julia
    @moo <x::Int> [<y::String> = "aaa"]
    ```""")
    @test contains(moduledoc, """```julia
    struct Foo <: Real
    ```""")
    @test contains(moduledoc, """```julia
    foo(x::Float64) -> Int
    ```""")
    @test contains(moduledoc, """```julia
    foo(x::Float64) -> Int
    ```""")

    docstr = sprint(show, @doc(Example.Multi.jieko))
    @test contains(docstr, "jieko(x::Real) -> Int")
    @test contains(docstr, "jieko(x::MyAliasName) -> Complex")
end # @testset "doc"

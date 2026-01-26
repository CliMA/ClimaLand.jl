@noarg mutable struct NoArg
end

@noarg mutable struct NoArg2
    no_type
    with_type::Int
    with_default = 1
    with_type_default::Int = 1
    with_tag &(xml=(key="with-tag",),)
    with_tag_type::Int &(xml=(key="with-tag-type",),)
    with_tag_default = 1 &(xml=(key="with-tag-default",),)
    with_tag_type_default::Int = 1 &(xml=(key="with-tag-default",),)
    @atomic no_type_atomic
    @atomic with_type_atomic::Int
    # atomic fields w/ defaults *must* use parens to disambiguate from atomic
    # variable assignment syntax: https://github.com/JuliaLang/julia/issues/53893
    @atomic(with_default_atomic) = 1
    @atomic(with_type_default_atomic::Int) = 1
    @atomic(with_tag_atomic) &(xml=(key="with-tag-atomic",),)
    @atomic(with_tag_type_atomic::Int) &(xml=(key="with-tag-type-atomic",),)
    @atomic(with_tag_default_atomic) = 1 &(xml=(key="with-tag-default-atomic",),)
    @atomic(with_tag_type_default_atomic::Int) = 1 &(xml=(key="with-tag-default-atomic",),)
end

abstract type AbstractNoArg end

@noarg mutable struct NoArg3{T, S <: IO} <: AbstractNoArg
    a::T = 10 * 20
    io::S
end

"""
Documentation for NoArg4
"""
@noarg mutable struct NoArg4
    a::String
end

@kwarg mutable struct KwDef1
    no_type
    with_type::Int
    with_default = 1
    with_type_default::Int = 1
    with_tag &(xml=(key="with-tag",),)
    with_tag_type::Int &(xml=(key="with-tag-type",),)
    with_tag_default = 1 &(xml=(key="with-tag-default",),)
    with_tag_type_default::Int = 1 &(xml=(key="with-tag-default",),)
    @atomic no_type_atomic
    @atomic with_type_atomic::Int
    @atomic(with_default_atomic) = 1
    @atomic(with_type_default_atomic::Int) = 1
    @atomic(with_tag_atomic) &(xml=(key="with-tag-atomic",),)
    @atomic(with_tag_type_atomic::Int) &(xml=(key="with-tag-type-atomic",),)
    @atomic(with_tag_default_atomic) = 1 &(xml=(key="with-tag-default-atomic",),)
    @atomic(with_tag_type_default_atomic::Int) = 1 &(xml=(key="with-tag-default-atomic",),)
    const no_type_const
    const with_type_const::Int
    const with_default_const = 1
    const with_type_default_const::Int = 1
    const with_tag_const &(xml=(key="with-tag-const",),)
    const with_tag_type_const::Int &(xml=(key="with-tag-type-const",),)
    const with_tag_default_const = 1 &(xml=(key="with-tag-default-const",),)
    const with_tag_type_default_const::Int = 1 &(xml=(key="with-tag-default-const",),)
end

@kwarg struct NoFields
end

@kwarg struct KwDef2{T, S <: IO} <: AbstractNoArg
    a::T = 10 * 20
    io::S
end

# tests from Base.@kwdef misc.jl
@kwarg struct Test27970Typed
    a::Int
    b::String = "hi"
end

@kwarg struct Test27970Untyped
    a
end

@kwarg struct Test27970Empty end

abstract type AbstractTest29307 end
@kwarg struct Test29307{T<:Integer} <: AbstractTest29307
    a::T=2
end

@kwarg struct TestInnerConstructor
    a = 1
    TestInnerConstructor(a::Int) = (@assert a>0; new(a))
    function TestInnerConstructor(a::String)
        @assert length(a) > 0
        new(a)
    end
end

const outsidevar = 7
@kwarg struct TestOutsideVar
    a::Int=outsidevar
end
@test TestOutsideVar() == TestOutsideVar(7)

@kwarg mutable struct Test_kwdef_const_atomic
    a
    b::Int
    c::Int = 1
    const d
    const e::Int
    const f = 1
    const g::Int = 1
    @atomic h::Int
end

@kwarg struct Test_kwdef_lineinfo
    a::String
end
# @testset "StructUtils.@kwdef constructor line info" begin
#     for method in methods(Test_kwdef_lineinfo)
#         @test method.file === Symbol(@__FILE__)
#         @test ((@__LINE__)-6) ≤ method.line ≤ ((@__LINE__)-5)
#     end
# end
@kwarg struct Test_kwdef_lineinfo_sparam{S<:AbstractString}
    a::S
end
# @testset "@kwdef constructor line info with static parameter" begin
#     for method in methods(Test_kwdef_lineinfo_sparam)
#         @test method.file === Symbol(@__FILE__)
#         @test ((@__LINE__)-6) ≤ method.line ≤ ((@__LINE__)-5)
#     end
# end

module KwdefWithEsc
    using StructUtils
    const Int1 = Int
    const val1 = 42
    macro define_struct()
        quote
            @kwarg struct $(esc(:Struct))
                a
                b = val1
                c::Int1
                d::Int1 = val1

                $(esc(quote
                    e
                    f = val2
                    g::Int2
                    h::Int2 = val2
                end))

                $(esc(:(i = val2)))
                $(esc(:(j::Int2)))
                $(esc(:(k::Int2 = val2)))

                l::$(esc(:Int2))
                m::$(esc(:Int2)) = val1

                n = $(esc(:val2))
                o::Int1 = $(esc(:val2))

                $(esc(:p))
                $(esc(:q)) = val1
                $(esc(:s))::Int1
                $(esc(:t))::Int1 = val1
            end
        end
    end
end

module KwdefWithEsc_TestModule
    using ..KwdefWithEsc
    const Int2 = Int
    const val2 = 42
    KwdefWithEsc.@define_struct()
end

@defaults struct Defaults1
    a::Int
    b::String
    c::Float64 = 3.14
    d::Int = 42
end

@defaults struct Defaults2{S <: IO, T} <: AbstractNoArg
    io::S
    a::T = 10 * 20
end

@defaults struct DefaultsAll
    id::Int = 0
    name::String = ""
    value::Float64 = 1.0
end

@defaults struct DefaultsPartial
    required::String
    optional1::Int = 10
    optional2::Bool = true
end

@tags struct TagsWithDefaults
    id::Int = 0 &(json=(name="identifier",),)
    name::String = "default" &(json=(name="full_name",),)
    score::Float64 = 0.0
end

@tags struct TagsPartialDefaults
    required::String &(json=(required=true,),)
    optional::Int = 42 &(json=(name="opt",),)
end

@testset "macros" begin

    @test NoArg() isa NoArg

    x = NoArg2()
    x.no_type = :hey
    @test x isa NoArg2
    @test StructUtils.noarg(StructUtils.DefaultStyle(), NoArg2)
    fd = StructUtils.fielddefaults(StructUtils.DefaultStyle(), NoArg2)
    @test fd.with_default == 1 && fd.with_type_default == 1 && fd.with_tag_default == 1 && fd.with_tag_type_default == 1 && fd.with_default_atomic == 1 && fd.with_type_default_atomic == 1 && fd.with_tag_default_atomic == 1 && fd.with_tag_type_default_atomic == 1
    ft = StructUtils.fieldtags(StructUtils.DefaultStyle(), NoArg2)
    @test ft.with_tag == (xml=(key="with-tag",),) && ft.with_tag_type == (xml=(key="with-tag-type",),) && ft.with_tag_default == (xml=(key="with-tag-default",),) && ft.with_tag_type_default == (xml=(key="with-tag-default",),) && ft.with_tag_atomic == (xml=(key="with-tag-atomic",),) && ft.with_tag_type_atomic == (xml=(key="with-tag-type-atomic",),) && ft.with_tag_default_atomic == (xml=(key="with-tag-default-atomic",),) && ft.with_tag_type_default_atomic == (xml=(key="with-tag-default-atomic",),)

    @static if VERSION >= v"1.11-DEV"
        @test @doc(NoArg4).text[1] == "Documentation for NoArg4\n"
    else
        @test string(@doc(NoArg4)) == "Documentation for NoArg4\n"
    end

    @test_throws ArgumentError @macroexpand StructUtils.@noarg struct NonMutableNoArg
        with_type::Int = 1
    end

    # can't use @noarg w/ struct w/ const fields
    @test_throws ArgumentError @macroexpand StructUtils.@noarg mutable struct NoArgConst
        const a::Int = 1
    end

    x = KwDef1(no_type=1, with_type=1, with_tag=1, with_tag_type=1, no_type_atomic=1, with_type_atomic=1, with_tag_atomic=1, with_tag_type_atomic=1, no_type_const=1, with_type_const=1, with_tag_const=1, with_tag_type_const=1)
    @test x.no_type == 1
    @test StructUtils.kwarg(StructUtils.DefaultStyle(), KwDef1)
    fd = StructUtils.fielddefaults(StructUtils.DefaultStyle(), KwDef1)
    @test fd.with_default == 1 && fd.with_type_default == 1 && fd.with_tag_default == 1 && fd.with_tag_type_default == 1 && fd.with_default_atomic == 1 && fd.with_type_default_atomic == 1 && fd.with_tag_default_atomic == 1 && fd.with_tag_type_default_atomic == 1 && fd.with_default_const == 1 && fd.with_type_default_const == 1 && fd.with_tag_default_const == 1 && fd.with_tag_type_default_const == 1
    ft = StructUtils.fieldtags(StructUtils.DefaultStyle(), KwDef1)
    @test ft.with_tag == (xml=(key="with-tag",),) && ft.with_tag_type == (xml=(key="with-tag-type",),) && ft.with_tag_default == (xml=(key="with-tag-default",),) && ft.with_tag_type_default == (xml=(key="with-tag-default",),) && ft.with_tag_atomic == (xml=(key="with-tag-atomic",),) && ft.with_tag_type_atomic == (xml=(key="with-tag-type-atomic",),) && ft.with_tag_default_atomic == (xml=(key="with-tag-default-atomic",),) && ft.with_tag_type_default_atomic == (xml=(key="with-tag-default-atomic",),) && ft.with_tag_const == (xml=(key="with-tag-const",),) && ft.with_tag_type_const == (xml=(key="with-tag-type-const",),) && ft.with_tag_default_const == (xml=(key="with-tag-default-const",),) && ft.with_tag_type_default_const == (xml=(key="with-tag-default-const",),)

    @test NoFields() isa NoFields

    @testset "No default values in @kwdef" begin
        @test Test27970Typed(a=1) == Test27970Typed(1, "hi")
        # Implicit type conversion (no assertion on kwarg)
        @test Test27970Typed(a=0x03) == Test27970Typed(3, "hi")
        @test_throws UndefKeywordError Test27970Typed()
    
        @test Test27970Untyped(a=1) == Test27970Untyped(1)
        @test_throws UndefKeywordError Test27970Untyped()
    
        # Just checking that this doesn't stack overflow on construction
        @test Test27970Empty() == Test27970Empty()
    end

    @testset "subtyped @kwdef" begin
        @test Test29307() == Test29307{Int}(2)
        @test Test29307(a=0x03) == Test29307{UInt8}(0x03)
        @test Test29307{UInt32}() == Test29307{UInt32}(2)
        @test Test29307{UInt32}(a=0x03) == Test29307{UInt32}(0x03)
    end

    @testset "@kwdef inner constructor" begin
        @test TestInnerConstructor() == TestInnerConstructor(1)
        @test TestInnerConstructor(a=2) == TestInnerConstructor(2)
        @test_throws AssertionError TestInnerConstructor(a=0)
        @test TestInnerConstructor(a="2") == TestInnerConstructor("2")
        @test_throws AssertionError TestInnerConstructor(a="")
    end

    @testset "const and @atomic fields in @kwdef" begin
        x = Test_kwdef_const_atomic(a = 1, b = 1, d = 1, e = 1, h = 1)
        for f in fieldnames(Test_kwdef_const_atomic)
            @test getfield(x, f) == 1
        end
        @testset "const fields" begin
            @test_throws ErrorException x.d = 2
            @test_throws ErrorException x.e = 2
            @test_throws MethodError x.e = "2"
            @test_throws ErrorException x.f = 2
            @test_throws ErrorException x.g = 2
        end
        @testset "atomic fields" begin
            @test_throws ConcurrencyViolationError x.h = 1
            @atomic x.h = 1
            @test @atomic(x.h) == 1
            @atomic x.h = 2
            @test @atomic(x.h) == 2
        end
    end

    @test isdefined(KwdefWithEsc_TestModule, :Struct)

    @test Defaults1(1, "hey") == Defaults1(1, "hey", 3.14, 42)
    io = IOBuffer()
    @test Defaults2(io) == Defaults2(io, 200)

    @test_throws ArgumentError @macroexpand @defaults struct FooFoo
        a::Int = 1
        b::String
    end

    # Test partial constructors with defaults
    @testset "@defaults partial constructors" begin
        # Test DefaultsAll - all fields have defaults
        @test DefaultsAll() == DefaultsAll(0, "", 1.0)
        @test DefaultsAll(5) == DefaultsAll(5, "", 1.0)
        @test DefaultsAll(5, "test") == DefaultsAll(5, "test", 1.0)
        @test DefaultsAll(5, "test", 2.5) == DefaultsAll(5, "test", 2.5)
        
        # Test DefaultsPartial - mixed required and default fields
        @test DefaultsPartial("required") == DefaultsPartial("required", 10, true)
        @test DefaultsPartial("required", 20) == DefaultsPartial("required", 20, true)
        @test DefaultsPartial("required", 20, false) == DefaultsPartial("required", 20, false)
    end

    @testset "@tags with defaults partial constructors" begin
        # Test TagsWithDefaults - all fields have defaults
        @test TagsWithDefaults() == TagsWithDefaults(0, "default", 0.0)
        @test TagsWithDefaults(1) == TagsWithDefaults(1, "default", 0.0)
        @test TagsWithDefaults(1, "custom") == TagsWithDefaults(1, "custom", 0.0)
        @test TagsWithDefaults(1, "custom", 99.9) == TagsWithDefaults(1, "custom", 99.9)

        # Test TagsPartialDefaults - mixed required and default fields
        @test TagsPartialDefaults("needed") == TagsPartialDefaults("needed", 42)
        @test TagsPartialDefaults("needed", 100) == TagsPartialDefaults("needed", 100)

        # Test that field tags are preserved
        ft = StructUtils.fieldtags(StructUtils.DefaultStyle(), TagsWithDefaults)
        @test ft.id == (json=(name="identifier",),)
        @test ft.name == (json=(name="full_name",),)

        # Test that field defaults are preserved
        fd = StructUtils.fielddefaults(StructUtils.DefaultStyle(), TagsWithDefaults)
        @test fd.id == 0
        @test fd.name == "default"
        @test fd.score == 0.0
    end

end # @testset "macros"
# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using Test
import MutableArithmetics as MA

@testset "Int" begin
    a = [1, 2]
    b = 3
    c = [3, 4]
    @test MA.broadcast!!(+, a, b) == [4, 5]
    @test a == [4, 5]
    # Need to have immutable structs that contain references to be allocated on
    # the stack: https://github.com/JuliaLang/julia/pull/33886
    alloc_test(() -> MA.broadcast!!(+, a, b), 0)
    alloc_test(() -> MA.broadcast!!(+, a, c), 0)
end

@testset "BigInt" begin
    x = BigInt(1)
    y = BigInt(2)
    a = [x, y]
    b = 3
    c = [2x, 3y]
    @test MA.broadcast!!(+, a, b) == [4, 5]
    @test a == [4, 5]
    @test x == 4
    @test y == 5
    # FIXME This should not allocate but I couldn't figure out where these
    #       allocations come from.
    n = (VERSION >= v"1.11" ? 42 : 30) * sizeof(Int)
    alloc_test(() -> MA.broadcast!!(+, a, b), n)
    alloc_test(() -> MA.broadcast!!(+, a, c), 0)
end

@testset "broadcast_issue_158" begin
    x, y = BigInt[2 3], BigInt[2 3; 3 4]
    @test MA.@rewrite(x .+ y) == x .+ y
    @test MA.@rewrite(x .- y) == x .- y
    @test MA.@rewrite(y .+ x) == y .+ x
    @test MA.@rewrite(y .- x) == y .- x
    @test MA.@rewrite(y .* x) == y .* x
    @test MA.@rewrite(x .* y) == x .* y
end

struct Struct221 <: AbstractArray{Int,1} end
struct BroadcastStyle221 <: Broadcast.BroadcastStyle end
Base.BroadcastStyle(::Type{Struct221}) = BroadcastStyle221()

@testset "promote_broadcast_for_new_style" begin
    @test MA.promote_broadcast(MA.add_mul, Vector{Int}, Struct221) === Any
end

@testset "broadcast_length_1_dimensions" begin
    A = rand(2, 1, 3)
    B = rand(2, 3)
    @test MA.broadcast!!(MA.sub_mul, A, B) ≈ A .- B
    @test MA.broadcast!!(MA.sub_mul, B, A) ≈ B .- A
end

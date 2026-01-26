# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using Test
import MutableArithmetics as MA

include("utilities.jl")

include("interface.jl")

include("range.jl")

include("copy.jl")

@testset "Int" begin
    include("int.jl")
end
@testset "BigInt" begin
    include("big.jl")
end
@testset "BigFloat negation and absolute value" begin
    include("bigfloat_neg_abs.jl")
end
@testset "BigFloat fma and muladd" begin
    include("bigfloat_fma.jl")
end
@testset "BigFloat dot" begin
    include("bigfloat_dot.jl")
end
@testset "evalpoly" begin
    include("evalpoly.jl")
end
@testset "Broadcast" begin
    include("broadcast.jl")
end
include("matmul.jl")
include("dispatch.jl")
include("rewrite.jl")
include("rewrite_generic.jl")

include("SparseArrays.jl")

# It is easy to introduce macro scoping issues into MutableArithmetics,
# particularly ones that rely on `MA` or `MutableArithmetics` being present in
# the current scope. To work around that, include the "hygiene" script in a
# clean module with no other scope.

function _include_sandbox(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

_include_sandbox("hygiene.jl")

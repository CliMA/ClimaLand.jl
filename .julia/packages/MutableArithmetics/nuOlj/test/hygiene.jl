# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# hygiene.jl
# Make sure that our macros have good hygiene

module M
using Test
# Import it as a different name so that we test whether MutableArithmetics is
# needed in the current scope.
import MutableArithmetics as NewSymbolMA

macro _rewrite(expr)
    variable, code = NewSymbolMA.rewrite(expr)
    return quote
        $code
        $variable
    end
end

# Don't include this in `@test` to make sure it is in global scope
x = NewSymbolMA.@rewrite sum(i for i in 1:10)
@test x == 55
x = @_rewrite sum(i for i in 1:10)
@test x == 55

x = NewSymbolMA.@rewrite sum(i for i in 1:10 if isodd(i))
@test x == 25
x = @_rewrite sum(i for i in 1:10 if isodd(i))
@test x == 25

x = NewSymbolMA.@rewrite sum(i * j for i in 1:4 for j in 1:4 if i == j)
@test x == 30
x = @_rewrite sum(i * j for i in 1:4 for j in 1:4 if i == j)
@test x == 30

x = big(1)
y = NewSymbolMA.@rewrite(x + (x + 1)^1)
@test y == 3

end

# Test the scoping outside the module. See also the note in runtests.jl.

using Test

@testset "test_scoping" begin
    @test M.@_rewrite(1 + (1 + 1)^1) == 3
end

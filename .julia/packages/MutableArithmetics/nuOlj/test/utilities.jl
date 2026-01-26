# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

include("dummy.jl")

# Allocating size for allocating a `BigInt`. Half size on 32-bit.
const BIGINT_ALLOC = @static if VERSION >= v"1.12-beta1"
    Sys.WORD_SIZE == 64 ? 72 : 36
elseif VERSION >= v"1.11"
    Sys.WORD_SIZE == 64 ? 56 : 28
else
    Sys.WORD_SIZE == 64 ? 48 : 24
end

function alloc_test(f, n)
    f() # compile
    @test n == @allocated f()
end

function alloc_test_le(f, n)
    f() # compile
    @test n >= @allocated f()
end

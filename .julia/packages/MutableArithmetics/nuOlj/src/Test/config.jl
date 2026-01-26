# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

macro test_suite(setname, subsets = false)
    testname = Symbol(string(setname) * "_test")
    testdict = Symbol(string(testname) * "s")
    if subsets
        runtest = :(f(args...; exclude = exclude))
    else
        runtest = :(f(args...))
    end
    return esc(
        :(function $testname(args...; exclude::Vector{String} = String[])
            for (name, f) in $testdict
                if name in exclude
                    continue
                end
                @testset "$name" begin
                    $runtest
                end
            end
        end),
    )
end

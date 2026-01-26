# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using Documenter, MutableArithmetics

makedocs(;
    modules = [MutableArithmetics],
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = String[],
    ),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = ["Home" => "index.md"],
    repo = "https://github.com/jump-dev/MutableArithmetics.jl/blob/{commit}{path}#L{line}",
    sitename = "MutableArithmetics",
    authors = "Gilles Peiffer, Benoît Legat, and Sascha Timme",
)

deploydocs(; repo = "github.com/jump-dev/MutableArithmetics.jl.git")

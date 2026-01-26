# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

module Test

import MutableArithmetics as MA

using LinearAlgebra, SparseArrays, Test

include("config.jl")

include("int.jl")
include("generic.jl")
include("scalar.jl")
include("quadratic.jl")
include("array.jl")
include("sparse.jl")

end # module

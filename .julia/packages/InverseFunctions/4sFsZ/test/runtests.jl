# This file is a part of InverseFunctions.jl, licensed under the MIT License (MIT).

import Test
import InverseFunctions
import Documenter

Test.@testset "Package InverseFunctions" begin
    include("test_functions.jl")
    include("test_inverse.jl")
    include("test_setinverse.jl")

    # doctests
    Documenter.DocMeta.setdocmeta!(
        InverseFunctions,
        :DocTestSetup,
        :(using InverseFunctions);
        recursive=true,
    )
    Documenter.doctest(InverseFunctions)
end # testset

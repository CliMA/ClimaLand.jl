using Test
using Documenter
import ClimaDiagnostics

@testset "Test docstrings" begin

    DocMeta.setdocmeta!(
        ClimaDiagnostics,
        :DocTestSetup,
        :(using Dates);
        recursive = true,
    )

    doctest(ClimaDiagnostics; manual = false)
end

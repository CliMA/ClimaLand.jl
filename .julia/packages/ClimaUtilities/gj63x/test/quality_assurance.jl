using Test
using ClimaUtilities
using Aqua

using Documenter

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    ua = Aqua.detect_unbound_args_recursively(ClimaUtilities)
    @test length(ua) == 0

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(ClimaUtilities; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("ClimaUtilities", pkgdir(last(x).module)), ambs)

    # Uncomment for debugging:
    # for method_ambiguity in ambs
    #     @show method_ambiguity
    # end
    @test length(ambs) == 0
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(ClimaUtilities)
    Aqua.test_stale_deps(ClimaUtilities)
    Aqua.test_deps_compat(ClimaUtilities)
    Aqua.test_project_extras(ClimaUtilities)
    Aqua.test_piracies(ClimaUtilities)
end

@testset "Test docstrings" begin

    DocMeta.setdocmeta!(
        ClimaUtilities,
        :DocTestSetup,
        :(using Dates);
        recursive = true,
    )

    doctest(ClimaUtilities; manual = false)
end

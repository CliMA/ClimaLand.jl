using Test
using ClimaLand
using Aqua

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    ua = Aqua.detect_unbound_args_recursively(ClimaLand)
    @test length(ua) == 0

    # See: https://github.com/SciML/SciMLBase.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(ClimaLand; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("ClimaLand", pkgdir(last(x).module)), ambs)

    # Uncomment for debugging:
    # for method_ambiguity in ambs
    #     @show method_ambiguity
    # end
    @test length(ambs) == 0
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(ClimaLand)
    # Insolation is a direct dependency for the package extension
    Aqua.test_stale_deps(ClimaLand; ignore = [:Insolation])
    Aqua.test_deps_compat(ClimaLand)
    Aqua.test_project_extras(ClimaLand)
    Aqua.test_piracies(ClimaLand)
end

nothing

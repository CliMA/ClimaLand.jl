using Test
using ClimaDiagnostics
using Aqua
import ExplicitImports

@testset "Aqua tests" begin
    Aqua.test_undefined_exports(ClimaDiagnostics)
    Aqua.test_stale_deps(ClimaDiagnostics)
    Aqua.test_deps_compat(ClimaDiagnostics)
    Aqua.detect_ambiguities(ClimaDiagnostics; recursive = true)
    Aqua.test_piracies(ClimaDiagnostics)
end

@testset "Explicit Imports" begin
    @test isnothing(ExplicitImports.check_no_implicit_imports(ClimaDiagnostics))

    # We import some variables to bring them to the top level, so that they are
    # easier to use in user packages. These are technically stale imports, so we
    # have to ignore them

    ignore = (
        :DiagnosticVariable,
        :DivisorSchedule,
        :ScheduledDiagnostic,
        :average_pre_output_hook!,
    )
    @test isnothing(
        ExplicitImports.check_no_stale_explicit_imports(
            ClimaDiagnostics;
            ignore,
        ),
    )

    # ClimaCore idiosyncrasies
    ignore = (:HDF5, :topology, :vertical_topology)
    @test isnothing(
        ExplicitImports.check_all_qualified_accesses_via_owners(
            ClimaDiagnostics;
            ignore,
        ),
    )
    @test isnothing(
        ExplicitImports.check_no_self_qualified_accesses(ClimaDiagnostics),
    )
end

using SymbolicIndexingInterface, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(SymbolicIndexingInterface)
    Aqua.test_ambiguities(SymbolicIndexingInterface, recursive = false)
    Aqua.test_deps_compat(SymbolicIndexingInterface)
    Aqua.test_piracies(SymbolicIndexingInterface)
    Aqua.test_project_extras(SymbolicIndexingInterface)
    Aqua.test_stale_deps(SymbolicIndexingInterface)
    Aqua.test_unbound_args(SymbolicIndexingInterface)
    Aqua.test_undefined_exports(SymbolicIndexingInterface)
end

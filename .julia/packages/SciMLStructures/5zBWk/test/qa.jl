using SciMLStructures, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(SciMLStructures)
    Aqua.test_ambiguities(SciMLStructures, recursive = false)
    Aqua.test_deps_compat(SciMLStructures)
    Aqua.test_piracies(SciMLStructures)
    Aqua.test_project_extras(SciMLStructures)
    Aqua.test_stale_deps(SciMLStructures)
    Aqua.test_unbound_args(SciMLStructures)
    Aqua.test_undefined_exports(SciMLStructures)
end

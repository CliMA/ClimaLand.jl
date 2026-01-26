using Aqua
Aqua.test_all(SortingAlgorithms, deps_compat=false)
# It's okay to use the latest version of test deps:
Aqua.test_deps_compat(SortingAlgorithms, check_extras=false)

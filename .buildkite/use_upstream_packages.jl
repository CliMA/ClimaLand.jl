# This script updates the .buildkite and test environments to use other versions
# of upstream packages. This functionality is primarily for testing by upstream
# packages.
#
# Following the same convention as ClimaCoupler's nightly pipeline:
#     Package        -> use the main branch
#     Package@1.2.3  -> use a registered version
#     Package#rev    -> use a branch or commit
#
# For example, you can run
#
#     julia .buildkite/use_upstream_packages.jl ClimaCore
#
# to use the main branch of ClimaCore.

import Pkg

"""
    package_spec(entry)

Return the `Pkg.PackageSpec` for `entry`, which is one of
- `Package`: the main branch,
- `Package@1.2.3`: the registered version 1.2.3,
- `Package#rev`: the branch or commit `rev`.
"""
function package_spec(entry)
    if occursin("#", entry)
        name, rev = split(entry, "#"; limit = 2)
        return Pkg.PackageSpec(; name = String(name), rev = String(rev))
    elseif occursin("@", entry)
        name, version = split(entry, "@"; limit = 2)
        return Pkg.PackageSpec(; name = String(name), version = String(version))
    else
        return Pkg.PackageSpec(; name = entry, rev = "main")
    end
end

specs = map(package_spec, ARGS)
upstream = [spec.name for spec in specs]
isempty(upstream) && error("No upstream packages given")

root = dirname(@__DIR__)
buildkite_project = joinpath(root, ".buildkite")
test_project = joinpath(root, "test")

# Relax compat bounds for upstream packages, so versions outside them still
# resolve
for project in (root, buildkite_project, test_project)
    Pkg.activate(project)
    deps = Pkg.project().dependencies
    for name in filter(name -> haskey(deps, name), upstream)
        @info "Relaxing [compat] entry" project name
        Pkg.compat(name, ">= 0")
    end
end

# Update both environments to use the requested upstream versions
for project in (buildkite_project, test_project)
    Pkg.activate(project)
    Pkg.add(specs)
end

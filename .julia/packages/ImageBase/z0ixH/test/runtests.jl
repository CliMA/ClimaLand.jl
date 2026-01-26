using ImageBase, OffsetArrays, StackViews
using ImageFiltering
using Test, TestImages, Aqua, Documenter

using OffsetArrays: IdentityUnitRange
include("testutils.jl")

@testset "ImageBase.jl" begin
    @testset "Project meta quality checks" begin
        if Base.VERSION >= v"1.3"
            # Not checking compat section for test-only dependencies
            Aqua.test_ambiguities(ImageBase)
            Aqua.test_all(ImageBase;
                ambiguities=false,
                project_extras=true,
                deps_compat=true,
                stale_deps=true,
                project_toml_formatting=true
            )
            doctest(ImageBase, manual = false)
        end
    end

    include("diff.jl")
    include("restrict.jl")
    include("statistics.jl")

    if Base.JLOptions().depwarn != 2
        @info "deprecations are expected"
        include("deprecated.jl")
    end
end

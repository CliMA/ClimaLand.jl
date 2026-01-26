using SafeTestsets, Pkg

const GROUP = get(ENV, "GROUP", "All")

@time begin
    if GROUP == "All" || GROUP == "Core"
        include("setup.jl")
        include("known_values.jl")
        include("array_index.jl")
        include("axes.jl")
        include("broadcast.jl")
        include("dimensions.jl")
        include("indexing.jl")
        include("ranges.jl")
        include("size.jl")
        include("stridelayout.jl")
        include("misc.jl")
        @time @safetestset "StaticArrays" begin include("staticarrays.jl") end
        @time @safetestset "OffsetArrays" begin include("offsetarrays.jl") end
    end
end
using CommonDataModel
using CommonDataModel: MemoryDataset
using Test

@testset "CommonDataModel" begin
    #include("test_conversion.jl")
    include("test_empty.jl")
    include("test_scaling.jl")
    include("test_variable.jl")
    include("test_attrib.jl")
    include("test_copy.jl")
end

@testset "CF conventions" begin
    include("test_cfconventions.jl")
    include("test_coord.jl")
    include("test_bounds.jl")
end

@testset "Multi-file" begin
    include("test_multifile.jl")
end

@testset "views" begin
    include("test_subvariable.jl")
end

@testset "@select macro" begin
    include("test_select.jl")
    include("test_multifile_select.jl")
end

@testset "groupby" begin
    include("test_groupby.jl")
    include("test_rolling.jl")
end

@testset "aqua checks" begin
    include("test_aqua.jl")
end

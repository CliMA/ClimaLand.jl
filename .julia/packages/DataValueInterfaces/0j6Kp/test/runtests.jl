using Test, DataValueInterfaces

@testset "DataValueInterfaces" begin

@testset "nondatavaluetype" begin

    @test DataValueInterfaces.nondatavaluetype(Int64) == Int64
    @test DataValueInterfaces.nondatavaluetype(Union{}) == Union{}

end

@testset "datavaluetype" begin

    @test DataValueInterfaces.datavaluetype(Int64) == Int64
    @test DataValueInterfaces.datavaluetype(Union{}) == Union{}

end

end # @testset "DataValueInterfaces"

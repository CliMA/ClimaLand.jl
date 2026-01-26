using Test
using PkgVersion
using Pkg

const VERSION = PkgVersion.@Version
const UUID = PkgVersion.@Uuid
const AUTHOR = PkgVersion.@Author

@testset "macro own package Version - Uuid - Author" begin
    @test VERSION isa VersionNumber
    @test UUID isa Base.UUID
    @test AUTHOR isa String
end

@testset "function other package Version - Uuid - Author" begin
    @test PkgVersion.Version(Pkg) isa VersionNumber
    @test PkgVersion.Uuid(Test) isa Base.UUID
    @test PkgVersion.Author(Pkg) isa String
    @test PkgVersion.Version(Pkg) == PkgVersion.Version(Pkg.Types)
    @test PkgVersion.Uuid(Core.Compiler, 99) == Base.UUID(99)
end

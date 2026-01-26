using RelocatableFolders, Test, JSON3

module M

using RelocatableFolders

const BARE_FILE = @path "bare-file.jl"
const DIR = @path "path"
const FILE = @path joinpath("path", "file.jl")
const OTHER = @path joinpath(DIR, "subfolder/other.jl") # issue 8
const IGNORE_RE = @path "path" r".md$"
const IGNORE_RE_REL = @path "path" r"^subfolder"
const IGNORE_RES = @path "path" [r".md$"]
const IGNORE_FN = @path "path" path -> endswith(path, ".md")
const IGNORE_FN_REL = @path "path" path -> !startswith(path, "file")
const JSON_FILE = @path joinpath("path", "file.json")

end

@testset "RelocatableFolders" begin
    tests = function ()
        @test isfile(M.BARE_FILE)
        @test read(M.BARE_FILE, String) == "# bare file.jl"

        @test isfile(M.FILE)
        @test read(M.FILE, String) == "# file.jl"

        @test isfile(M.OTHER)
        @test read(M.OTHER, String) == "# other.jl"

        @test isdir(M.DIR)
        @test sort(readdir(M.DIR)) == ["file.jl", "file.json", "subfolder", "text.md"]
        @test readdir(joinpath(M.DIR, "subfolder")) == ["other.jl"]

        @test read(joinpath(M.DIR, "file.jl"), String) == "# file.jl"
        @test read(joinpath(M.DIR, "text.md"), String) == "text.md"
        @test read(joinpath(M.DIR, "subfolder", "other.jl"), String) == "# other.jl"

        @test length(M.DIR.files) == 4
        @test M.DIR.mod == M
        @test M.DIR.path == joinpath(@__DIR__, "path")

        @test length(M.IGNORE_RE.files) == 3
        @test !haskey(M.IGNORE_RE.files, joinpath("path", "text.md"))
        @test length(M.IGNORE_RE_REL.files) == 3
        @test !haskey(M.IGNORE_RE_REL.files, joinpath("path", "text.md"))
        @test !haskey(M.IGNORE_RE_REL.files, joinpath("path", "file.jl"))
        @test length(M.IGNORE_RES.files) == 3
        @test !haskey(M.IGNORE_RES.files, joinpath("path", "text.md"))
        @test length(M.IGNORE_FN.files) == 3
        @test !haskey(M.IGNORE_FN.files, joinpath("path", "text.md"))
        @test length(M.IGNORE_FN_REL.files) == 2
        @test !haskey(M.IGNORE_FN_REL.files, joinpath("path", "file.jl"))

        # JSON3.read has a different interface on versions that support Julia < 1.6.
        if VERSION >= v"1.6"
            @test JSON3.read(M.JSON_FILE).key == "value"
            @test JSON3.read(M.JSON_FILE, Dict{String,Any}) == Dict{String,Any}("key" => "value")
        end
    end
    from, to = joinpath.(Ref(@__DIR__), ("path", "moved"))
    file_from, file_to = joinpath.(Ref(@__DIR__), ("bare-file.jl", "moved-bare-file.jl"))
    try
        tests()
        @test String(M.DIR) == joinpath(@__DIR__, "path")
        @test String(M.FILE) == joinpath(@__DIR__, "path", "file.jl")

        # Remove the referenced folder `DIR`.
        @test isdir(from)
        @test !isdir(to)
        @test isfile(file_from)
        @test !isfile(file_to)
        mv(from, to)
        mv(file_from, file_to)
        @test isdir(to)
        @test !isdir(from)
        @test isfile(file_to)
        @test !isfile(file_from)

        let path = String(M.DIR)
            rm(path; recursive=true)
            @test !isdir(path)
        end
        let file = String(M.FILE)
            @test isfile(file)
            @test file != joinpath(@__FILE__, "path", "file.jl")
        end
        tests()
        @test String(M.DIR) != joinpath(@__DIR__, "path")
        @test String(M.FILE) != joinpath(@__DIR__, "path", "file.jl")

        # Test that `getroot` returns the correct values and does not generate files
        let path = String(M.DIR)
            rm(path, recursive=true)
            @test getroot(M.DIR) == joinpath(@__DIR__, "path")
            @test !isdir(getroot(M.DIR))
            @test getroot(M.DIR, path) == path
            @test !isdir(getroot(M.DIR, path))
        end
        let file = String(M.FILE)
            rm(file)
            root = dirname(file)
            @test getroot(M.FILE) == joinpath(@__DIR__, "path", "file.jl")
            @test !isfile(getroot(M.FILE))
            @test getroot(M.FILE, root) == file
            @test !isfile(getroot(M.FILE, root))
        end

        # Return the referenced folder `DIR`.
        @test !isdir(from)
        @test isdir(to)
        mv(to, from)
        @test isdir(from)
        @test !isdir(to)

        tests()
        @test String(M.DIR) == joinpath(@__DIR__, "path")
        @test String(M.FILE) == joinpath(@__DIR__, "path", "file.jl")

        @test_throws(ErrorException, @path("path", 123))
    finally
        isdir(from) || mv(to, from)
        isfile(file_from) || mv(file_to, file_from)
    end
end

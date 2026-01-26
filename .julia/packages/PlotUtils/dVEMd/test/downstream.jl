using Pkg, PlotUtils, Test

LibGit2 = Pkg.GitTools.LibGit2
TOML = Pkg.TOML

failsafe_clone_checkout(path, url) = begin
    local repo
    for i ∈ 1:6
        try
            repo = Pkg.GitTools.ensure_clone(stdout, path, url)
            break
        catch err
            @warn err
            sleep(20i)
        end
    end

    @assert isfile(joinpath(path, "Project.toml")) "spurious network error: clone failed, bailing out"

    name, _ = splitext(basename(url))
    registries = joinpath(first(DEPOT_PATH), "registries")
    general = joinpath(registries, "General")
    versions = joinpath(general, name[1:1], name, "Versions.toml")
    if !isfile(versions)
        mkpath(general)
        run(setenv(`tar xf $general.tar.gz`; dir = general))
    end
    @assert isfile(versions)

    stable = maximum(VersionNumber.(keys(TOML.parse(read(versions, String)))))
    tag = LibGit2.GitObject(repo, "v$stable")
    hash = string(LibGit2.target(tag))
    LibGit2.checkout!(repo, hash)
    nothing
end

fake_supported_version!(path) = begin
    toml = joinpath(path, "Project.toml")
    # fake the supported PlotUtils version for testing (for `Pkg.develop`)
    PlotUtils_version =
        Pkg.Types.read_package(normpath(@__DIR__, "..", "Project.toml")).version
    parsed_toml = TOML.parse(read(toml, String))
    parsed_toml["compat"]["PlotUtils"] = string(PlotUtils_version)
    open(toml, "w") do io
        TOML.print(io, parsed_toml)
    end
    nothing
end

develop_stable_Plots() = begin
    tmpd = mktempdir()
    Plots_jl = joinpath(tmpd, "Plots.jl")

    failsafe_clone_checkout(Plots_jl, "https://github.com/JuliaPlots/Plots.jl")
    fake_supported_version!(Plots_jl)

    Pkg.develop(path = Plots_jl)
    Pkg.status(["PlotUtils", "Plots"])
    nothing
end

develop_stable_Makie(extended = false) = begin
    tmpd = mktempdir()
    Makie_jl = joinpath(tmpd, "Makie.jl")

    failsafe_clone_checkout(Makie_jl, "https://github.com/MakieOrg/Makie.jl")
    fake_supported_version!(Makie_jl)

    Pkg.develop(path = joinpath(tmpd, "Makie.jl", "MakieCore"))
    Pkg.develop(path = joinpath(tmpd, "Makie.jl"))
    if extended  # too costly ?
        Pkg.develop(path = joinpath(tmpd, "Makie.jl", "ReferenceTests"))
        Pkg.develop(path = joinpath(tmpd, "Makie.jl", "CairoMakie"))
        # Pkg.develop(path = joinpath(tmpd, "Makie.jl", "GLMakie"))
    end
    Pkg.status(["PlotUtils", "MakieCore", "Makie"])
    nothing
end

develop_stable_Plots()
using Plots

@testset "downstream Plots" begin
    # test basic plots creation & display (Plots tests are too long to run)
    withenv("GKSwstype" => "nul") do
        @time for i ∈ 1:length(Plots._examples)
            i ∈ Plots._backend_skips[:gr] && continue  # skip unsupported examples
            Plots._examples[i].imports ≡ nothing || continue  # skip examples requiring optional test deps
            show(devnull, Plots.test_examples(:gr, i; disp = false))  # trigger display logic
        end
    end
end

extended = tryparse(Bool, get(ENV, "CI", "false")) === true  # extended test in CI

develop_stable_Makie(extended)
@testset "downstream Makie" begin
    Pkg.test("Makie")
    extended && Pkg.test("CairoMakie")
end

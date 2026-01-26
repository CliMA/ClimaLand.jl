using BaseDirs
using Documenter
using Org

orgfiles = filter(f -> endswith(f, ".org"),
                  readdir(joinpath(@__DIR__, "src"), join=true))

for orgfile in orgfiles
    mdfile = replace(orgfile, r"\.org$" => ".md")
    read(orgfile, String) |>
        c -> Org.parse(OrgDoc, c) |>
        o -> sprint(markdown, o) |>
        s -> replace(s, r"\.org]" => ".md]") |>
        m -> write(mdfile, m)
end

makedocs(;
    modules=[BaseDirs],
    authors="TEC <git@tecosaur.net> and contributors",
    repo="https://github.com/tecosaur/BaseDirs.jl/blob/{commit}{path}#{line}",
    sitename="BaseDirs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tecosaur.github.io/BaseDirs.jl",
        edit_link="main",
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "Introduction" => "index.md",
        "Usage" => "usage.md",
        "Defaults" => "defaults.md",
        "Prior Art" => "others.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(;
    repo="github.com/tecosaur/BaseDirs.jl",
    devbranch="main",
)

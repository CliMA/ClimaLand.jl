using Documenter, TableTraits

makedocs(
	modules = [TableTraits],
	sitename = "TableTraits.jl",
	analytics="UA-132838790-1",
	pages = [
        "Introduction" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/queryverse/TableTraits.jl.git"
)

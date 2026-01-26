using Documenter, IteratorInterfaceExtensions

makedocs(
	modules = [IteratorInterfaceExtensions],
	sitename = "IteratorInterfaceExtensions.jl",
	analytics="UA-132838790-1",
	pages = [
        "Introduction" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/queryverse/IteratorInterfaceExtensions.jl.git"
)

using Documenter

include(joinpath("..", "test", "test_and_analyze.jl"))

comparison_tables_file = joinpath("docs", "src", "comparison_tables.md")
preamble_file = joinpath("docs", "src", "comparison_tables_preamble.md")
cp(preamble_file, comparison_tables_file; force = true)
open(comparison_tables_file, "a") do io
    for (title, comparison_table_dict) in comparison_table_dicts
        print_comparison_table(title, comparison_table_dict, io)
    end
end

makedocs(;
    sitename = "UnrolledUtilities.jl",
    modules = [UnrolledUtilities],
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "User Guide" => "user_guide.md",
        "Developer Guide" => "developer_guide.md",
        "Comparison Tables" => basename(comparison_tables_file),
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = false,
        size_threshold_ignore = [
            "introduction.md",
            basename(comparison_tables_file),
        ],
    ),
    clean = true,
)

rm(comparison_tables_file)

deploydocs(
    repo = "github.com/CliMA/UnrolledUtilities.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
    forcepush = true,
)

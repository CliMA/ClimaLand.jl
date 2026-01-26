using Markdown
using Printf
using DelimitedFiles
using Statistics
using Glob

table = Any[["Module", "median", "minimum", "mean", "std. dev."]]

fmt(x) = @sprintf("%4.3f", x)
fmt(x) = @sprintf("%6.5f", x)

for f in glob("*.txt")
    data = readdlm(f, ' ', Float64)

    push!(
        table, [
            replace(f, ".txt" => ""),
            fmt(median(data)), fmt(minimum(data)),
            fmt(mean(data)), fmt(std(data)),
        ]
    )
end

print(Markdown.plain(Markdown.Table(table, [:l, :r, :r, :r, :r])))

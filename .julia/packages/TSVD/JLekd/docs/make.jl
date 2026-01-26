using Documenter, TSVD

DocMeta.setdocmeta!(
    TSVD,
    :DocTestSetup,
    :(using TSVD, MatrixDepot; matrixdepot("LPnetlib/lp_osa_30"));
    recursive=true)

makedocs(sitename="TSVD Documentation")

deploydocs(
    repo = "github.com/JuliaLinearAlgebra/TSVD.jl.git",
)

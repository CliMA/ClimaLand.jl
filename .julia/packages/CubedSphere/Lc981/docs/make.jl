using
  Documenter,
  DocumenterCitations,
  Literate,
  GLMakie,
  CubedSphere

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")
  
to_be_literated = [
    "compute_taylor_coefficients.jl"
]

for file in to_be_literated
    filepath = joinpath(EXAMPLES_DIR, file)
    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor(), execute = true)
    end
end


#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://clima.github.io/CubedSphere.jl/stable/",
)

pages = [
    "Home" => "index.md",
    "Conformal Cubed Sphere" => "conformal_cubed_sphere.md",
    "Advanced" => [
      "Reproduce Taylor coefficients for conformal mapping" => "literated/compute_taylor_coefficients.md",
      ],
    "References" => "references.md",
    "Library" => [ 
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],
]

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

makedocs(
   sitename = "CubedSphere.jl",
    modules = [CubedSphere],
    plugins = [bib],
     format = format,
      pages = pages,
    doctest = true,
      clean = true,
  checkdocs = :exports
)

deploydocs(        repo = "github.com/CliMA/CubedSphere.jl.git",
               versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
              forcepush = true,
              devbranch = "main",
           push_preview = true
)

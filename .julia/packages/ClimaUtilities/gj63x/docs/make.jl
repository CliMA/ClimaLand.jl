using Documenter
using ClimaUtilities

# Load everything to load extensions
import Interpolations
import ClimaComms
import ClimaCore
import NCDatasets
import ClimaCoreTempestRemap

DocMeta.setdocmeta!(
    ClimaUtilities,
    :DocTestSetup,
    :(using Dates);
    recursive = true,
)

pages = [
    "Overview" => "index.md",
    "Data Structures" => "datastructures.md",
    "ClimaArtifacts" => "climaartifacts.md",
    "Space and Time Inputs" => "inputs.md",
    "FileReaders" => "filereaders.md",
    "DataHandling" => "datahandling.md",
    "Regridders" => "regridders.md",
    "OnlineLogging" => "onlinelogging.md",
    "OutputPathGenerator" => "outputpathgenerator.md",
    "TimeManager" => "timemanager.md",
    "Frequently Asked Questions" => "faqs.md",
]

mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ),
)
format = Documenter.HTML(
    prettyurls = !isempty(get(ENV, "CI", "")),
    collapselevel = 1,
    mathengine = mathengine,
)

makedocs(
    sitename = "ClimaUtilities.jl",
    authors = "CliMA Utilities Developers",
    format = format,
    pages = pages,
    checkdocs = :exports,
    doctest = true,
    strict = false,
    clean = true,
    modules = [
        ClimaUtilities,
        Base.get_extension(ClimaUtilities, :ClimaUtilitiesClimaCoreExt),
        Base.get_extension(
            ClimaUtilities,
            :ClimaUtilitiesClimaCoreInterpolationsExt,
        ),
        Base.get_extension(
            ClimaUtilities,
            :ClimaUtilitiesClimaCoreNCDatasetsExt,
        ),
        Base.get_extension(ClimaUtilities, :ClimaUtilitiesNCDatasetsExt),
        Base.get_extension(
            ClimaUtilities,
            :ClimaUtilitiesClimaCoreTempestRemapExt,
        ),
    ],
)

deploydocs(
    repo = "github.com/CliMA/ClimaUtilities.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

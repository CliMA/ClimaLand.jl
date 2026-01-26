module UnicodeFun

using REPL

# Here, we collect and reverse the REPL's latex autocompletion map.
const symbols_unsorted = [k[2:end] => v[1] for (k, v) in REPL.REPLCompletions.latex_symbols]

const latex_symbol_map = sort!(symbols_unsorted, by=(x)-> length(x[1]), rev=true)

include("sub_super_scripts.jl")
export to_superscript, to_subscript
export to_fraction, to_fraction_nl

include("fontstyles.jl")
export to_blackboardbold
export to_boldface
export to_italic
export to_caligraphic
export to_frakture
export to_underline
export to_overline

include("latex.jl")
export to_latex

include("roots.jl")
export to_root

end # module

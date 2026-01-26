# UnicodeFun
[![Build Status](https://travis-ci.org/SimonDanisch/UnicodeFun.jl.svg?branch=master)](https://travis-ci.org/SimonDanisch/UnicodeFun.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ri3ybegh0ffwyq0n/branch/master?svg=true)](https://ci.appveyor.com/project/SimonDanisch/unicodefun-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/github/SimonDanisch/UnicodeFun.jl/badge.svg?branch=master)](https://coveralls.io/github/SimonDanisch/UnicodeFun.jl?branch=master)
[![codecov](https://codecov.io/gh/SimonDanisch/UnicodeFun.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SimonDanisch/UnicodeFun.jl)


unicode transformation library offering e.g. latex --> unicode:

```Julia
latexstring = "\\alpha^2 \\cdot \\alpha^{2+3} \\equiv \\alpha^7"
to_latex(latexstring) == "α² ⋅ α²⁺³ ≡ α⁷"
latexstring = "\\itA \\in \\bbR^{nxn}, \\bfv \\in \\bbR^n, \\lambda_i \\in \\bbR: \\itA\\bfv = \\lambda_i\\bfv"
==> "𝐴 ∈ ℝⁿˣⁿ, 𝐯 ∈ ℝⁿ, λᵢ ∈ ℝ: 𝐴𝐯 = λᵢ𝐯"
latexstring = "\\bf{boldface} \\it{italic} \\bb{blackboard} \\cal{calligraphic} \\frak{fraktur} \\mono{monospace}"
==> "𝐛𝐨𝐥𝐝𝐟𝐚𝐜𝐞 𝑖𝑡𝑎𝑙𝑖𝑐 𝕓𝕝𝕒𝕔𝕜𝕓𝕠𝕒𝕣𝕕 𝓬𝓪𝓵𝓵𝓲𝓰𝓻𝓪𝓹𝓱𝓲𝓬 𝔣𝔯𝔞𝔨𝔱𝔲𝔯 𝚖𝚘𝚗𝚘𝚜𝚙𝚊𝚌𝚎"
to_fraction_nl("α² ⋅ α²⁺³ ≡ α⁷", "ℝ: 𝐴𝐯 = λᵢ𝐯")
==>"
α̲²̲ ̲⋅̲ ̲α̲²̲⁺̲³̲ ̲≡̲ ̲α̲⁷̲
ℝ: 𝐴𝐯 = λᵢ𝐯"
```

Currently offered functions:
```Julia
to_superscript
to_subscript
to_fraction
to_fraction_nl
to_blackboardbold
to_boldface
to_italic
to_caligraphic
to_frakture
to_latex
```


Lookup tables taken from:
https://github.com/ypsu/latex-to-unicode/tree/master/data

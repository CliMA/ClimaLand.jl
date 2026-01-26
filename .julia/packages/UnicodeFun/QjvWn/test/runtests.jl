using UnicodeFun
using Test


# write your own tests here
@test UnicodeFun.to_superscript(-1234567890) == "⁻¹²³⁴⁵⁶⁷⁸⁹⁰"
@test UnicodeFun.to_subscript(-1234567890) == "₋₁₂₃₄₅₆₇₈₉₀"

latexstring = "\\alpha^2 \\cdot \\alpha^{2+3} \\equiv \\alpha^7"
@test to_latex(latexstring) == "α² ⋅ α²⁺³ ≡ α⁷"
latexstring = "\\itA \\in \\bbR^{nxn}, \\bfv \\in \\bbR^n, \\lambda_i \\in \\bbR: \\itA\\bfv = \\lambda_i\\bfv"
@test to_latex(latexstring) == "𝐴 ∈ ℝⁿˣⁿ, 𝐯 ∈ ℝⁿ, λᵢ ∈ ℝ: 𝐴𝐯 = λᵢ𝐯"
@test to_latex("\\sqrt\\cbrt") == "√∛"

latexstring = "\\bf{boldface} \\it{italic} \\bb{blackboard} \\cal{calligraphic} \\frak{fraktur} \\mono{monospace}"
@test to_latex(latexstring) == "𝐛𝐨𝐥𝐝𝐟𝐚𝐜𝐞 𝑖𝑡𝑎𝑙𝑖𝑐 𝕓𝕝𝕒𝕔𝕜𝕓𝕠𝕒𝕣𝕕 𝓬𝓪𝓵𝓵𝓲𝓰𝓻𝓪𝓹𝓱𝓲𝓬 𝔣𝔯𝔞𝔨𝔱𝔲𝔯 𝚖𝚘𝚗𝚘𝚜𝚙𝚊𝚌𝚎"

@test to_fraction("a-123", 392) == "ᵃ⁻¹²³⁄₃₉₂"

@test to_fraction_nl("α² ⋅ α²⁺³ ≡ α⁷", "ℝ: 𝐴𝐯 = λᵢ𝐯") == "α̲²̲ ̲⋅̲ ̲α̲²̲⁺̲³̲ ̲≡̲ ̲α̲⁷̲
ℝ: 𝐴𝐯 = λᵢ𝐯"
@test to_overline("abc") == "a̅b̅c̅"
@test to_underline("abc") == "a̲b̲c̲"

@test to_root("542") == "√5̅4̅2̅"
@test to_root(3,"542") == "∛5̅4̅2̅"
@test to_root(4,"542") == "∜5̅4̅2̅"
@test to_root(17,"542") == "¹⁷√5̅4̅2̅"
@test to_root(1, "1") == "1"
@test to_root(1, "-1") == "-1"

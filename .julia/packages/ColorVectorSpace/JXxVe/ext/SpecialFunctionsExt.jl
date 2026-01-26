module SpecialFunctionsExt
    import ColorVectorSpace
    using ColorTypes

    isdefined(Base, :get_extension) ? (import SpecialFunctions) : (import ..SpecialFunctions)

    const unaryops = (
        :gamma, :logfactorial, :erf, :erfc, :erfcx, :erfi, :dawson,
        :airyai, :airyaiprime, :airybi, :airybiprime,
        :besselj0, :besselj1, :bessely0, :bessely1,
        :eta, :zeta, :digamma
    )

    for op in unaryops
        @eval SpecialFunctions.$op(c::AbstractGray) = Gray(SpecialFunctions.$op(gray(c)))
    end

    function SpecialFunctions.logabsgamma(c::AbstractGray)
        lagc, s = SpecialFunctions.logabsgamma(gray(c))
        return Gray(lagc), s
    end
end

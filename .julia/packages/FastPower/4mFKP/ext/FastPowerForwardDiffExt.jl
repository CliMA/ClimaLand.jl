module FastPowerForwardDiffExt

import FastPower
using ForwardDiff

@inline FastPower.fastpower(x::ForwardDiff.Dual, y::ForwardDiff.Dual) = x^y

end

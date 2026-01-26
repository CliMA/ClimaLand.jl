module FastPowerReverseDiffExt

using FastPower, ReverseDiff
FastPower.fastpower(x::ReverseDiff.TrackedReal, y::ReverseDiff.TrackedReal) = x^y

end

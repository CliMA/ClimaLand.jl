module FastPowerTrackerExt

using FastPower, Tracker
FastPower.fastpower(x::Tracker.TrackedReal, y::Tracker.TrackedReal) = x^y

end

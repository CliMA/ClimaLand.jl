module FastPowerMeasurementsExt

using FastPower
using Measurements

@inline FastPower.fastpower(x::Measurements.Measurement, y::Measurements.Measurement) = x^y

end

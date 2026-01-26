module FastPowerMonteCarloMeasurementsExt

using FastPower
using MonteCarloMeasurements

@inline function FastPower.fastpower(x::MonteCarloMeasurements.AbstractParticles,
        y::MonteCarloMeasurements.AbstractParticles)
    x^y
end

end

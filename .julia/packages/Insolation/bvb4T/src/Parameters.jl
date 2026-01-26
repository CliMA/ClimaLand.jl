module Parameters

abstract type AbstractInsolationParams end
const AIP = AbstractInsolationParams

import Dates: DateTime

Base.@kwdef struct InsolationParameters{FT} <: AbstractInsolationParams
    year_anom::FT
    day::FT
    orbit_semimaj::FT
    tot_solar_irrad::FT
    epoch::DateTime
    mean_anom_epoch::FT
    eccentricity_epoch::FT
    obliq_epoch::FT
    lon_perihelion_epoch::FT
end

# Method wrappers
for var in fieldnames(InsolationParameters)
    @eval $var(ps::AIP) = ps.$var
end

end # module

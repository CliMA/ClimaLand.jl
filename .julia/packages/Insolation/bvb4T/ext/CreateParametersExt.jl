module CreateParametersExt

import Insolation.Parameters.InsolationParameters

import ClimaParams as CP

InsolationParameters(
    ::Type{FT},
    overrides = NamedTuple(),
) where {FT <: AbstractFloat} =
    InsolationParameters(CP.create_toml_dict(FT), overrides)

function InsolationParameters(
    toml_dict::CP.ParamDict{FT},
    overrides = NamedTuple(),
) where {FT}
    name_map = (;
        :epoch_time => :epoch,
        :day => :day,
        :anomalistic_year_length => :year_anom,
        :length_orbit_semi_major => :orbit_semimaj,
        :orbit_eccentricity_at_epoch => :eccentricity_epoch,
        :total_solar_irradiance => :tot_solar_irrad,
        :orbit_obliquity_at_epoch => :obliq_epoch,
        :mean_anomaly_at_epoch => :mean_anom_epoch,
        :longitude_perihelion_at_epoch => :lon_perihelion_epoch,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Insolation")
    parameters = merge(parameters, overrides)
    return InsolationParameters{FT}(; parameters...)
end

end

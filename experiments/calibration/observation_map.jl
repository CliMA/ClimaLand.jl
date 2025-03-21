using ClimaAnalysis
using Statistics
import EnsembleKalmanProcesses as EKP

"""
"""
function process_member_data(simdir, training_locations, variable_list; spinup_period = Year(1))

    # For later, for more flexible simulation times,
    # vardates = ClimaAnalysis.Utils.time_to_date.(Ref(start_date), ClimaAnalysis.times(lhf))
    # vardates[findlast(x->month(x)==11, vardates - Month(1))-1]

    start_date = DateTime(get(simdir, variable_list[1]).attributes["start_date"])
    left = start_date+spinup_period
    right = start_date+spinup_period+Year(1)
    data = Dict(var => window(shift_to_start_of_previous_month(get(simdir, var)), "time"; left, right) for var in variable_list)

    # We chose December as a starting date so that the default spinup_period (1 year) makes it easy to extract seasons (starting with DJF).
    @assert month(start_date) == 12 "The start date of the simulation should be December"

    obs = []

    # TODO: Function below inefficient, clean up later
    for (lon, lat) in training_locations
        for var in variable_list
            loc_data = ClimaAnalysis.slice(data[var]; lon, lat)
            seasonal_data = split_by_season_across_time(loc_data)
            seasonal = [mean(seasonal_data[i].data) for i in 1:length(seasonal_data)]
            push!(obs, seasonal)
        end
    end

    return vcat(obs...)
end

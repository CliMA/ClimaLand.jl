using Insolation
using Plots
using Dates
using Statistics
using Formatting

import ClimaParams as CP
import Insolation.Parameters as IP

"""
    calc_day_lat_insolation(od::Insolation.OrbitalData,
                            n_days::I,
                            n_lats::I,
                            param_set::IP.AIP) where {I<:Int}
"""
function calc_day_lat_insolation(
    od,
    n_days::I,
    n_lats::I,
    param_set::IP.AIP,
) where {I <: Int}
    d_arr = Array{I}(round.(collect(range(0, stop = 365, length = n_days))))
    l_arr = collect(range(-90, stop = 90, length = n_lats))
    F_arr = zeros(n_days, n_lats)
    # loop over days
    for (i, d) in enumerate(d_arr)
        for (j, lat) in enumerate(l_arr)
            date = Dates.DateTime(2000, 1, 1) + Dates.Day(d)
            θ, dist = daily_zenith_angle(
                date,
                od,
                lat,
                param_set,
                milankovitch = false,
            )
            F_arr[i, j] = insolation(θ, dist, param_set)
        end
    end
    return d_arr, l_arr, F_arr
end

"""
    plot_day_lat_insolation(n_days::I,
                            n_lats::I,
                            F_arr::Array{FT},
                            scmap,
                            stitle,
                            file_name) where {FT<:AbstractFloat,I<:Int}
"""
function plot_day_lat_insolation(
    d_arr::Array{I},
    l_arr::Array{FT},
    F_arr::Array{FT},
    scmap,
    stitle,
    file_name,
) where {FT <: AbstractFloat, I <: Int}
    if scmap == "YlOrRd"
        cmap = :YlOrRd
        vmin, vmax = 0, ceil(max(F_arr...) / 100) * 100
    elseif scmap == "PRGn"
        cmap = :PRGn
        vmin, vmax = ceil(max(abs.(F_arr)...) / 10) * -10,
        ceil(max(abs.(F_arr)...) / 10) * 10
    end

    p1 = contourf(
        d_arr,
        l_arr,
        F_arr',
        c = cmap,
        clims = (vmin, vmax),
        title = stitle,
        xlabel = "Days since Jan 1 2000",
        ylabel = "Latitude",
        colorbar_title = "ToA Insolation [W/m2]",
    )
    meanF = mean(F_arr, dims = 1)[1, :]
    x1, x2 = [floor(minimum(meanF)), ceil(maximum(meanF))]
    p2 = plot(
        meanF,
        l_arr,
        xlabel = "Annual-mean TOA \nInsolation [W/m2]",
        ylims = [-90, 90],
        yticks = [-90, -60, -30, 0, 30, 60, 90],
        xlims = [x1, x2],
        legend = false,
    )
    plot(
        p1,
        p2,
        size = (800, 400),
        layout = grid(1, 2, widths = (0.8, 0.2)),
        dpi = 250,
        left_margin = 4Plots.mm,
        bottom_margin = 6Plots.mm,
        right_margin = 6Plots.mm,
    )

    savefig(file_name)
end

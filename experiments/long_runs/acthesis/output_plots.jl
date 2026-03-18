my_param_dir = "global_all_models_f64"
old_param_dir = nothing
htessel_param_dir = nothing
era5_land_dir = nothing

base_dir = "/home/acharbon/thesis_outputs"
output_dir = "/home/acharbon/thesis_plots"

using NCDatasets, DataFrames, Statistics

function get_agg_ids(dates::Vector{DateTime}, period::Symbol, bin::Bool = false)
    if period == :year
        bin && error("Cannot bin yearly data")
        return year.(dates)
    elseif period == :season
        mon = month.(dates)
        seas = fill(:DJF, length(dates))
        seas[in.(mon, [[3, 4, 5]])] .= :MAM
        seas[in.(mon, [[6, 7, 8]])] .= :JJA
        seas[in.(mon, [[9, 10, 11]])] .= :SON
        bin && return seas
        yr = ifelse.(mon .== 12, year.(dates) + 1, year.(dates))
        return collect(zip(yr, seas))
    elseif period == :month
        return bin ? month.(dates) : collect(zip(year.(dates), month.(dates)))
    elseif period == :week
        wk = floor.(Int, dayofyear.(dates) ./ 7) .+ 1
        wk[wk .== 53] .= 52
        return bin ? wk : collect(zip(year.(dates), wk))
    elseif period == :day
        doy = dayofyear.(dates)
        doy[doy .== 366] .= 365
        return bin ? doy : nothing
    end
end

function agg_indices(data::Array, dates::Vector{DateTime}, period::Symbol, agg; bin::Bool = false)
    agg_ids = get_agg_ids(dates, period, bin)
    if isnothing(agg_ids)
        return data, dates
    else
        idset = unique(agg_ids)
        dat = Vector{Array}(undef, length(idset))
        for (i, id) in enumerate(idset)
            idxs = findall(==(id), idset)
            dat[i] = agg(data[idxs, :, :], dims = 1)
        end
        return cat(dat..., dims = 1), idset
    end
end

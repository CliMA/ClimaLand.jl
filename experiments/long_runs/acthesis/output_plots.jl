using NCDatasets, DataFrames, Statistics, GeoMakie
using ClimaLand

my_param_dir = "/home/acharbon/thesis_outputs/global_all_models_f64"
old_param_dir = nothing
htessel_param_dir = nothing
output_dir = "/home/acharbon/thesis_plots"
era5_land_obs = joinpath(ClimaLand.Artifacts.era5_surface_data_1979_2024_path(), "era5_monthly_averages_surface_single_level_197901-202410.nc")
#^^ this does not have everything, you need to get your own.

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

#mapslices(f, data, dims = 1)

function globe_heatmap(lat, lon, data)
    fig = Figure()
    viz_mask = GeoMakie.NaturalEarth.bathymetry(0).geometry
    ax = GeoMakie.GeoAxis(fig[1,1]; title = "plot title")
    plot = GeoMakie.surface!(ax, lon, lat, data; shading = NoShading)
    GeoMakie.poly!(ax, viz_mask;)
    GeoMakie.lines!(ax, GeoMakie.coastlines(); color = "black")
    GeoMakie.Colorbar(
            fig[1,2],
            plot,
            label = "variable name";
            height = Relative(0.65)
        )
    return fig
end
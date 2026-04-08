using NCDatasets, DataFrames, Statistics, CairoMakie, GeoMakie, Dates, JLD2, Statistics, HypothesisTests
using ClimaLand

function combine_clima_data(path)
    filenames = readdir(path)
    data = Dict()
    for file in filenames
        if file[end-2:end] == ".nc"
            varname = Symbol(split(file, "_")[1])
            data[varname] = NCDataset(joinpath(path, "$(string(varname))_1d_average.nc"))[varname]
        end
    end
    return data
end

const my_param_dir = "/home/acharbon/outputs/all_global"
const old_param_dir = "/home/acharbon/outputs/none_global"
const htessel_param_dir = "/home/acharbon/outputs/others_global"
const obs_param_dir = "/home/acharbon/outputs/era5_global"
const clm5_dir = "/home/acharbon/outputs/clm5_global"
const output_dir = "/home/acharbon/plots"

const me_outputs = combine_clima_data(my_param_dir)
const old_outputs = combine_clima_data(old_param_dir)
const others_outputs = combine_clima_data(htessel_param_dir)
const clm5_outputs = NCDataset(joinpath(clm5_dir, "all_clm_outputs.nc"))
const era5_rad_outputs = NCDataset(joinpath(obs_param_dir, "era5_monthly_rad.nc"))
const era5_snow_outputs = NCDataset(joinpath(obs_param_dir, "all_snow_vars.nc"))
era5_outputs(type) = type == :rad ? era5_rad_outputs : era5_snow_outputs

const elevs = NCDataset(joinpath(obs_param_dir, "elevation.nc"))[:elev][:,:]
const glacial_base = NCDataset(joinpath(obs_param_dir, "glacialfrac.nc"))[:glf][:,:]
const clima_land_mask = JLD2.load(joinpath(my_param_dir, "landmask.jld2"))["mask"]
const clima_lat = collect(range(-90, 89))
const clima_lon = collect(range(-180, 179))

function glacial_mask(thresh)
    mask = (glacial_base .<= thresh) .* clima_land_mask
    ret_mask = Matrix{Float64}(deepcopy(mask))
    ret_mask[ret_mask .== 1] .= 1.0
    ret_mask[ret_mask .< 1] .= NaN
    #add the four bad points:
    for (p1, p2) in [(111, 57), (261, 133), (316, 153), (302, 158)]
        ret_mask[p1, p2] = NaN
    end
    return ret_mask
end

function apply_mask(arr, mask::Matrix)
    return arr .* reshape(mask, 1, size(mask)...)
end

myargs(; type = :snow) = (clm5 = clm5_outputs, era5 = era5_outputs(type), me = me_outputs, others = others_outputs, old = old_outputs)
const snow_args = myargs()
const rad_args = myargs(type = :rad)

function get_data(field)
    dims = dimnames(field)
    lon_idx = findfirst(==("lon"), dims)
    lat_idx = findfirst(==("lat"), dims)
    time_idx = in("time", dims) ? findfirst(==("time"), dims) : findfirst(==("date"), dims)
    if any(isnothing.([lat_idx, lon_idx, time_idx]))
        error("No time or coordinate dimensions on provided field.")
    end
    dat = Array(permutedims(field, (time_idx, lon_idx, lat_idx)))
    return dat
end

function field_data(clima_name, era_name, clm5_name; args)
    return (
        me = get_data(args.me[clima_name]),
        old = get_data(args.old[clima_name]),
        others = get_data(args.others[clima_name]),
        era5 = get_data(args.era5[era_name]),
        clm5 = get_data(args.clm5[clm5_name]),
    ),
    (
        clima = args.me[clima_name][:date][:],
        era5 = args.era5[era_name][:date][:],
        clm5 = args.clm5[clm5_name][:date][:]
    )
end

function get_agg_ids(dates::Vector{<:Union{Date, DateTime}}, period::Symbol, bin::Bool = false)
    if period == :year
        bin && error("Cannot bin yearly data")
        return year.(dates)
    elseif period == :snowyear
        bin && error("Cannot bin yearly data")
        return year.(dates .+ Dates.Day(122))
    elseif period == :season
        mon = month.(dates)
        seas = fill(:DJF, length(dates))
        seas[in.(mon, Ref([3, 4, 5]))] .= :MAM
        seas[in.(mon, Ref([6, 7, 8]))] .= :JJA
        seas[in.(mon, Ref([9, 10, 11]))] .= :SON
        bin && return seas
        yr = ifelse.(mon .== 12, year.(dates) + 1, year.(dates))
        return collect(zip(yr, seas))
    elseif period == :month
        return bin ? month.(dates) : collect(zip(year.(dates), month.(dates)))
    elseif period == :week
        wk = ceil.(Int, dayofyear.(dates) ./ 7)
        wk = clamp.(wk, 1, 52)
        return bin ? wk : collect(zip(year.(dates), wk))
    elseif period == :day
        doy = dayofyear.(dates)
        doy[doy .== 366] .= 365
        return bin ? doy : nothing
    end
end

function agg_indices(data::Array, dates::Vector{<:Union{Date, DateTime}}, period::Symbol, agg; bin::Bool = false)
    agg_ids = get_agg_ids(dates, period, bin)
    if isnothing(agg_ids)
        return data, dates
    else
        idset = sort(unique(agg_ids))
        dat = Vector{Array}(undef, length(idset))
        for (i, id) in enumerate(idset)
            idxs = findall(==(id), agg_ids)
            dat[i] = agg(data[idxs, :, :], dims = 1)
        end
        return cat(dat..., dims = 1), idset
    end
end

is_data(x) = !ismissing(x) && !isnan(x)
valid_pairs(sim::Union{Vector, SubArray,}, obs::Union{Vector, SubArray}) = (
    (s, o) for (s, o) in zip(sim, obs)
    if is_data(s) && is_data(o)
)
function sumf(v::Union{Vector, SubArray})
    xs = [x for x in v if is_data(x)]
    isempty(xs) && return NaN
    return sum(xs)
end
function meanf(v::Union{Vector, SubArray}) #defined twice to make sure I am not collecting an empty array when not using over model outputs
    xs = [x for x in v if is_data(x)]
    isempty(xs) && return NaN
    return mean(xs)
end
function stdf(v::Union{Vector, SubArray}) #defined twice to make sure I am not collecting an empty array when not using over model outputs
    xs = [x for x in v if is_data(x)]
    isempty(xs) && return NaN
    return std(xs, corrected = false)
end
stdf(v) = std([x for x in v if is_data(x)], corrected = false)
maxf(v) = maximum([x for x in v if is_data(x)])
minf(v) = minimum([x for x in v if is_data(x)])
meanf(v) = mean([x for x in v if is_data(x)])
medianf(v) = median([x for x in v if is_data(x)])
quantilef(v, q) = quantile([x for x in v if is_data(x)], q)
mae(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(abs.(sim .- obs))
rmse(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = sqrt(meanf(abs2.(sim .- obs)))
bias(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(sim .- obs)
rmse_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = sqrt(meanf(abs2.(sim .- obs) ./ obs))
bias_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf((sim .- obs) ./ obs)
function nse(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray})
    p = collect(valid_pairs(sim, obs))
    s = first.(p)
    o = last.(p)
    μ = mean(o)
    return 1 - sum(abs2.(s .- o))/sum(abs2.(o .- μ))
end
function r2(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray})
    p = collect(valid_pairs(sim, obs))
    length(p) <= 2 && return NaN
    return cor(first.(p), last.(p))^2
end
function spatial_r2(sim::Union{Matrix, SubArray}, obs::Union{Matrix, SubArray})
    p = collect(valid_pairs(sim[:], obs[:]))
    length(p) <= 2 && return NaN
    return cor(first.(p), last.(p))^2
end

function integral_image(A)
    H, W = size(A)
    I = zeros(eltype(A), H+1, W+1)
    I[2:end, 2:end] .= cumsum(cumsum(A, dims=1), dims=2)
    return I
end
@inline function boxsum(I, x1, y1, x2, y2)
    return I[x2+1,y2+1] - I[x1,y2+1] - I[x2+1,y1] + I[x1,y1]
end
function ssim(sim, obs; r=5, C1=1e-4, C2=9e-4)
    H, W = size(sim)
    @assert size(obs) == (H, W)
    valid_data = is_data.(sim) .& is_data.(obs)
    simv = ifelse.(valid_data, sim, 0.0)
    obsv = ifelse.(valid_data, obs, 0.0)
    Ix   = integral_image(simv)
    Iy   = integral_image(obsv)
    Ix2  = integral_image(simv.^2)
    Iy2  = integral_image(obsv.^2)
    Ixy  = integral_image(simv .* obsv)
    In   = integral_image(Float64.(valid_data))
    total = 0.0
    count = 0
    for j in 1+r:H-r
        for i in 1+r:W-r
            x1, x2 = j-r, j+r
            y1, y2 = i-r, i+r

            n = boxsum(In, x1, y1, x2, y2)
            n < 2 && continue

            sx  = boxsum(Ix,  x1, y1, x2, y2)
            sy  = boxsum(Iy,  x1, y1, x2, y2)
            sx2 = boxsum(Ix2, x1, y1, x2, y2)
            sy2 = boxsum(Iy2, x1, y1, x2, y2)
            sxy = boxsum(Ixy, x1, y1, x2, y2)

            μx = sx / n
            μy = sy / n

            σx2 = sx2/n - μx^2
            σy2 = sy2/n - μy^2
            σxy = sxy/n - μx*μy

            num = (2*μx*μy + C1) * (2*σxy + C2)
            den = (μx^2 + μy^2 + C1) * (σx2 + σy2 + C2)

            if den != 0
                total += num / den
                count += 1
            end
        end
    end
    return count == 0 ? NaN : total / count
end
function time_stat(f, arr)
    T, LON, LAT = size(arr)
    out = fill(NaN, LON, LAT)
    for i in 1:LON
        for j in 1:LAT
            out[i, j] = f(view(arr, :, i, j))
        end
    end
    return out
end
function time_stat(f, arr1, arr2)
    T1, LON1, LAT1 = size(arr1)
    T2, LON2, LAT2 = size(arr2)
    @assert (LAT1 == LAT2) && (LON1 == LON2) && (T1 == T2)
    out = fill(NaN, LON1, LAT1)
    for i in 1:LON1
        for j in 1:LAT1
            out[i, j] = f(view(arr1, :, i, j), view(arr2, :, i, j))
        end
    end
    return out
end
function space_stat(f, arr)
    T, LON, LAT = size(arr)
    out = fill(NaN, T)
    for i in 1:T
        out[i] = f(view(arr, i, :, :))
    end
    return out
end
function space_stat(f, arr1, arr2)
    T1, LON1, LAT1 = size(arr1)
    T2, LON2, LAT2 = size(arr2)
    @assert (LAT1 == LAT2) && (LON1 == LON2) && (T1 == T2)
    out = fill(NaN, T1)
    for i in 1:T1
        out[i] = f(view(arr1, i, :, :), view(arr2, i, :, :))
    end
    return out
end

function stat_spread(v)
    return Dict(
        "min" => minf(v),
        "max" => maxf(v),
        "median" => medianf(v),
        "mean" => meanf(v),
        "q25" => quantilef(v, 0.25),
        "q75" => quantilef(v, 0.75),
        "std" => stdf(v)
    )
end

function make_nosnow_mask(; z_thresh = 0, glacial_thresh = 0)
    era5_snow_depth = get_data(era5_snow_outputs[:snd])
    max_snow_cell = maximum(era5_snow_depth, dims = 1)[1, :, :]
    max_snow_cell[ismissing.(max_snow_cell)] .= NaN
    no_snow = findall(<(z_thresh), max_snow_cell)
    mask = glacial_mask(glacial_thresh)
    mask[no_snow] .= NaN
    return mask
end

function depth_analysis(; args = snow_args)
    @info "Running Depth Analysis..."
    print("   Extracting data...\n")
    z, dates = field_data(:snd, :snd, :SNOW_DEPTH_month; args = args)
    calcs = [(time_stat, bias), (time_stat, rmse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    #make data mask: glacial mask, plus spots where there is no snow
    mask = make_nosnow_mask() #leaves 8659 columns

    z_era_m = apply_mask(agg_indices(z.era5, dates.era5, :month, mean)[1], mask)
    z_era_a = apply_mask(agg_indices(z.era5, dates.era5, :snowyear, maximum)[1], mask)

    metrics = Dict()
    std_era = time_stat(stdf, z_era_a)
    metrics["era5"] = Dict("IAV" => Dict("stats" => stat_spread(std_era), "vals" => std_era))

    #plot_min, plot_max = -2, 25
    for tag in [:me, :old, :others, :clm5]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        if tag != :clm5
            m_tidx, era_tidx = 1:243, 1:243
            print("   Aggregating Data...\n")
            z_month = apply_mask(agg_indices(z[tag], dates.clima, :month, mean)[1], mask)
            z_annual = apply_mask(agg_indices(z[tag], dates.clima, :year, maximum)[1], mask)
            std_yr = time_stat(stdf, z_annual)
            print("   Assessing IAV...\n")
            metrics[string(tag)]["IAV"] = Dict("stats" => stat_spread(std_yr), "vals" => std_yr)
        else
            m_tidx, era_tidx = Colon(), 71:130
            z_month = z.clm5 #no annual because we don't have daily outputs to work with; no additional monthly agg step needed
        end
        for (stat_f, f) in calcs
            print("   Assessing $(string(f))...\n")
            calc_d = stat_f(f, z_month[m_tidx, :, :], z_era_m[era_tidx, :, :])
            metrics[string(tag)][string(f)] = Dict("stats" => stat_spread(calc_d), "vals" => calc_d)
            #Hypothesis Test:
            if tag != :me
                print("    + hypothesis test\n")
                my_vals = metrics["me"][string(f)]["vals"][:]
                string(stat_f) == "space_stat" && (my_vals = my_vals[era_tidx]) #get the smaller temporal subset for comparison if this is clm5
                these_vals = calc_d[:]
                filt = is_data.(my_vals) .& is_data.(these_vals)
                test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                metrics[string(tag)][string(f)]["p"] = pvalue(test)
            end
        end
    end
    # add a mask excluding points where there is no snowpack, or not?
    # split up these errors by season? by year?
    # update your model with GPU fixes - should we start all from the era5Land initial state?

    # check anderson parameters again (which to use? or don't compare to it; just compare to CLM5)
    # should we be fine-tuning/calibrating? see what happens when you don't use an over-fit neural depth model?

    # do we need to divide by the snow cover fraction here? what's the right comparison?
    # make a histogram of the rmse errors by scheme betweeen 0.004 and 1 m for paper? or just you
    return metrics
end

function globe_heatmap(data; lat = clima_lat, lon = clima_lon)
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

#########################################################

#Utilities for printing in the terminal on the cluster:

function normalize(val; vmin, vmax, vcenter=nothing, stretch=:linear, gamma=1.0)
    if vcenter !== nothing
        dmax = max(abs(vmax - vcenter), abs(vmin - vcenter))
        x = clamp((val - vcenter) / dmax, -1, 1)
        centered = true
    else
        x = clamp((val - vmin) / (vmax - vmin), 0, 1)
        centered = false
    end

    # apply stretch
    if stretch == :linear
        y = x
    elseif stretch == :power
        y = centered ? sign(x) * abs(x)^gamma : x^gamma
    elseif stretch == :symlog
        y = sign(x) * log1p(abs(x)) / log(2)
    else
        error("Unknown stretch: $stretch")
    end

    return centered ? 0.5 * (y + 1) : clamp(y, 0, 1)
end


function colorprint(text, t; vmin=0.0, vmax=1.0, cmap=cgrad(:viridis))
    t = clamp(t, 0, 1)
    col = cmap[t]  # returns RGB color
    r = round(Int, 255 * col.r)
    g = round(Int, 255 * col.g)
    b = round(Int, 255 * col.b)
    print("\e[38;2;$(r);$(g);$(b)m$text\e[0m")
end

function plotm(m; mask = nothing,
    c=:viridis,
    vcenter=nothing,
    stretch=:linear,
    gamma=1.0,
    vmin = nothing,
    vmax = nothing,
    )
    to_plot = isnothing(mask) ? reverse(m; dims = 2) : reverse(m .* mask, dims = 2) #reverse along latitude
    data_max = maximum(x for x in skipmissing(to_plot) if !isnan(x))
    data_min = minimum(x for x in skipmissing(to_plot) if !isnan(x))

    cmin = isnothing(vmin) ? data_min : vmin
    cmax = isnothing(vmax) ? data_max : vmax

    cmap = cgrad(c)

    for i in 1:size(to_plot)[2]
        for j in 1:size(to_plot)[1]
            if ismissing(to_plot[j, i]) || isnan(to_plot[j, i])
                print("  ")
            else
                val = to_plot[j, i]
                t = normalize(
                    val;
                    vmin=cmin,
                    vmax=cmax,
                    vcenter=vcenter,
                    stretch=stretch,
                    gamma=gamma
                )
                colorprint("██", t; cmap=cmap)
            end
        end
        println()
    end
    println()
    print_colorbar(cmin, cmax;
        cmap=cmap,
        vcenter=vcenter,
        stretch=stretch,
        gamma=gamma
    )
    print("\n\n")
end

function print_colorbar(vmin, vmax; cmap = cgrad(:viridis), n=100, ticks=5, vcenter=nothing, stretch=:linear, gamma=1.0)
    vals = range(vmin, vmax; length=n)

    for val in vals
        t = normalize(val;
            vmin=vmin,
            vmax=vmax,
            vcenter=vcenter,
            stretch=stretch,
            gamma=gamma
        )
        col = cmap[t]
        r = round(Int, 255 * col.r)
        g = round(Int, 255 * col.g)
        b = round(Int, 255 * col.b)
        print("\e[48;2;$(r);$(g);$(b)m \e[0m")
    end
    println()

    tick_positions = round.(Int, LinRange(1, n, ticks))
    tick_values = round.(LinRange(vmin, vmax, ticks), digits=2)

    tick_line = [" " for _ in 1:n]
    for pos in tick_positions
        tick_line[pos] = "|"
    end

    println(join(tick_line))
    println(join([rpad(string(val), 24) for val in tick_values]))
end
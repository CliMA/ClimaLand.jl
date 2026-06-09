using NCDatasets, DataFrames, Statistics, CairoMakie, GeoMakie, Dates, JLD2, HypothesisTests
using ClimaLand

const my_param_dir = "/home/acharbon/outputs/all_global_daily"
const old_param_dir = "/home/acharbon/outputs/none_global_daily"
const htessel_param_dir = "/home/acharbon/outputs/others_global_daily"
const obs_param_dir = "/home/acharbon/outputs/era5_global"
const clm5_dir = "/home/acharbon/outputs/clm5_global"
const output_dir = "/home/acharbon/stats"

function combine_clima_data(path)
    filenames = readdir(path)
    data = Dict()
    for file in filenames
        if file[end-2:end] == ".nc"
            varname = Symbol(split(file, "_")[1])
            data[varname] = NCDataset(joinpath(path, file))[varname]
        end
    end
    return data
end

const me_outputs = combine_clima_data(my_param_dir)
const old_outputs = combine_clima_data(old_param_dir)
const others_outputs = combine_clima_data(htessel_param_dir)
const clm5_outputs = NCDataset(joinpath(clm5_dir, "all_clm_outputs.nc"))
const era5_rad_outputs = NCDataset(joinpath(obs_param_dir, "era_means_monthly_rad.nc")) #"/home/acharbon/my_clima_api_rad.nc/era_means_monthly_rad.nc")
const era5_swu = permutedims(era5_rad_outputs[:swd][:,:,:] .- era5_rad_outputs[:swn][:,:,:], (3, 1, 2)); #circshift (180, 0, 0) away from clima artifact
const era5_lwu = permutedims(era5_rad_outputs[:lwd][:,:,:] .- era5_rad_outputs[:lwn][:,:,:], (3, 1, 2)); #circshift (180, 0, 0) away from clima artifact
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
    #add the four glacial points not in mask:
    for (p1, p2) in [(111, 57), (261, 133), (316, 153), (302, 158)]
        ret_mask[p1, p2] = NaN
    end
    return ret_mask
end

function missing2NaN(arr)
    out = copy(arr)
    out[ismissing.(out)] .= NaN
    return out
end

function apply_mask(arr, mask::Array)
    if length(size(mask)) == 2
        @assert size(arr)[2:3] == size(mask)
        return missing2NaN(arr) .* reshape(mask, 1, size(mask)...)
    elseif length(size(mask)) == 3
        @assert size(mask) == size(arr)
        return missing2NaN(arr) .* mask
    end
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
    if "nbnd" in dims #albedo for clm5
        dat_avg = 0.5 .* (field[:, :, 1, :] .+ field[:, :, 2, :])
        return permutedims(dat_avg, (3, 1, 2))
    else
        dat = Array(permutedims(field, (time_idx, lon_idx, lat_idx)))
        return dat
    end
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
        return year.(dates .+ Dates.Day(122)) #Sept 1 of year 200X is snow year 200X + 1
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

function do_time_agg(agg, m)
    nt, n1, n2 = size(m)
    out = fill(NaN, n1, n2)
    @inbounds for i in 1:n1
        @inbounds for j in 1:n2
            out[i, j] = agg(view(m, :, i, j))
        end
    end
    return out
end

function agg_time_indices(data::Array, dates::Vector{<:Union{Date, DateTime}}, period::Symbol, agg; bin::Bool = false)
    nt, nx, ny = size(data)
    agg_ids = get_agg_ids(dates, period, bin)
    if isnothing(agg_ids)
        return data, dates
    else
        idset = sort(unique(agg_ids))
        dat = [fill(NaN, nx, ny) for _ in 1:length(idset)]
        for (i, id) in enumerate(idset)
            idxs = findall(==(id), agg_ids)
            dat[i] .= do_time_agg(agg, view(data, idxs, :, :))
        end
        return stack(dat, dims = 1), idset
    end
end

is_data(x) = !ismissing(x) && !isnan(x)
function meanf(v)
    n, tot = 0, 0.0
    @inbounds for x in v
        if is_data(x)
            n += 1
            tot += x
        end
    end
    n == 0 && return NaN
    return tot/n
end
function wmeanf(v, w)
    num, denom = 0.0, 0.0
    @inbounds for i in eachindex(v, w)
        xi = v[i]
        wi = w[i]
        if is_data(xi)
            num += xi*wi
            denom += wi
        end
    end
    return denom == 0 ? NaN : num/denom
end
function nonzero_meanf(v)
    tot, n = 0.0, 0
    @inbounds for x in v
        if is_data(x) && x >= 0.001
            n += 1
            tot += x
        end
    end
    n == 0 && return NaN
    return tot/n
end
function stdf(v)
    n, mean, M2 = 0, 0.0, 0.0
    @inbounds for x in v
        if is_data(x)
            n += 1
            delta = x - mean
            mean += delta / n
            M2 += delta * (x - mean)
        end
    end
    return n > 1 ? sqrt(M2 / n) : NaN
end
function maxf(v)
    m = -Inf
    found = false
    @inbounds for x in v
        if is_data(x) && x > m
            m = x
            found = true
        end
    end
    !found && return NaN
    return m
end
function minf(v)
    m = Inf
    found = false
    @inbounds for x in v
        if is_data(x) && x < m
            m = x
            found = true
        end
    end
    !found && return NaN
    return m
end
function medianf(v)
    xs = [x for x in v if is_data(x)]
    isempty(xs) && return NaN
    return Statistics.median(xs)
end
function quantilef(v, q::AbstractFloat)
    xs = [x for x in v if is_data(x)]
    isempty(xs) && return NaN
    return Statistics.quantile(xs, q)
end
mae(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(abs.(sim .- obs))
mse(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(abs2.(sim .- obs))
rmse(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = sqrt(meanf(abs2.(sim .- obs)))
bias(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(sim .- obs)
bias_sq(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = abs2(meanf(sim .- obs))
mae_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(abs.((sim .- obs) ./ obs))
rel_mae(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = mae(sim, obs) / nonzero_meanf(obs)
rmse_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = sqrt(meanf(abs2.((sim .- obs) ./ obs)))
mse_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(abs2.((sim .- obs) ./ obs))
bias_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf((sim .- obs) ./ obs)
bias_perc_sq(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = abs2(meanf((sim .- obs) ./ obs))
function nse(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray})
    idxs = is_data.(sim) .& is_data.(obs)
    count(idxs) == 0 && return NaN
    μ = mean(obs[idxs])
    return 1 - sum(abs2.(sim[idxs] .- obs[idxs]))/sum(abs2.(obs[idxs] .- μ))
end
function r2(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray})
    idxs = is_data.(sim) .& is_data.(obs)
    count(idxs) <= 2 && return NaN
    return cor(sim[idxs], obs[idxs])^2
end
spatial_r2(sim::Union{Matrix, SubArray}, obs::Union{SubArray,Matrix}) = r2(sim[:], obs[:])
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
function confusion_matrix(arr1, arr2, thresh = 0.5)
    @assert size(arr1) == size(arr2)
    T, LON, LAT = size(arr1)
    cm = fill(0, T, 2, 2)
    for i in 1:T
        for j in 1:LON
            for k in 1:LAT
                x1 = arr1[i, j, k]
                x2 = arr2[i, j, k]
                if is_data(x1) && is_data(x2)
                    if x1 >= thresh && x2 >= thresh
                        cm[i, 1, 1] += 1
                    elseif x1 >= thresh && x2 < thresh
                        cm[i, 1, 2] += 1
                    elseif x1 < thresh && x2 >= thresh
                        cm[i, 2, 1] += 1
                    else
                        cm[i, 2, 2] += 1
                    end
                end
            end
        end
    end
    return cm
end
function binary_stats(sim, obs; thresh = 0.5)
    cm = confusion_matrix(sim, obs, thresh)
    out = Dict()
    TP = view(cm, :, 1, 1)
    FP = view(cm, :, 1, 2)
    TN = view(cm, :, 2, 2)
    FN = view(cm, :, 2, 1)
    out["accuracy"] = Dict()
    out["accuracy"]["vals"] = (TP .+ TN) ./ (TP .+ FP .+ TN .+ FN)
    out["accuracy"]["stats"] = stat_spread(out["accuracy"]["vals"])

    out["precision"] = Dict()
    out["precision"]["vals"] = TP ./ (TP .+ FP)
    out["precision"]["stats"] = stat_spread(out["precision"]["vals"])

    out["recall"] = Dict()
    out["recall"]["vals"] = TP ./ (TP .+ FN)
    out["recall"]["stats"] = stat_spread(out["recall"]["vals"])

    out["f1"] = Dict()
    out["f1"]["vals"] = 2 .* TP ./ (2 .* TP .+ FP .+ FN)
    out["f1"]["stats"] = stat_spread(out["f1"]["vals"])

    out["csi"] = Dict()
    out["csi"]["vals"] = TP ./ (TP .+ FN .+ FP)
    out["csi"]["stats"] = stat_spread(out["csi"]["vals"])
    return out
end
function time_stat(f, arr)
    T, LON, LAT = size(arr)
    out = fill(NaN, LON, LAT)
    @inbounds for i in 1:LON
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
    @inbounds for i in 1:LON1
        for j in 1:LAT1
            out[i, j] = f(view(arr1, :, i, j), view(arr2, :, i, j))
        end
    end
    return out
end
function space_stat(f, arr)
    T, LON, LAT = size(arr)
    out = fill(NaN, T)
    @inbounds for i in 1:T
        out[i] = f(view(arr, i, :, :))
    end
    return out
end
function space_stat(f, arr1, arr2)
    T1, LON1, LAT1 = size(arr1)
    T2, LON2, LAT2 = size(arr2)
    @assert (LAT1 == LAT2) && (LON1 == LON2) && (T1 == T2)
    out = fill(NaN, T1)
    @inbounds for i in 1:T1
        out[i] = f(view(arr1, i, :, :), view(arr2, i, :, :))
    end
    return out
end

const grid_weights = collect(transpose(hcat([cosd.((-90:89)) for _ in 1:360]...)))
function stat_spread(v; globe_int = false)
    data = Dict(
        "min" => minf(v),
        "max" => maxf(v),
        "median" => medianf(v),
        "mean" => meanf(v),
        "q25" => quantilef(v, 0.25),
        "q75" => quantilef(v, 0.75),
        "std" => stdf(v)
    )
    if globe_int
        data["grid_mean"] = data["mean"]
        delete!(data, "mean")
        data["global_int_mean"] = wmeanf(v, grid_weights)
    end
    return data
end

function make_nosnow_mask(; swe_thresh = 0.001, glacial_thresh = 0) #what threshold? I think this is fine? 0.004 to get rid of 0-peak
    era5_snow_depth = get_data(era5_snow_outputs[:swe])
    max_snow_cell = maximum(era5_snow_depth, dims = 1)[1, :, :]
    max_snow_cell[ismissing.(max_snow_cell)] .= NaN
    no_snow = findall(<(swe_thresh), max_snow_cell)
    mask = glacial_mask(glacial_thresh)
    mask[no_snow] .= NaN
    return mask
end

function depth_analysis(; args = snow_args)
    @info "Running Depth Analysis..."
    print("   Extracting data...\n")
    z, dates = field_data(:snd, :snd, :SNOWDP_month; args = args)
    calcs = [(time_stat, rel_mae), (time_stat, bias), (time_stat, bias_sq), (time_stat, rmse), (time_stat, mse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    mask = make_nosnow_mask()

    z_era_m = apply_mask(agg_time_indices(z.era5, dates.era5, :month, meanf)[1], mask)
    z_era_a = apply_mask(agg_time_indices(z.era5, dates.era5, :snowyear, maxf)[1], mask)

    metrics = Dict()
    std_era = time_stat(stdf, z_era_a)
    metrics["era5"] = Dict("IAV" => Dict("stats" => stat_spread(std_era), "vals" => std_era))

    for tag in [:me, :old, :others, :clm5]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        if tag != :clm5
            m_tidx, era_tidx = 4:243, 4:243
            print("   Aggregating Data...\n")
            z_month = apply_mask(agg_time_indices(z[tag], dates.clima, :month, meanf)[1], mask)
            z_annual = apply_mask(agg_time_indices(z[tag], dates.clima, :snowyear, maxf)[1], mask)
            std_yr = time_stat(stdf, z_annual)
            print("   Assessing IAV...\n")
            metrics[string(tag)]["IAV"] = Dict("stats" => stat_spread(std_yr), "vals" => std_yr)
        else
            m_tidx, era_tidx = Colon(), 71:130
            z_month = z.clm5
        end
        for (stat_f, f) in calcs
            print("   Assessing $(string(f))...\n")
            calc_d = stat_f(f, z_month[m_tidx, :, :], z_era_m[era_tidx, :, :])
            sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true) : stat_spread(calc_d)
            metrics[string(tag)][string(f)] = Dict("stats" => sspread, "vals" => calc_d)
            if tag != :me
                print("    + hypothesis test\n")
                my_vals = metrics["me"][string(f)]["vals"][:]
                if string(stat_f) == "space_stat" && tag == :clm5
                    my_vals = my_vals[68:127]
                end
                these_vals = calc_d[:]
                filt = is_data.(my_vals) .& is_data.(these_vals)
                test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                metrics[string(tag)][string(f)]["p"] = pvalue(test)
            end
        end
    end
    metrics["upgrade_diff"] = Dict()
    model_diff = time_stat(meanf, apply_mask(z.me .- z.old, mask))
    metrics["upgrade_diff"]["vals"] = model_diff
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff, globe_int = true)
    return metrics
end

function swe_analysis(; args = snow_args)
    @info "Running SWE Analysis..."
    print("   Extracting data...\n")
    swe, dates = field_data(:swe, :swe, :SNOWICE_month; args = args)
    swe.clm5 .= (swe.clm5 .+ get_data(args.clm5[:SNOWLIQ_month])) ./ 1000
    calcs = [(time_stat, rel_mae), (time_stat, bias), (time_stat, bias_sq), (time_stat, rmse), (time_stat, mse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    mask = make_nosnow_mask()

    swe_era_m = apply_mask(agg_time_indices(swe.era5, dates.era5, :month, meanf)[1], mask)
    swe_era_a = apply_mask(agg_time_indices(swe.era5, dates.era5, :snowyear, maxf)[1], mask)

    metrics = Dict()
    std_era = time_stat(stdf, swe_era_a)
    metrics["era5"] = Dict("IAV" => Dict("stats" => stat_spread(std_era), "vals" => std_era))

    for tag in [:me, :old, :others, :clm5]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        if tag != :clm5
            m_tidx, era_tidx = 4:243, 4:243
            print("   Aggregating Data...\n")
            swe_month = apply_mask(agg_time_indices(swe[tag], dates.clima, :month, meanf)[1], mask)
            swe_annual = apply_mask(agg_time_indices(swe[tag], dates.clima, :snowyear, maxf)[1], mask)
            std_yr = time_stat(stdf, swe_annual)
            print("   Assessing IAV...\n")
            metrics[string(tag)]["IAV"] = Dict("stats" => stat_spread(std_yr), "vals" => std_yr)
        else
            m_tidx, era_tidx = Colon(), 71:130
            swe_month = swe.clm5
        end
        for (stat_f, f) in calcs
            print("   Assessing $(string(f))...\n")
            calc_d = stat_f(f, swe_month[m_tidx, :, :], swe_era_m[era_tidx, :, :])
            sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true) : stat_spread(calc_d)
            metrics[string(tag)][string(f)] = Dict("stats" => sspread, "vals" => calc_d)
            if tag != :me
                print("    + hypothesis test\n")
                my_vals = metrics["me"][string(f)]["vals"][:]
                if string(stat_f) == "space_stat" && tag == :clm5
                    my_vals = my_vals[68:127]
                end
                these_vals = calc_d[:]
                filt = is_data.(my_vals) .& is_data.(these_vals)
                test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                metrics[string(tag)][string(f)]["p"] = pvalue(test)
            end
        end
    end
    metrics["upgrade_diff"] = Dict()
    model_diff = time_stat(meanf, apply_mask(swe.me .- swe.old, mask))
    metrics["upgrade_diff"]["vals"] = model_diff
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff, globe_int = true)
    return metrics
end

function snow_albedo_masks(scf_data; scf_thresh = 0.05)
    masks = Dict(
        name => fill(NaN, size(scf_data[name])...) for name in propertynames(scf_data)
    )
    is_snow(v) = !ismissing(v) && !isnan(v) && v >= scf_thresh
    for name in propertynames(scf_data)
        masks[name][findall(is_snow, scf_data[name])] .= 1.0
    end
    return NamedTuple(masks)
end

function snow_alb_analysis(; args = snow_args)
    @info "Running Snow Albedo Analysis..."
    print("   Extracting data...\n")
    snalb, dates = field_data(:snalb, :snalb, :ALBSN_HIST_month; args = args)
    snalb = (snalb..., modis_1 = get_data(args.era5["bsa_1"]), modis_2 = get_data(args.era5["bsa_2"]))
    calcs = [(time_stat, rel_mae), (time_stat, bias), (time_stat, bias_sq), (time_stat, rmse), (time_stat, mse), (time_stat, mae), (time_stat, mae_perc), (time_stat, nse), (time_stat, bias_perc), (time_stat, bias_perc_sq), (time_stat, rmse_perc), (time_stat, mse_perc), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    
    #what threshold of scf to compare snow albedo? 0.95?
    scf, _ = field_data(:snowc, :snowc, :FSNO_EFF_month, args = args) #FNO or FNO_EFF?
    scf.era5 ./= 100
    scf = (scf..., modis_1 = get_data(args.era5["scf_1"]) ./ 100, modis_2 = get_data(args.era5["scf_2"]) ./ 100)
    masks = snow_albedo_masks(scf, scf_thresh = 0.95)
    for name in propertynames(snalb)
        snalb[name] .= apply_mask(snalb[name], masks[name])
    end

    snalb_era_m = agg_time_indices(snalb.era5, dates.era5, :month, meanf)[1]
    modis_1_m = agg_time_indices(snalb.modis_1, dates.era5, :month, meanf)[1]
    modis_2_m = agg_time_indices(snalb.modis_2, dates.era5, :month, meanf)[1]

    metrics = Dict()
    print("   Gathering IAV stats for all comparison data...\n")
    for tag in [:era5, :modis_1, :modis_2]
        snalb_a = agg_time_indices(snalb[tag], dates.era5, :snowyear, maxf)[1] #maximum for this one, or a different one?
        std_d = time_stat(stdf, snalb_a)
        metrics[string(tag)] = Dict("IAV" => Dict("stats" => stat_spread(std_d), "vals" => std_d))
    end

    for tag in [:me, :old, :others, :clm5]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        if tag != :clm5
            m_tidx, era_tidx = 4:243, 4:243
            print("   Aggregating Data...\n")
            snalb_month = agg_time_indices(snalb[tag], dates.clima, :month, meanf)[1]
            snalb_annual = agg_time_indices(snalb[tag], dates.clima, :snowyear, maxf)[1] #maximum for this one, or a different one?
            std_yr = time_stat(stdf, snalb_annual)
            print("   Assessing IAV...\n")
            metrics[string(tag)]["IAV"] = Dict("stats" => stat_spread(std_yr), "vals" => std_yr)
        else
            m_tidx, era_tidx = Colon(), 71:130
            snalb_month = snalb.clm5
        end
        for (dlabel, compare_data) in [("era5", snalb_era_m), ("modis_1", modis_1_m), ("modis_2", modis_2_m)]
            metrics[string(tag)][string(dlabel)] = Dict()
            for (stat_f, f) in calcs
                print("   Assessing $(string(f)) in $(dlabel)...\n")
                calc_d = stat_f(f, snalb_month[m_tidx, :, :], compare_data[era_tidx, :, :])
                sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true) : stat_spread(calc_d)
                metrics[string(tag)][string(dlabel)][string(f)] = Dict("stats" => sspread, "vals" => calc_d)
                if tag != :me
                    print("    + hypothesis test\n")
                    my_vals = metrics["me"][string(dlabel)][string(f)]["vals"][:]
                    if string(stat_f) == "space_stat" && tag == :clm5
                        my_vals = my_vals[68:127]
                    end
                    these_vals = calc_d[:]
                    filt = is_data.(my_vals) .& is_data.(these_vals)
                    test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                    metrics[string(tag)][string(dlabel)][string(f)]["p"] = pvalue(test)
                end
            end
        end
    end
    metrics["upgrade_diff"] = Dict()
    model_diff = time_stat(meanf, snalb.me .- snalb.old)
    metrics["upgrade_diff"]["vals"] = model_diff
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff, globe_int = true)
    return metrics
end

function surface_alb_analysis(; args = snow_args)
    @info "Running Surface Albedo Analysis..."
    print("   Extracting data...\n")
    swa, dates = field_data(:rpar, :swa, :ALB_HIST_month; args = args)
    swa = (swa..., modis_1 = get_data(args.era5["bsa_1"]), modis_2 = get_data(args.era5["bsa_2"]))
    swa.me .= (swa.me .+ get_data(args.me[:rnir])) ./ 2
    swa.old .= (swa.old .+ get_data(args.old[:rnir])) ./ 2
    swa.others .= (swa.others .+ get_data(args.others[:rnir])) ./ 2

    calcs = [(time_stat, rel_mae), (time_stat, bias), (time_stat, bias_sq), (time_stat, rmse), (time_stat, mse), (time_stat, mae), (time_stat, mae_perc), (time_stat, nse), (time_stat, bias_perc), (time_stat, bias_perc), (time_stat, rmse_perc), (time_stat, mse_perc), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    
    swa_era_m = agg_time_indices(swa.era5, dates.era5, :month, meanf)[1]
    modis_1_m = agg_time_indices(swa.modis_1, dates.era5, :month, meanf)[1]
    modis_2_m = agg_time_indices(swa.modis_2, dates.era5, :month, meanf)[1]

    metrics = Dict()
    print("   Gathering IAV stats for all comparison data...\n")
    for tag in [:era5, :modis_1, :modis_2]
        swa_a = agg_time_indices(swa[tag], dates.era5, :snowyear, maxf)[1] #maximum/:snowyear for this one, something different?
        std_d = time_stat(stdf, swa_a)
        metrics[string(tag)] = Dict("IAV" => Dict("stats" => stat_spread(std_d), "vals" => std_d))
    end

    for tag in [:me, :old, :others, :clm5]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        if tag != :clm5
            m_tidx, era_tidx = 4:243, 4:243
            print("   Aggregating Data...\n")
            swa_month = agg_time_indices(swa[tag], dates.clima, :month, meanf)[1]
            swa_annual = agg_time_indices(swa[tag], dates.clima, :snowyear, maxf)[1] #maximum/:snowyear for this one, something different?
            std_yr = time_stat(stdf, swa_annual)
            print("   Assessing IAV...\n")
            metrics[string(tag)]["IAV"] = Dict("stats" => stat_spread(std_yr), "vals" => std_yr)
        else
            m_tidx, era_tidx = Colon(), 71:130
            swa_month = swa.clm5
        end
        for (dlabel, compare_data) in [("era5", swa_era_m), ("modis_1", modis_1_m), ("modis_2", modis_2_m)]
            metrics[string(tag)][string(dlabel)] = Dict()
            for (stat_f, f) in calcs
                print("   Assessing $(string(f)) in $(dlabel)...\n")
                calc_d = stat_f(f, swa_month[m_tidx, :, :], compare_data[era_tidx, :, :])
                sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true) : stat_spread(calc_d)
                metrics[string(tag)][string(dlabel)][string(f)] = Dict("stats" => sspread, "vals" => calc_d)
                if tag != :me
                    print("    + hypothesis test\n")
                    my_vals = metrics["me"][string(dlabel)][string(f)]["vals"][:]
                    if string(stat_f) == "space_stat" && tag == :clm5
                        my_vals = my_vals[68:127]
                    end
                    these_vals = calc_d[:]
                    filt = is_data.(my_vals) .& is_data.(these_vals)
                    test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                    metrics[string(tag)][string(dlabel)][string(f)]["p"] = pvalue(test)
                end
            end
        end
    end
    metrics["upgrade_diff"] = Dict()
    model_diff = time_stat(meanf, swa.me .- swa.old)
    metrics["upgrade_diff"]["vals"] = model_diff
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff, globe_int = true)
    return metrics
end

function rad_data(rad_type)
    if !(rad_type in [:swu, :lwu, :shf, :lhf])
        error("this rad type not defined")
    end
    rad_me = get_data(me_outputs[rad_type])
    rad_old = get_data(old_outputs[rad_type])
    rad_others = get_data(others_outputs[rad_type])

    if rad_type == :swu
        rad_era5 = era5_swu
    elseif rad_type == :lwu
        rad_era5 = era5_lwu
    else
        rad_era5 = get_data(era5_rad_outputs[rad_type])
    end
    
    scale = rad_type in [:swu, :lwu] ? 1 : -1
    data = (me = rad_me, old = rad_old, others = rad_others, era5 = scale .* rad_era5)
    dates = (
        clima = me_outputs[rad_type][:date][:],
        era5 = era5_rad_outputs[:date][:],
    )
    return data, dates
end

function radiation_analysis(rad_type)
    @info "Running Flux Analysis ($(rad_type))..."
    
    print("   Extracting data...\n")
    rad, dates = rad_data(rad_type)

    calcs = [(time_stat, rel_mae), (time_stat, bias), (time_stat, bias_sq), (time_stat, rmse), (time_stat, mse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]

    metrics = Dict()
    for tag in [:me, :old, :others]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        m_tidx, era_tidx = 4:243, 6:245
        print("   Aggregating Data...\n")
        rad_month = agg_time_indices(rad[tag], dates.clima, :month, meanf)[1]
        for (stat_f, f) in calcs
            print("   Assessing $(string(f))...\n")
            calc_d = stat_f(f, rad_month[m_tidx, :, :], rad.era5[era_tidx, :, :])
            sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true) : stat_spread(calc_d)
            metrics[string(tag)][string(f)] = Dict("stats" => sspread, "vals" => calc_d)
            if tag != :me
                print("    + hypothesis test\n")
                my_vals = metrics["me"][string(f)]["vals"][:]
                these_vals = calc_d[:]
                filt = is_data.(my_vals) .& is_data.(these_vals)
                test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                metrics[string(tag)][string(f)]["p"] = pvalue(test)
            end
        end
    end
    metrics["upgrade_diff"] = Dict()
    model_diff = time_stat(meanf, rad.me .- rad.old)
    metrics["upgrade_diff"]["vals"] = model_diff
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff, globe_int = true)
    return metrics
end

function scf_analysis(; args = snow_args)
    @info "Running Snow Cover Analysis..."
    print("   Extracting data...\n")
    
    scf, dates = field_data(:snowc, :snowc, :FSNO_EFF_month, args = args) #FNO or FNO_EFF?
    scf.era5 ./= 100
    scf = (scf..., modis_1 = get_data(args.era5["scf_1"]) ./ 100, modis_2 = get_data(args.era5["scf_2"]) ./ 100)
    
    mask = make_nosnow_mask()
    calcs = [(time_stat, rel_mae), (time_stat, bias), (time_stat, bias_sq), (time_stat, rmse), (time_stat, mse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]

    scf_era_m = apply_mask(agg_time_indices(scf.era5, dates.era5, :month, meanf)[1], mask)
    modis_1_m = apply_mask(agg_time_indices(scf.modis_1, dates.era5, :month, meanf)[1], mask)
    modis_2_m = apply_mask(agg_time_indices(scf.modis_2, dates.era5, :month, meanf)[1], mask)

    metrics = Dict()
    print("   Gathering IAV stats for all comparison data...\n")
    for tag in [:era5, :modis_1, :modis_2]
        scf_a = apply_mask(agg_time_indices(scf[tag], dates.era5, :snowyear, maxf)[1], mask)
        std_d = time_stat(stdf, scf_a)
        metrics[string(tag)] = Dict("IAV" => Dict("stats" => stat_spread(std_d), "vals" => std_d))
    end

    for tag in [:me, :old, :others, :clm5]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        if tag != :clm5
            m_tidx, era_tidx = 4:243, 4:243
            print("   Aggregating Data...\n")
            scf_month = apply_mask(agg_time_indices(scf[tag], dates.clima, :month, meanf)[1], mask)
            scf_annual = apply_mask(agg_time_indices(scf[tag], dates.clima, :snowyear, maxf)[1], mask)
            std_yr = time_stat(stdf, scf_annual)
            print("   Assessing IAV...\n")
            metrics[string(tag)]["IAV"] = Dict("stats" => stat_spread(std_yr), "vals" => std_yr)
        else
            m_tidx, era_tidx = Colon(), 71:130
            scf_month = scf.clm5
        end
        for (dlabel, compare_data) in [("era5", scf_era_m), ("modis_1", modis_1_m), ("modis_2", modis_2_m)]
            metrics[string(tag)][string(dlabel)] = Dict()
            for (stat_f, f) in calcs
                print("   Assessing $(string(f)) in $(dlabel)...\n")
                calc_d = stat_f(f, scf_month[m_tidx, :, :], compare_data[era_tidx, :, :])
                sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true) : stat_spread(calc_d)
                metrics[string(tag)][string(dlabel)][string(f)] = Dict("stats" => sspread, "vals" => calc_d)
                if tag != :me
                    print("    + hypothesis test\n")
                    my_vals = metrics["me"][string(dlabel)][string(f)]["vals"][:]
                    if string(stat_f) == "space_stat" && tag == :clm5
                        my_vals = my_vals[68:127]
                    end
                    these_vals = calc_d[:]
                    filt = is_data.(my_vals) .& is_data.(these_vals)
                    test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                    metrics[string(tag)][string(dlabel)][string(f)]["p"] = pvalue(test)
                end
            end
            print("   Getting binary statistics...\n")
            metrics[string(tag)][string(dlabel)]["binary"] = binary_stats(scf_month[m_tidx, :, :], compare_data[era_tidx, :, :])
            if tag != :me
                print("    + hypothesis tests\n")
                for stat in keys(metrics[string(tag)][string(dlabel)]["binary"])
                    use_idxs = tag == :clm5 ? (68:127) : Colon()
                    my_vals = metrics["me"][string(dlabel)]["binary"][stat]["vals"][use_idxs]
                    these_vals = metrics[string(tag)][string(dlabel)]["binary"][stat]["vals"][:]
                    test = HypothesisTests.SignedRankTest(my_vals, these_vals);
                    metrics[string(tag)][string(dlabel)]["binary"][stat]["p"] = pvalue(test)
                end
            end
        end
    end
    metrics["upgrade_diff"] = Dict()
    model_diff = time_stat(meanf, apply_mask(scf.me .- scf.old, mask))
    metrics["upgrade_diff"]["vals"] = model_diff
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff, globe_int = true)
    return metrics
end

function get_blocks(cond; length_req = 0)
    all_indices = findall(cond)
    if isempty(all_indices) return UnitRange{Int}[] end
    gap_indices = findall(diff(all_indices) .> 1)
    starts = vcat(1, gap_indices .+ 1)
    ends = vcat(gap_indices, length(all_indices))
    return [all_indices[s]:all_indices[e] for (s, e) in zip(starts, ends) if all_indices[e] - all_indices[s] >= length_req - 1]
end

z_season_filter(x, z_thresh) = is_data(x) && x >= z_thresh
no_snow_filter(x, z_thresh) = is_data(x) && x < z_thresh

#=
function get_snow_timing(z, dates; z_thresh::AbstractFloat, ref_date::DateTime)
    snow_cycle = year.(dates .- Day(210))
    nyr = length(sort(unique(snow_cycle)))
    starts = fill(NaN, nyr)
    finishes = fill(NaN, nyr)
    maxdates = fill(NaN, nyr)
    maxzs = fill(NaN, nyr)
    for (i, yr) in enumerate(sort(unique(snow_cycle)))
        subdates = snow_cycle .== yr
        snow_idxs = findall(x -> z_season_filter(x, z_thresh), z[subdates])
        if isempty(snow_idxs) || length(snow_idxs) > 0.95*count(subdates)
            continue
        end
        first_snowday, last_snowday = snow_idxs[1], snow_idxs[end]
        starts[i] = Dates.value(Day(dates[subdates][first_snowday] - ref_date))
        finishes[i] = Dates.value(Day(dates[subdates][last_snowday] + Day(1) - ref_date))
        (zmax, maxidx) = findmax(z[subdates])
        maxdates[i] = Dates.value(Day(dates[subdates][maxidx] - ref_date))
        maxzs[i] = zmax
    end
    return starts, finishes, maxdates, maxzs
end

function process_snowpacks(z, dates; z_thresh = 0.05, min_days = 14, ignore_south = true)
    T, LON, LAT = size(z)
    snow_cycles = sort(unique(year.(dates .- Day(210))))
    NYR = length(snow_cycles)
    ref_date = DateTime("1999-01-01")
    start_dates = fill(NaN, NYR, LON, LAT)
    max_dates = fill(NaN, NYR, LON, LAT)
    max_z = fill(NaN, NYR, LON, LAT)
    end_dates = fill(NaN, NYR, LON, LAT)
    lat_range = ignore_south ? (91:LAT) : (1:LAT)
    for i in 1:LON
        for j in lat_range
            starts, finishes, maxdates, maxzs = get_snow_timing(view(z, :, i, j), dates; z_thresh = z_thresh, ref_date = ref_date)
            start_dates[:, i, j] .= starts
            max_dates[:, i, j] .= maxdates
            max_z[:, i, j] .= maxzs
            end_dates[:, i, j] .= finishes
        end
    end
    return start_dates, end_dates, max_dates, max_z, ref_date
end
=#

function get_snow_timing(z, dates; z_thresh::AbstractFloat, min_days::Int, ref_date::DateTime, day_shift::Int = 210)
    snow_free_periods = get_blocks(no_snow_filter.(z, Ref(z_thresh)), length_req = 60)
    coverage_fraction = count(skipmissing(z .>= z_thresh)) / length(z)
    nyr = length(unique(year.(dates))) - 1
    if length(snow_free_periods) < nyr/2 || coverage_fraction > 0.95
        return nothing
    end
    snow_idxs = get_blocks(z_season_filter.(z, Ref(z_thresh)), length_req = min_days)
    isempty(snow_idxs) && return nothing
    maxs = [findmax(view(z, s)) for s in snow_idxs]
    out = [(
        start_date = Dates.value(Day(dates[s.start] - ref_date)),
        start_year = year(dates[s.start] - Day(day_shift)),
        end_date = Dates.value(Day(dates[s.stop] + Day(1) - ref_date)),
        end_year = year(dates[s.stop] - Day(day_shift - 1)),
        max_date = Dates.value(Day(dates[s][maxs[i][2]] - ref_date)),
        snow_cycle = year(dates[s][maxs[i][2]] - Day(day_shift)),
        max_z = maxs[i][1],
    ) for (i, s) in enumerate(snow_idxs)]
    return sort!(out, by = x -> x.end_date - x.start_date)
end

function process_snowpacks(z, dates; z_thresh = 0.05, min_days = 14, ignore_south = true) #better threshold values?
    T, LON, LAT = size(z)
    snow_cycles = sort(unique(year.(dates .- Day(210))))
    NYR = length(snow_cycles)
    ref_date = DateTime("1999-01-01")
    start_dates = fill(Inf, NYR, LON, LAT)
    max_dates = fill(NaN, NYR, LON, LAT)
    max_z = fill(-Inf, NYR, LON, LAT)
    end_dates = fill(-Inf, NYR, LON, LAT)
    yr2index = Dict(y => i for (i, y) in enumerate(snow_cycles))
    lat_range = ignore_south ? (91:LAT) : (1:LAT)
    for i in 1:LON
        for j in lat_range
            snowpacks = get_snow_timing(view(z, :, i, j), dates, z_thresh = z_thresh, min_days = min_days, ref_date = ref_date, day_shift = 210)
            isnothing(snowpacks) && continue
            for s in snowpacks
                yi = yr2index[s.snow_cycle]
                #yr_start = yr2index[s.start_year]
                #yr_end = yr2index[s.end_year]
                sd = start_dates[yi, i, j]
                ed = end_dates[yi, i, j]
                mz = max_z[yi, i, j]
                #if !isinf(ed) && s.start_date - ed > 60
                #    #more likely a starter snow for the next year's cycle:
                #    yi += 1
                #    yi > NYR && continue
                #    sd = start_dates[yi, i, j]
                #    end_dates[yi, i, j]
                #    max_z[yi, i, j]
                #end
                if s.start_date < sd #&& (sd - s.end_date) < 30
                    start_dates[yi, i, j] = s.start_date
                end
                if s.end_date > ed #&& (s.start_date - ed) < 30
                    end_dates[yi, i, j] = s.end_date
                end
                if s.max_z > mz
                    max_dates[yi, i, j] = s.max_date
                    max_z[yi, i, j] = s.max_z
                end
            end
        end
    end
    start_dates[isinf.(start_dates)] .= NaN
    end_dates[isinf.(end_dates)] .= NaN
    max_z[isinf.(max_z)] .= NaN
    return start_dates, end_dates, max_dates, max_z, ref_date
end

function phenology_analysis2(z, dates; z_thresh = 0.05, min_days = 14)  
    T, LON, LAT = size(z)
    ref_date = DateTime("1999-01-01")
    start_dates = Dict(l => [] for l in clima_lat)
    max_dates = Dict(l => [] for l in clima_lat)
    end_dates = Dict(l => [] for l in clima_lat)
    durations_strg = fill([], LON, LAT)
    durations = fill(NaN, LON, LAT)
    for i in 1:LON
        for j in 1:LAT
            snowpacks = get_snow_timing(view(z, :, i, j), dates, z_thresh = z_thresh, min_days = min_days, ref_date = ref_date, day_shift = 0)
            isnothing(snowpacks) && continue
            for s in snowpacks
                sd = Dates.dayofyear(ref_date + Day(s.start_date))
                md = Dates.dayofyear(ref_date + Day(s.max_date))
                ed = Dates.dayofyear(ref_date + Day(s.end_date))
                push!(start_dates[clima_lat[j]], sd)
                push!(max_dates[clima_lat[j]], md)
                push!(end_dates[clima_lat[j]], ed)
                push!(durations_strg[i, j], s.end_date - s.start_date)
            end
            length(durations_strg[i, j]) > 0 && (durations[i, j] = median(durations_strg[i, j]))
        end
    end
    return start_dates, end_dates, max_dates, durations
end

function phenology_analysis(; args = snow_args)
    @info "Running Snow Phenology Analysis..."
    print("   Extracting data...\n")

    z, dates = field_data(:snd, :snd, :SNOWDP_month; args = args)
    calcs = [(time_stat, bias), (time_stat, bias_sq), (time_stat, rmse), (time_stat, mse), (time_stat, mae), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    mask = make_nosnow_mask()

    for tag in propertynames(z)
        z[tag] .= apply_mask(z[tag], mask)
    end

    nt_me = size(z.me)[1]
    era_start, era_end, era_maxd, era_maxz, ref_date = process_snowpacks(
        z.era5[1:nt_me, :, :],
        dates.era5[1:nt_me];
        z_thresh = 0.05,
        min_days = 14,
        ignore_south = true
    )

    fields = ["start_date", "end_date", "max_date", "max_val", "duration"]
    metrics = Dict()
    metrics["ref_date"] = ref_date
    metrics["start_date"] = Dict("era5" => era_start)
    metrics["end_date"] = Dict("era5" => era_end)
    metrics["max_date"] = Dict("era5" => era_maxd)
    metrics["max_val"] = Dict("era5" => era_maxz)
    metrics["duration"] = Dict("era5" => era_end .- era_start)

    for tag in [:me, :old, :others]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        z_start, z_end, z_maxd, z_maxz, _ = process_snowpacks(z[tag], dates.clima; z_thresh = 0.05, min_days = 14, ignore_south = true)
        metrics["start_date"][string(tag)] = z_start
        metrics["end_date"][string(tag)] = z_end
        metrics["max_date"][string(tag)] = z_maxd
        metrics["max_val"][string(tag)] = z_maxz
        metrics["duration"][string(tag)] = z_end .- z_start
        for (feature, sim, obs) in [("start_date", z_start, era_start), ("end_date", z_end, era_end), ("max_date", z_maxd, era_maxd)]
            print("   Assessing $(string(feature))...\n")
            metrics[string(tag)][feature] = Dict()
            for (stat_f, f) in calcs
                print("    Assessing $(string(f))...\n")
                calc_d = stat_f(f, sim, obs)
                sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true) : stat_spread(calc_d)
                metrics[string(tag)][feature][string(f)] = Dict("stats" => sspread, "vals" => calc_d)
                if tag != :me
                    print("    + hypothesis test\n")
                    my_vals = metrics["me"][feature][string(f)]["vals"][:]
                    these_vals = calc_d[:]
                    filt = is_data.(my_vals) .& is_data.(these_vals)
                    test = HypothesisTests.SignedRankTest(disallowmissing(my_vals[filt]),disallowmissing(these_vals[filt]));
                    metrics[string(tag)][feature][string(f)]["p"] = pvalue(test)
                end
            end
        end
    end
    metrics["upgrade_diff"] = Dict()
    for field in fields
        model_diff = time_stat(meanf, metrics[field]["me"] .- metrics[field]["old"])
        metrics["upgrade_diff"][field] = Dict("vals" => model_diff, "stats" => stat_spread(model_diff, globe_int = true))
    end
    return metrics
end

function get_direct_analyses()
    return Dict(
        "z" => depth_analysis(),
        "swe" => swe_analysis(),
        "snalb" => snow_alb_analysis(),
        "swa" => surface_alb_analysis(),
        "scf" => scf_analysis(),
        "rad" => Dict(
            "lhf" => radiation_analysis(:lhf),
            "shf" => radiation_analysis(:shf),
            "swu" => radiation_analysis(:swu),
            "lwu" => radiation_analysis(:lwu),
        ),
        "phenom" => phenology_analysis(),
    )
end

function full_analysis(output_path; args = snow_args)
    analysis = get_direct_analyses()
    # - land-type? Is there a easy-to-grab PFT map?
    @info "Adding Meta-analysis information"
    data = Dict()
    data["elevation"] = elevs
    data["land_mask"] = clima_land_mask
    data["glacial_mask"] = glacial_mask(0)
    data["nosnow_mask"] = make_nosnow_mask()
    data["dates"] = Dict()
    print("  Grabbing z data to get nonzero snowpack means per pixel...\n")
    z, dates = field_data(:snd, :snd, :SNOWDP_month, args = args)
    data["dates"]["clm5"] = dates.clm5
    data["dates"]["clima"] = dates.clima
    data["dates"]["era5"] = dates.era5
    data["avg_z"] = Dict() 
    data["med_szn_dur"] = Dict()
    data["maximum_zs"] = Dict()
    mask = data["nosnow_mask"]
    for (tag1, tag2) in [(:me, :clima), (:old, :clima), (:others, :clima), (:era5, :era5), (:clm5, :clm5)]
        print("    - :$(tag1)\n")
        z_avg_nonzero = time_stat(nonzero_meanf, z[tag1])
        data["avg_z"][string(tag1)] = z_avg_nonzero
        z_annual_max = agg_time_indices(apply_mask(z[tag1], mask), dates[tag2], :year, maxf)[1]
        data["maximum_zs"][string(tag1)] = z_annual_max
    end
    print("  Grabbing scf data to get gantt chart timing...\n")
    scf, dates = field_data(:snowc, :snowc, :FSNO_EFF_month, args = args) #FNO or FNO_EFF?
    scf.era5 ./= 100
    scf = (scf..., modis_1 = get_data(args.era5["scf_1"]) ./ 100, modis_2 = get_data(args.era5["scf_2"]) ./ 100)
    data["scf_gantt"] = Dict()
    for (tag1, tag2) in [(:me, :clima), (:old, :clima), (:others, :clima), (:era5, :era5), (:modis_1, :era5), (:modis_2, :era5)]
        print("    - :$(tag1)\n")
        scf_reduced, _ = agg_time_indices(scf[tag1], dates[tag2], :week, meanf, bin = true)
        scf_gantt = fill(NaN, 52, 180)
        for i in 1:52
            for j in 1:180
                scf_gantt[i, j] = meanf(scf_reduced[i, :, j][:])
            end
        end
        data["scf_gantt"][string(tag1)] = scf_gantt
    end
    print("  getting phenology binning data:\n")
    data["lon_bin_dates"] = Dict()
    for (tag1, tag2) in [(:me, :clima), (:old, :clima), (:others, :clima), (:era5, :era5)]
        print("    - :$(tag1)\n")
        sd, ed, md, dur = phenology_analysis2(z[tag1], dates[tag2])
        data["lon_bin_dates"][string(tag1)] = Dict(
            "start_dates" => sd,
            "end_dates" => ed,
            "max_dates" => md,
        )
        data["med_szn_dur"][string(tag1)] = dur
    end
    full_data = Dict(
        "scoring" => analysis,
        "meta" => data,
    )
    print("   getting snow age data:\n")
    print("Writing file (size: $(round(Base.summarysize(full_data) / 1e6, digits = 2)) MB)\n")
    jldsave(output_path, data = full_data)
end

#full_analysis("./stats/first_stats.jld2")

# see what happens when you don't use an over-fit neural depth model?

# split up these errors by season? by year? by land type?
# update your model with GPU fixes

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

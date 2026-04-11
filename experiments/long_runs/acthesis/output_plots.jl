using NCDatasets, DataFrames, Statistics, CairoMakie, GeoMakie, Dates, JLD2, HypothesisTests
using ClimaLand

const my_param_dir = "/home/acharbon/outputs/all_global"
const old_param_dir = "/home/acharbon/outputs/none_global"
const htessel_param_dir = "/home/acharbon/outputs/others_global"
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
const era5_rad_outputs = NCDataset("/home/acharbon/my_clima_api_rad.nc/era_means_monthly_rad.nc") ##NCDataset("/home/acharbon/my_clima_api_rad.nc/era_means_monthly_rad.nc")
#^^SW_up is SW_down - SW_net, LW_up is LW_down - LW_net 
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
        return missing2NaN(arr .* reshape(mask, 1, size(mask)...))
    elseif length(size(mask)) == 3
        @assert size(mask) == size(arr)
        return missing2NaN(arr .* mask)
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
            dat[i] .= dat[i] = do_time_agg(agg, view(data, idxs, :, :))
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
rmse(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = sqrt(meanf(abs2.(sim .- obs)))
bias(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(sim .- obs)
mae_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf(abs.((sim .- obs) ./ obs))
rmse_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = sqrt(meanf(abs2.((sim .- obs) ./ obs)))
bias_perc(sim::Union{Vector, SubArray}, obs::Union{Vector, SubArray}) = meanf((sim .- obs) ./ obs)
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

function stat_spread(v; globe_int = false, weights = nothing)
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
        data["global_int_mean"] = wmeanf(v, weights)
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
    calcs = [(time_stat, bias), (time_stat, rmse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    #make data mask: glacial mask, plus spots where there is no snow
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
            m_tidx, era_tidx = 1:243, 1:243
            print("   Aggregating Data...\n")
            z_month = apply_mask(agg_time_indices(z[tag], dates.clima, :month, meanf)[1], mask)
            z_annual = apply_mask(agg_time_indices(z[tag], dates.clima, :year, maxf)[1], mask)
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
    #diff between old and new model:
    metrics["upgrade_diff"] = Dict()
    model_diff = time_stat(meanf, apply_mask(z.me .- z.old, mask))
    metrics["upgrade_diff"]["vals"] = model_diff
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff)
    # split up these errors by season? by year? see katherine's suggested plot (gantt timeline chart or stacked range bar chart)
    # update your model with GPU fixes - should we start all from the era5Land initial state? (no just skip a burn-in year or start at the summer when snow disappears)

    # check anderson parameters again (which to use? or don't compare to it; just compare to CLM5)
    # should we be fine-tuning/calibrating? see what happens when you don't use an over-fit neural depth model?
    return metrics
end

function swe_analysis(; args = snow_args)
    @info "Running SWE Analysis..."
    print("   Extracting data...\n")
    swe, dates = field_data(:swe, :swe, :SNOWICE_month; args = args)
    swe.clm5 .= (swe.clm5 .+ get_data(args.clm5[:SNOWLIQ_month])) ./ 1000
    calcs = [(time_stat, bias), (time_stat, rmse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
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
            m_tidx, era_tidx = 1:243, 1:243
            print("   Aggregating Data...\n")
            swe_month = apply_mask(agg_time_indices(swe[tag], dates.clima, :month, meanf)[1], mask)
            swe_annual = apply_mask(agg_time_indices(swe[tag], dates.clima, :year, maxf)[1], mask)
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
            metrics[string(tag)][string(f)] = Dict("stats" => stat_spread(calc_d), "vals" => calc_d)
            if tag != :me
                print("    + hypothesis test\n")
                my_vals = metrics["me"][string(f)]["vals"][:]
                string(stat_f) == "space_stat" && (my_vals = my_vals[era_tidx])
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
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff)
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
    calcs = [(time_stat, bias), (time_stat, rmse), (time_stat, mae), (time_stat, mae_perc), (time_stat, nse), (time_stat, bias_perc), (time_stat, rmse_perc), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    
    #make/apply data masks: only want to compare data at places where scf is above some thresh (0.95?):
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
            m_tidx, era_tidx = 1:243, 1:243
            print("   Aggregating Data...\n")
            snalb_month = agg_time_indices(snalb[tag], dates.clima, :month, meanf)[1]
            snalb_annual = agg_time_indices(snalb[tag], dates.clima, :year, maxf)[1] #maximum for this one, or a different one?
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
                metrics[string(tag)][string(dlabel)][string(f)] = Dict("stats" => stat_spread(calc_d), "vals" => calc_d)
                if tag != :me
                    print("    + hypothesis test\n")
                    my_vals = metrics["me"][string(dlabel)][string(f)]["vals"][:]
                    string(stat_f) == "space_stat" && (my_vals = my_vals[era_tidx])
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
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff)
    return metrics
end

function surface_alb_analysis(; args = snow_args)
    @info "Running Surface Albedo Analysis..."
    print("   Extracting data...\n")
    swa, dates = field_data(:rpar, :swa, :ALB_HIST_month; args = args)
    swa = (swa..., modis_1 = get_data(args.era5["bsa_1"]), modis_2 = get_data(args.era5["bsa_2"]))
    swa.me .= (swa.me .+ get_data(args.me[:rnir])) ./ 2

    calcs = [(time_stat, bias), (time_stat, rmse), (time_stat, mae), (time_stat, mae_perc), (time_stat, nse), (time_stat, bias_perc), (time_stat, rmse_perc), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]
    
    #make/apply data masks: no mask for this one, full globe:
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
            m_tidx, era_tidx = 1:243, 1:243
            print("   Aggregating Data...\n")
            swa_month = agg_time_indices(swa[tag], dates.clima, :month, meanf)[1]
            swa_annual = agg_time_indices(swa[tag], dates.clima, :year, maxf)[1] #maximum/:snowyear for this one, something different?
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
                metrics[string(tag)][string(dlabel)][string(f)] = Dict("stats" => stat_spread(calc_d), "vals" => calc_d)
                if tag != :me
                    print("    + hypothesis test\n")
                    my_vals = metrics["me"][string(dlabel)][string(f)]["vals"][:]
                    string(stat_f) == "space_stat" && (my_vals = my_vals[era_tidx])
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
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff)
    return metrics
end
    
function rad_data(rad_type)
    if !(rad_type in [:swn, :lwn, :shf, :lhf])
        error("this rad type not defined")
    end
    rad_me = get_data(me_outputs[rad_type])
    rad_old = get_data(old_outputs[rad_type])
    rad_others = get_data(others_outputs[rad_type])
    
    rad_era5 = get_data(era5_rad_outputs[rad_type])
    #need to multiply era5 lhf/shf data by -1 to get direction right
    scale = rad_type in [:swn, :lwn] ? 1 : -1
    data = (me = rad_me, old = rad_old, others = rad_others, era5 = scale .* rad_era5)
    dates = (
        clima = me_outputs[rad_type][:date][:],
        era5 = era5_rad_outputs[rad_type][:date][:],
    )
    return data, dates
end

# save the SW_down/SW_up and LW_down/LW_up diagnostics to compare nets to data
# Your comparison file is identical to within 3e-12 per element... the difference must manifest in the model outputs themselves!
function radiation_analysis(rad_type)
    @info "Running Flux Analysis ($(rad_type))..."
    
    print("   Extracting data...\n")
    rad, dates = rad_data(rad_type)
    #no mask on this data, whole globe
    #era5 already at monthly levels here, so no IAV from daily data
    #we will need to weight the output matrices by the cosine of their latitutde, as grid cells have different areas on globe
    grid_weights = collect(transpose(hcat([cosd.((-90:89)) for _ in 1:360]...)))
    calcs = [(time_stat, bias), (time_stat, rmse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]

    metrics = Dict()
    for tag in [:me, :old, :others]
        print("  Analyzing Model of class :$(tag)\n")
        metrics[string(tag)] = Dict()
        m_tidx, era_tidx = 1:243, 3:245
        print("   Aggregating Data...\n")
        rad_month = agg_time_indices(rad[tag], dates.clima, :month, meanf)[1]
        #no IAV to compare to, leave out annual agg for this one
        for (stat_f, f) in calcs
            print("   Assessing $(string(f))...\n")
            calc_d = stat_f(f, rad_month[m_tidx, :, :], rad.era5[era_tidx, :, :])
            sspread = string(stat_f) == "time_stat" ? stat_spread(calc_d, globe_int = true, weights = grid_weights) : stat_spread(calc_d)
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
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff, globe_int = true, weights = grid_weights)
    return metrics
end

function scf_analysis(; args = snow_args)
    @info "Running Snow Cover Analysis..."
    print("   Extracting data...\n")
    
    scf, dates = field_data(:snowc, :snowc, :FSNO_EFF_month, args = args) #FNO or FNO_EFF?
    scf.era5 ./= 100
    scf = (scf..., modis_1 = get_data(args.era5["scf_1"]) ./ 100, modis_2 = get_data(args.era5["scf_2"]) ./ 100)
    
    mask = make_nosnow_mask()
    calcs = [(time_stat, bias), (time_stat, rmse), (time_stat, mae), (time_stat, nse), (time_stat, r2), (space_stat, ssim), (space_stat, spatial_r2)]

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
            m_tidx, era_tidx = 1:243, 1:243
            print("   Aggregating Data...\n")
            scf_month = apply_mask(agg_time_indices(scf[tag], dates.clima, :month, meanf)[1], mask)
            scf_annual = apply_mask(agg_time_indices(scf[tag], dates.clima, :year, maxf)[1], mask)
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
                metrics[string(tag)][string(dlabel)][string(f)] = Dict("stats" => stat_spread(calc_d), "vals" => calc_d)
                if tag != :me
                    print("    + hypothesis test\n")
                    my_vals = metrics["me"][string(dlabel)][string(f)]["vals"][:]
                    string(stat_f) == "space_stat" && (my_vals = my_vals[era_tidx])
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
                    my_vals = metrics["me"][string(dlabel)]["binary"][stat]["vals"][era_tidx]
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
    metrics["upgrade_diff"]["stats"] = stat_spread(model_diff)
    return metrics
end

# phenology
function phenology_analysis(; args = snow_args)
    #use Z threshold to determine snow presence...how much? how many cm? 5cm?

end


# metaanalyses plots/calculations
# siteplots/stats

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

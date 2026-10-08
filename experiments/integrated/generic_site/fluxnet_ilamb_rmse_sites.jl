# Run the default ClimaLand model at FLUXNET2015 sites, one after another, and
# write the monthly means of modeled and observed LE, H, SWup and LWup of each
# site to `out/ilamb_rmse/sites/<SITE_ID>.csv`. Each site spins up for one year
# on its own forcing (up to the next month boundary) and is scored on the rest
# of its record.
# `fluxnet_ilamb_rmse_plot.jl` turns these files into ILAMB RMSEs.
#
# Without arguments, the sites are those of `ilamb_land_hist_site_rmse.csv`
# that are in the `fluxnet2015` artifact, which is too large to be
# downloadable (see ClimaArtifacts). Under Buildkite `parallelism`, the sites
# are split over the jobs by record length.
#
# Usage:
#   julia --project=.buildkite experiments/integrated/generic_site/fluxnet_ilamb_rmse_sites.jl [SITE_ID ...]

import ClimaComms
ClimaComms.@import_required_backends
using Dates
import DelimitedFiles
import NCDatasets
import ClimaDiagnostics
import ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations

include(joinpath(@__DIR__, "list_fluxnet_sites.jl"))

const FT = Float64
const SPINUP = Year(1)
# Sites with a shorter record after the spin-up are skipped
const MIN_SCORED = Year(1)
# ClimaLand short name => (ILAMB variable, column of the FLUXNET2015 monthly
# file); ILAMB scores these columns
const VARIABLES = (
    lhf = ("hfls", "LE_F_MDS"),
    shf = ("hfss", "H_F_MDS"),
    swu = ("rsus", "SW_OUT"),
    lwu = ("rlus", "LW_OUT"),
)
const OUTDIR = joinpath(@__DIR__, "out", "ilamb_rmse", "sites")

"""
    site_directory(site_ID)

Return the directory of `site_ID` in the `fluxnet2015` artifact.
"""
function site_directory(site_ID)
    root = ClimaLand.Artifacts.fluxnet2015_data_path()
    return joinpath(
        root,
        only(filter(startswith("FLX_$(site_ID)_"), readdir(root))),
    )
end

"""
    monthly_observations(site_ID)

Return a `Dict` from FLUXNET2015 column name to a `Dict` from `(year, month)`
to the value of the site's FULLSET monthly file, the source of ILAMB's
FLUXNET2015 benchmark. Gaps (-9999) are left out.
"""
function monthly_observations(site_ID)
    site_dir = site_directory(site_ID)
    parts = split(basename(site_dir), "FULLSET")
    path = joinpath(site_dir, string(parts[1], "FULLSET_MM", parts[2], ".csv"))
    isfile(path) || error("$site_ID has no monthly FULLSET file $path")
    data, header = DelimitedFiles.readdlm(path, ','; header = true)
    header = vec(header)
    stamps = Int.(data[:, findfirst(==("TIMESTAMP"), header)])
    obs = Dict{String, Dict{Tuple{Int, Int}, Float64}}()
    for (_, column) in VARIABLES
        j = findfirst(==(column), header)
        isnothing(j) && continue
        obs[column] = Dict(
            (s ÷ 100, s % 100) => Float64(v) for
            (s, v) in zip(stamps, data[:, j]) if v != -9999
        )
    end
    return obs
end

"""
    site_lai(lat, long, start_date, stop_date)

Return a `TimeVaryingInput` of the MODIS LAI of the 1° cell nearest to the site,
in seconds from `start_date`, covering `start_date` to `stop_date`. Years with a
file in the `modis_lai` artifact (2000-2020) use it; the other years use the
MODIS climatology.
"""
function site_lai(lat, long, start_date, stop_date)
    modis_dir = ClimaLand.Artifacts.modis_lai_forcing_data_path()
    climatology = ClimaLand.Artifacts.modis_lai_climatology_data_path()
    times = DateTime[]
    values = Float64[]
    for y in (year(start_date) - 1):(year(stop_date) + 1)
        path = joinpath(modis_dir, "Yuan_et_al_$(y)_1x1.nc")
        observed = isfile(path)
        NCDatasets.NCDataset(observed ? path : climatology) do ds
            i = argmin(abs.(mod.(ds["lon"][:] .- long .+ 180, 360) .- 180))
            j = argmin(abs.(ds["lat"][:] .- lat))
            # The climatology is stored as the year 2000
            shift = observed ? Year(0) : Year(y - 2000)
            append!(times, ds["time"][:] .+ shift)
            append!(values, Float64.(ds["lai"][i, j, :]))
        end
    end
    seconds = [Second(t - start_date).value for t in times]
    return TimeVaryingInput(Float64.(seconds), values)
end

"""
    run_site(site_ID)

Run the default `LandModel` at `site_ID` over the span in which every forcing
variable is observed, and return the diagnostics writer, the start of the scored
period (UTC: the first local month boundary at least `SPINUP` after the start)
and the site's data timestep and UTC offset.
"""
function run_site(site_ID)
    toml_dict = LP.create_toml_dict(FT)
    site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
    (start_date, stop_date) = FluxnetSimulations.get_data_dates(
        site_ID,
        time_offset;
        required_columns = FluxnetSimulations.FLUXNET_FORCING_COLUMNS,
    )
    data_dt = FluxnetSimulations.get_data_dt(site_ID)
    # Dates are the midpoints of the averaging windows of the data
    half_dt = Second(data_dt ÷ 2)
    utc_to_local = Minute(round(Int, 60 * time_offset))
    scored_start =
        ceil(start_date - half_dt + SPINUP + utc_to_local, Month(1)) -
        utc_to_local
    scored_start + MIN_SCORED - half_dt <= stop_date || error(
        "$site_ID: the record ends at $stop_date, less than $MIN_SCORED after the spin-up ends at $scored_start",
    )

    Δt = 900.0
    domain = Column(;
        zlim = (FT(-15), FT(0)),
        nelements = 15,
        dz_tuple = (FT(3), FT(0.05)),
        longlat = (long, lat),
    )
    forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict,
        FT,
    )
    LAI = site_lai(lat, long, start_date, stop_date)
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    )
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )
    writer = ClimaDiagnostics.Writers.DictWriter()
    diags = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = writer,
        output_vars = collect(String.(keys(VARIABLES))),
        reduction_period = data_dt == 3600 ? :hourly : :halfhourly,
    )
    simulation = LandSimulation(
        start_date,
        stop_date,
        Δt,
        land;
        set_ic!,
        updateat = Second(data_dt),
        diagnostics = diags,
    )
    solve!(simulation)
    return (; writer, scored_start, data_dt, time_offset)
end

"""
    monthly_model(writer, short_name, scored_start, data_dt, time_offset)

Return a `Dict` from `(year, month)` to the mean of the diagnostic over that
calendar month in local standard time, the time of the FLUXNET2015 monthly
files, for the complete months from `scored_start` on.
"""
function monthly_model(writer, short_name, scored_start, data_dt, time_offset)
    name =
        only(filter(startswith("$(short_name)_"), collect(keys(writer.dict))))
    times, values = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, name)
    any(isnan, values) && error("$short_name has NaNs")
    utc_to_local = Minute(round(Int, 60 * time_offset))
    # Diagnostic times are the end of each averaging window
    midpoints = ClimaUtilities.TimeManager.date.(times) .- Second(data_dt ÷ 2)
    sums = Dict{Tuple{Int, Int}, Tuple{Float64, Int}}()
    for (t, v) in zip(midpoints, values)
        t < scored_start && continue
        local_t = t + utc_to_local
        key = (year(local_t), month(local_t))
        s, n = get(sums, key, (0.0, 0))
        sums[key] = (s + v, n + 1)
    end
    complete(key) = daysinmonth(key...) * 86400 ÷ data_dt
    return Dict(key => s / n for (key, (s, n)) in sums if n == complete(key))
end

function score_site(site_ID)
    obs = monthly_observations(site_ID)
    (; writer, scored_start, data_dt, time_offset) = run_site(site_ID)
    open(joinpath(OUTDIR, "$site_ID.csv"), "w") do io
        println(io, "variable,year,month,model,obs")
        for (short_name, (ilamb_name, column)) in pairs(VARIABLES)
            haskey(obs, column) || continue
            model = monthly_model(
                writer,
                short_name,
                scored_start,
                data_dt,
                time_offset,
            )
            for key in sort(collect(keys(model)))
                o = get(obs[column], key, NaN)
                println(io, join((ilamb_name, key..., model[key], o), ","))
            end
        end
    end
end

function default_sites()
    cohort, _ = DelimitedFiles.readdlm(
        joinpath(@__DIR__, "ilamb_land_hist_site_rmse.csv"),
        ',',
        String;
        header = true,
    )
    # These four resolve to the single-year `fluxnet_sites` artifact
    curated = ("US-MOz", "US-Var", "US-NR1", "US-Ha1")
    return filter(
        s -> s in cohort[:, 1] && !(s in curated),
        list_fluxnet_sites(),
    )
end

"""
    job_sites(sites, job, n_jobs)

Return the sites of job `job` (0-based) of `n_jobs`. Sites are taken longest
record first, each by the job with the fewest record years so far.
"""
function job_sites(sites, job, n_jobs)
    n_jobs == 1 && return sites
    years = Dict(map(sites) do site
        m = match(r"FULLSET_(\d{4})-(\d{4})", site_directory(site))
        site => parse(Int, m[2]) - parse(Int, m[1]) + 1
    end)
    load = zeros(Int, n_jobs)
    assigned = [String[] for _ in 1:n_jobs]
    for site in sort(sites; by = s -> years[s], rev = true)
        k = argmin(load)
        push!(assigned[k], site)
        load[k] += years[site]
    end
    return assigned[job + 1]
end

function main(sites)
    job = parse(Int, get(ENV, "BUILDKITE_PARALLEL_JOB", "0"))
    n_jobs = parse(Int, get(ENV, "BUILDKITE_PARALLEL_JOB_COUNT", "1"))
    sites = job_sites(sites, job, n_jobs)
    mkpath(OUTDIR)
    @info "Running $(length(sites)) sites" sites
    for site_ID in sites
        @info "Site $site_ID"
        try
            t = @elapsed score_site(site_ID)
            @info "Site $site_ID done in $(round(t; digits = 1)) s"
        catch e
            msg = sprint(showerror, e)
            @warn "Site $site_ID failed" exception = msg
            write(joinpath(OUTDIR, "$site_ID.failed"), msg)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? default_sites() : ARGS)
end

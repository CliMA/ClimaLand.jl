# Phase 2 gate: 21-site padded forcing is faithful, ordered, and integrable.
#
#   1. Padding correctness on the full union axis: native samples are
#      untouched; recycled samples equal the mapped interior-year sample
#      (rule recomputed independently here); everything is finite.
#   2. A 21-column ColumnEnsemble built from a cropped window reproduces every
#      site's raw driver values at sample times, in column order.
#   3. The default LandModel integrates 5 days on all 21 columns.
#
# Run:  julia --project=.buildkite experiments/callmip_phase1b/verify_21site_forcing.jl
# GPU:  CLIMACOMMS_DEVICE=CUDA julia --project=.buildkite ...

import ClimaComms
ClimaComms.@import_required_backends

using Test
using Dates
import ClimaCore
import Random

include(joinpath(@__DIR__, "forcing_phase1b.jl"))

FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
toml_dict = LP.create_toml_dict(FT)

@info "Loading all 21 calibration sites (full union axis)"
raw = Dict(id => load_site(id) for id in CALIBRATION_SITE_IDS)
padded = load_calibration_sites()
axis = padded[1].utc_dates
@info "Union axis" first(axis) last(axis) length(axis)

@testset "Padding correctness (full union axis)" begin
    rng = Random.MersenneTwister(2026)
    for (j, id) in enumerate(CALIBRATION_SITE_IDS)
        r, p = raw[id], padded[j]
        @test p.utc_dates == axis
        @test count(p.native) == length(r.utc_dates)
        @test all(isfinite, p.Tair) && all(isfinite, p.SWdown)
        @test isnothing(p.LAI) || all(isfinite, p.LAI)
        ridx = Dict(t => i for (i, t) in enumerate(r.utc_dates))
        yr_lo = Dates.year(first(r.utc_dates)) + 1
        yr_hi = Dates.year(last(r.utc_dates)) - 1
        native_ks = findall(p.native)
        recycled_ks = findall(.!p.native)
        for k in rand(rng, native_ks, 12)
            i = ridx[axis[k]]
            @test p.Tair[k] == r.Tair[i] && p.SWdown[k] == r.SWdown[i]
            @test isnothing(p.LAI) || p.LAI[k] == r.LAI[i]
        end
        isempty(recycled_ks) && continue
        for k in rand(rng, recycled_ks, 12)
            t = axis[k]
            ysrc = yr_lo + mod(Dates.year(t) - yr_lo, yr_hi - yr_lo + 1)
            d = Dates.day(t)
            Dates.month(t) == 2 &&
                d == 29 &&
                !Dates.isleapyear(ysrc) &&
                (d = 28)
            i = ridx[DateTime(
                ysrc,
                Dates.month(t),
                d,
                Dates.hour(t),
                Dates.minute(t),
            )]
            @test p.Tair[k] == r.Tair[i] && p.SWdown[k] == r.SWdown[i]
        end
    end
end

start_date = DateTime(2010, 7, 1)
stop_date = start_date + Day(5)
Δt = 450.0
window = (start_date - Day(1), stop_date + Day(1))
sites_w = load_calibration_sites(; window)

@testset "21-column ensemble drivers in column order" begin
    domain = ColumnEnsemble(;
        zlim = FT.((-2, 0)),
        nelements = 10,
        longlat = [FT.((s.long, s.lat)) for s in sites_w],
    )
    surface_space = domain.space.surface
    forcing = prescribed_forcing_callmip(
        sites_w,
        surface_space,
        start_date,
        toml_dict,
        FT,
    )
    utc_dates, aligned = align_sites(sites_w)
    seconds = Float64.(Dates.value.(utc_dates .- start_date)) ./ 1000
    dest = ClimaCore.Fields.zeros(surface_space)
    n = length(utc_dates)
    for k in unique(round.(Int, range(1, n, length = 7)))
        for (tvi, raw_of) in (
            (forcing.atmos.T, s -> s.Tair),
            (forcing.radiation.SW_d, s -> s.SWdown),
            (forcing.atmos.u, s -> s.Wind),
        )
            evaluate!(dest, tvi, seconds[k])
            vals = Array(ClimaCore.Fields.field2array(dest))
            expected = [raw_of(s)[k] for s in aligned]
            @test all(isapprox.(vals, expected; rtol = 1e-12))
        end
    end
    heights = Array(ClimaCore.Fields.field2array(forcing.atmos.h))
    expected_h = [site_info(id).zref for id in CALIBRATION_SITE_IDS]
    @test all(isapprox.(heights, expected_h; atol = 1e-3))
end

@info "Forward run: 21 columns, 5 days, default LandModel" Δt
simulation =
    build_callmip_simulation(sites_w, start_date, stop_date, Δt, toml_dict, FT)
stats = @timed solve!(simulation)
sypd = (5 / 365.25) / (stats.time / 86400)
@info "21-column integration done" wall_seconds = round(stats.time; digits = 1) approx_SYPD =
    round(sypd; digits = 1)

@testset "21-column forward run is finite per column" begin
    Y = simulation._integrator.u
    p = simulation._integrator.p
    for (i, id) in enumerate(CALIBRATION_SITE_IDS)
        for f in (Y.canopy.energy.T, Y.soil.ϑ_l, p.soil.T, p.drivers.T)
            a = parent(f)
            col = vec(Array(selectdim(a, ndims(a), i)))
            @test all(isfinite, col)
        end
    end
    T_canopy = [
        vec(
            Array(
                selectdim(
                    parent(Y.canopy.energy.T),
                    ndims(parent(Y.canopy.energy.T)),
                    i,
                ),
            ),
        )[1] for i in 1:21
    ]
    @info "Final canopy T per site" extrema(T_canopy)
    @test all(t -> 240 < t < 330, T_canopy)
end
@info "verify_21site_forcing.jl: all testsets finished"

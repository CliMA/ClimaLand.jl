# Compares two saved states of run_trendy.jl on CPU, layer by layer: the soil O2
# (area-weighted and SOC-weighted means, and the land share where it fell by more than
# 10% and 50%), and the change of SOC (Pg C) by layer and latitude band. Used to test
# whether a faster respiration lowers O2 and so Sm below the curves of an earlier cycle.
#
# With `pre` (the state the SOC solve started from), the column drift is also given without
# the solve's no-equilibrium layers (where A kept the SOC of `pre`), as G1 leaves them out.
#
# Usage (CPU job): julia --project=.buildkite o2_soc_change.jl <state_a.hdf5> <state_b.hdf5> [years] [pre.hdf5]
import ClimaComms
import ClimaCore
using ClimaLand
using Printf

const FT = Float64
include(joinpath(@__DIR__, "carbon_budget.jl"))
include(joinpath(@__DIR__, "state_io.jl"))
const BANDS = (("global", -60, 90), ("north_50N", 50, 90), ("30N_50N", 30, 50), ("tropics_30S_30N", -30, 30), ("south_60S_30S", -60, -30))

a_path, b_path = ARGS[1], ARGS[2]
years = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1.0
pre_path = length(ARGS) >= 4 ? ARGS[4] : nothing
context = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
domain = ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
nodes, area = land_nodes_and_areas(domain.space.surface)
lat = vec(Array(parent(ClimaCore.Fields.coordinate_field(domain.space.surface).lat)))[nodes]
z = Array(parent(ClimaCore.Fields.coordinate_field(domain.space.subsurface).z))[:, 1, 1, 1, nodes]
dz = Array(parent(ClimaCore.Fields.Δz_field(domain.space.subsurface)))[:, 1, 1, 1, nodes]
nz = size(z, 1)
field(path, name) = read_dataset(path, "Y", "soilco2", name)[:, 1, 1, 1, nodes]
o2a, o2b = field(a_path, "O2"), field(b_path, "O2")
soca, socb = field(a_path, "SOC"), field(b_path, "SOC")
ok = [all(isfinite, view(o2a, :, i)) && all(isfinite, view(o2b, :, i)) && all(isfinite, view(soca, :, i)) && all(isfinite, view(socb, :, i)) for i in eachindex(nodes)]
println("A = ", a_path, "\nB = ", b_path, "\nland nodes ", length(nodes), ", non-finite left out ", count(!, ok))

trop = abs.(lat) .< 30
println("\nB/A(SOC_B sums): Σ O2_B·w / Σ O2_A·w with w = area × SOC of B, over all land and the tropics")
println("\nlayer  depth_m   O2_A    O2_B   B/A(area) B/A(SOC-wtd) share<0.9 share<0.5  dSOC_PgC/yr  B/A(SOC_B sums) tropics")
for l in 1:nz
    w = area .* ok
    ws = w .* soca[l, :]
    mA, mB = sum(o2a[l, :] .* w) / sum(w), sum(o2b[l, :] .* w) / sum(w)
    rs = sum(o2b[l, :] ./ max.(o2a[l, :], 1e-12) .* ws) / sum(ws)
    r = o2b[l, :] ./ max.(o2a[l, :], 1e-12)
    s9 = sum(w .* (r .< 0.9)) / sum(w)
    s5 = sum(w .* (r .< 0.5)) / sum(w)
    dsoc = sum((socb[l, :] .- soca[l, :]) .* dz[l, :] .* w) * 1e-12 / years
    wb = w .* socb[l, :]
    rsum(m) = sum((o2b[l, :] .* wb)[m]) / sum((o2a[l, :] .* wb)[m])
    @printf("%3d  %7.2f  %.2e %.2e  %7.3f   %7.3f      %6.3f   %6.3f   %+8.3f      %6.3f   %6.3f\n", l, -z[l, 1], mA, mB, mB / mA, rs, s9, s5, dsoc, rsum(trues(length(w))), rsum(trop))
end

println("\nband             dSOC_0-1m dSOC_1m+ (Pg C/yr)  area_Mkm2  column drift share<2 and <10 g/m2/yr")
above(D) = max.(min.(-z .+ dz ./ 2, D) .- max.(-z .- dz ./ 2, 0.0), 0.0)
w1 = above(1.0)
dcol = vec(sum((socb .- soca) .* dz; dims = 1)) ./ years      # kg C m-2 yr-1
dtop = vec(sum((socb .- soca) .* w1; dims = 1)) ./ years
function band_table(dcol, dtop)
    for (name, lo, hi) in BANDS
        m = ok .& (lo .<= lat .< hi)
        A = sum(area[m])
        top = sum(dtop[m] .* area[m]) * 1e-12
        deep = sum((dcol[m] .- dtop[m]) .* area[m]) * 1e-12
        s2 = sum(area[m] .* (abs.(dcol[m]) .< 2e-3)) / A
        s10 = sum(area[m] .* (abs.(dcol[m]) .< 10e-3)) / A
        @printf("%-16s %+8.3f  %+8.3f            %7.2f    %6.3f  %6.3f\n", name, top, deep, A * 1e-12, s2, s10)
    end
end
band_table(dcol, dtop)

if !isnothing(pre_path)
    eq = soca .!= field(pre_path, "SOC")
    d = (socb .- soca) .* eq
    println("\nwithout the solve's no-equilibrium layers (", count(.!eq), " of ", length(eq), " layer-nodes; pre = ", pre_path, ")")
    band_table(vec(sum(d .* dz; dims = 1)) ./ years, vec(sum(d .* w1; dims = 1)) ./ years)
end

# Phase 1a, the offline choice of the DAMM parameters, on CPU, from the sums of a spin-up
# cycle run with alternative curve sets (SOC_ALTS; spinup_tools.jl).
#
# Litter: the equilibrium litter of the cycle's pool sums (Canopy.equilibrium_carbon_pools
# with the parameters in use), on the above-ground and root profiles. The target is
# SoilGrids as the model reads it (`set_soilco2_SOC_from_soilgrids!` on the model's
# layers and land nodes), with the same depth weights as the solve (the part of each layer
# above 1 m). For every curve set (the model's, and each alternative (ea, kmo2, doa)) and
# V_ref multiplier, the SOC scale s (kM multiplier) that matches the global top metre is
# found by bisection; then the capped litter (top metre and below), the SOC by band
# against SoilGrids, and the turnover. In the solve, layers with no equilibrium keep the
# SOC of the state given (not a cap). For the best candidates, the steady state of the
# leaderboard, SOC × litter/Sm at the state's SOC, is compared with the exact solve.
#
# Usage (CPU job):
#   julia --project=.buildkite damm_choice.jl <accum.hdf5> <state.hdf5> <targets.csv> <out_dir> \
#       [<set> <vref_mult> <s>]
# With the last three arguments (set 0 = the model's curves, j = alternative j), writes
# <out_dir>/state_damm.hdf5: <state.hdf5> with the SOC of that solve, and
# <out_dir>/respiration_override.toml with the matching parameters.
import ClimaComms
import ClimaCore
using ClimaLand
import ClimaLand.Canopy
import ClimaLand.Parameters as LP
using Printf

const FT = Float64
include(joinpath(@__DIR__, "carbon_budget.jl"))
include(joinpath(@__DIR__, "state_io.jl"))
include(joinpath(@__DIR__, "spinup_tools.jl"))
const YEAR = 365.25 * 86400
const BANDS = (("global", -60, 90), ("north_50N", 50, 90), ("30N_50N", 30, 50), ("tropics_30S_30N", -30, 30), ("south_60S_30S", -60, -30))
key(k) = "c" * lpad(k, 2, '0')

accum_path, state_file, targets_path, out_dir = ARGS[1:4]
mkpath(out_dir)
context = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
domain = ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
toml_dict = LP.create_toml_dict(FT)
surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
biomass = Canopy.PrognosticCarbonModel{FT}(Canopy.ZhouOptimalLAIModel{FT}(surface_domain, toml_dict), toml_dict)
params = biomass.parameters

nodes, area = land_nodes_and_areas(domain.space.surface)
lat = vec(Array(parent(ClimaCore.Fields.coordinate_field(domain.space.surface).lat)))[nodes]
nl = length(nodes)
z = Array(parent(ClimaCore.Fields.coordinate_field(domain.space.subsurface).z))[:, 1, 1, 1, nodes]
dz = Array(parent(ClimaCore.Fields.Δz_field(domain.space.subsurface)))[:, 1, 1, 1, nodes]
nz = size(z, 1)
rooting_depth = vec(Array(parent(biomass.rooting_depth)))[nodes]
# thickness (m) of each layer above the depth D
above(D) = max.(min.(-z .+ dz ./ 2, D) .- max.(-z .- dz ./ 2, 0.0), 0.0)
w1, w2 = above(1.0), above(2.0)
top1 = w1 .> 0
sg_field = ClimaCore.Fields.zeros(domain.space.subsurface)
ClimaLand.Simulations.set_soilco2_SOC_from_soilgrids!(sg_field)
soilgrids = Array(parent(sg_field))[:, 1, 1, 1, nodes]

T_pool = read_root_attribute(accum_path, "pool_seconds")
T_soc = read_root_attribute(accum_path, "soc_seconds")
alts = parse_alts(read_root_attribute(accum_path, "alts"))
surf(name) = vec(read_dataset(accum_path, "acc", "pools", name))[nodes] ./ T_pool
GPP, Rd, f_T, MAT, tree_share = surf("GPP"), surf("Rd"), surf("f_T"), surf("T_air"), surf("tree_share")
sub(chain...) = read_dataset(accum_path, "acc", chain...)[:, 1, 1, 1, nodes] ./ T_soc
litter_run = sub("litter")
curve_set(j) = cat((j == 0 ? sub("curves", key(k)) : sub("alts", "a$j", key(k)) for k in eachindex(SOC_GRID))...; dims = 3)
soc_state = read_dataset(state_file, "Y", "soilco2", "SOC")[:, 1, 1, 1, nodes]

finite_node = [isfinite(GPP[i] + Rd[i] + f_T[i] + MAT[i] + tree_share[i]) && all(isfinite, view(litter_run, :, i)) && all(isfinite, view(soc_state, :, i)) for i in 1:nl]
println("land nodes ", nl, "; left out (non-finite) ", count(!, finite_node), "; alternative sets: ", alts)
band_masks = [(name, finite_node .& (lo .<= lat .< hi)) for (name, lo, hi) in BANDS]
integral(x, mask) = sum(x[i] * area[i] for i in 1:nl if mask[i]; init = 0.0)
column(x) = vec(sum(x .* dz; dims = 1))

targets = Dict{Tuple{String, String}, Float64}()
for line in readlines(targets_path)[2:end]
    f = split(line, ",")
    targets[(f[1], f[2])] = parse(Float64, f[3])
end
model_grid(name, mask, w) = integral(vec(sum(soilgrids .* w; dims = 1)), mask) * 1e-12
for (name, mask) in band_masks
    @printf("SoilGrids 0-1 m, %s: on its own grid %.0f Pg C; on the model's layers and land %.0f Pg C (0-2 m %.0f)\n",
        name, targets[("SOC_0_1m_soilgrids", name)], model_grid(name, mask, w1), model_grid(name, mask, w2))
    targets[("SOC_0_1m_soilgrids", name)] = model_grid(name, mask, w1)
    targets[("SOC_0_2m_soilgrids", name)] = model_grid(name, mask, w2)
end

# equilibrium litter of the pool sums
L_eq = zeros(nz, nl)
for i in 1:nl
    finite_node[i] || continue
    pools = Canopy.equilibrium_carbon_pools(params, GPP[i], Rd[i], f_T[i], MAT[i], tree_share[i])
    τ_stem = params.τ_stem0 * Canopy.tau_stem_scale(MAT[i], params.T_ref_τ_stem, params.q_τ_stem)
    above = pools.C_leaf / params.τ_leaf + pools.C_stem / τ_stem
    root = pools.C_root / params.τ_root
    s = [exp(z[l, i] / params.soil_litter_depth) for l in 1:nz]
    r = [exp(z[l, i] / rooting_depth[i]) for l in 1:nz]
    s ./= sum(s .* dz[:, i])
    r ./= sum(r .* dz[:, i])
    @. L_eq[:, i] = above * s + root * r
end
all_land = band_masks[1][2]
@printf("litter: cycle mean %.2f Pg C/yr, equilibrium of the pool sums %.2f Pg C/yr\n",
    integral(column(litter_run), all_land) * YEAR * 1e-12, integral(column(L_eq), all_land) * YEAR * 1e-12)

# Between candidate SOC values, log R is interpolated linearly in log c: exact where R is a
# power of c (as in the linear regime), where linear R in log c overestimates R by up to 1%
# (curve_check.jl on C4's year 1870: +0.3% globally against the run's Sm; log-log −0.14%).
log_root(L, Rk, Rk1, lo, hi) =
    Rk > 0 && Rk1 > Rk ? exp(lo + log(L / Rk) / log(Rk1 / Rk) * (hi - lo)) :
    exp(lo + (L - Rk) / (Rk1 - Rk) * (hi - lo))
log_interp(f, a, b) = a > 0 && b > 0 ? exp(log(a) + f * log(b / a)) : a + f * (b - a)

"""SOC (kg C m⁻³) and the no-equilibrium mask for curves `C` × `rate`, SOC scale `s`."""
function solve(C, L, rate, s)
    soc = copy(soc_state)
    capped = falses(nz, nl)
    logc = log.(SOC_GRID)
    K = length(SOC_GRID)
    for i in 1:nl, l in 1:nz
        finite_node[i] || continue
        Li = L[l, i]
        if !(Li > 0)
            soc[l, i] = 0.0
            continue
        end
        R1, RK = rate * C[l, i, 1], rate * C[l, i, K]
        if Li >= RK
            capped[l, i] = true
            continue
        end
        c = if Li < R1
            SOC_GRID[1] * Li / R1
        else
            k = 1
            while k < K - 1 && rate * C[l, i, k + 1] <= Li
                k += 1
            end
            Rk, Rk1 = rate * C[l, i, k], rate * C[l, i, k + 1]
            log_root(Li, Rk, Rk1, logc[k], logc[k + 1])
        end
        soc[l, i] = c * s
    end
    return soc, capped
end

top_metre(soc, mask) = integral(vec(sum(soc .* w1; dims = 1)), mask) * 1e-12

function report(C, L, rate, s)
    soc, capped = solve(C, L, rate, s)
    lit = integral(column(L), all_land) * YEAR * 1e-12
    cap_top = integral(column(L .* capped .* top1), all_land) * YEAR * 1e-12
    cap_deep = integral(column(L .* capped .* .!top1), all_land) * YEAR * 1e-12
    bands = [(name, top_metre(soc, m), targets[("SOC_0_1m_soilgrids", name)]) for (name, m) in band_masks]
    misfit = sum(abs(log(b[2] / b[3])) for b in bands[2:end])
    lit1 = integral(column(L .* top1), all_land) * YEAR * 1e-12
    return (; soc, capped, lit, cap_top, cap_deep, bands, misfit,
        soc2 = integral(vec(sum(soc .* w2; dims = 1)), all_land) * 1e-12,
        turnover = bands[1][2] / lit1)
end

target_1m = targets[("SOC_0_1m_soilgrids", "global")]
rows = String["set,ea,kmo2,doa,vref_mult,s,litter,cap_top1m,cap_below1m,cap_share,SOC_0_1m,SOC_0_2m,turnover_1m_yr,misfit,north_50N,30N_50N,tropics,south_60S_30S"]
println("\nset (ea, kmo2, doa) | V× | s | capped top/below (Pg C/yr) share | SOC 0-1 m, 0-2 m | τ 1 m | band ratios model/SoilGrids (N50, 30-50N, tropics, S30) | misfit")
for j in 0:length(alts)
    C = curve_set(j)
    a = j == 0 ? (; ea = NaN, kmo2 = 1.0, doa = 1.0) : alts[j]
    for rate in (1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0)
        lo, hi = 0.01, 100.0
        for _ in 1:30
            mid = sqrt(lo * hi)
            soc, _ = solve(C, L_eq, rate, mid)
            top_metre(soc, all_land) < target_1m ? (lo = mid) : (hi = mid)
        end
        s = sqrt(lo * hi)
        r = report(C, L_eq, rate, s)
        ratios = [b[2] / b[3] for b in r.bands[2:end]]
        share = (r.cap_top + r.cap_deep) / r.lit
        push!(rows, join((j, a.ea, a.kmo2, a.doa, rate, s, r.lit, r.cap_top, r.cap_deep, share, r.bands[1][2], r.soc2, r.turnover, r.misfit, ratios...), ","))
        @printf("%d (%s, %.2f, %.1f) | %4.0f | %6.2f | %5.2f/%5.2f %.3f | %5.0f %5.0f | %5.1f | %.2f %.2f %.2f %.2f | %.2f\n",
            j, isnan(a.ea) ? "model" : string(a.ea), a.kmo2, a.doa, rate, s, r.cap_top, r.cap_deep, share, r.bands[1][2], r.soc2, r.turnover, ratios..., r.misfit)
    end
end
write(joinpath(out_dir, "damm_choice.csv"), join(rows, "\n") * "\n")

"""Respiration (kg C m⁻³ s⁻¹) at SOC `c` from curves `C` × `rate` with SOC scale `s`."""
function respiration_at(C, l, i, c, rate, s)
    x = c / s
    K = length(SOC_GRID)
    x <= SOC_GRID[1] && return rate * C[l, i, 1] * x / SOC_GRID[1]
    x >= SOC_GRID[K] && return rate * C[l, i, K]
    k = searchsortedlast(SOC_GRID, x)
    f = (log(x) - log(SOC_GRID[k])) / (log(SOC_GRID[k + 1]) - log(SOC_GRID[k]))
    return rate * log_interp(f, C[l, i, k], C[l, i, k + 1])
end

# the leaderboard's steady state, SOC × litter/Sm (kept where litter/Sm > 100), from the
# state's SOC and from SoilGrids, against the exact solve
println("\nleaderboard estimate against the exact solve (0-1 m, Pg C): global and by band")
best = sort(rows[2:end]; by = r -> parse(Float64, split(r, ",")[14]))
for r in vcat(filter(r -> split(r, ",")[1] == "0" && parse(Float64, split(r, ",")[5]) in (12.0, 20.0), rows[2:end]), best[1:3])
    f = split(r, ",")
    j, rate, s = parse(Int, f[1]), parse(Float64, f[5]), parse(Float64, f[6])
    C = curve_set(j)
    exact = report(C, L_eq, rate, s)
    for (label, start) in (("state", soc_state), ("SoilGrids", soilgrids))
        est = copy(start)
        for i in 1:nl, l in 1:nz
            finite_node[i] || continue
            R = respiration_at(C, l, i, start[l, i], rate, s)
            ratio = L_eq[l, i] / R
            est[l, i] = ratio <= 100 ? start[l, i] * ratio : start[l, i]
        end
        @printf("  set %d V× %g s %.2f from %-9s: exact %s | estimate %s\n", j, rate, s, label,
            join([@sprintf("%.0f", top_metre(exact.soc, m)) for (_, m) in band_masks], " "),
            join([@sprintf("%.0f", top_metre(est, m)) for (_, m) in band_masks], " "))
    end
end

if length(ARGS) >= 7
    j, rate, s = parse(Int, ARGS[5]), parse(Float64, ARGS[6]), parse(Float64, ARGS[7])
    a = j == 0 ? (; ea = NaN, kmo2 = 1.0, doa = 1.0) : alts[j]
    r = report(curve_set(j), L_eq, rate, s)
    full = read_dataset(state_file, "Y", "soilco2", "SOC")
    full[:, 1, 1, 1, nodes] .= r.soc
    out_state = joinpath(out_dir, "state_damm.hdf5")
    cp(state_file, out_state; force = true)
    HDF5.h5open(out_state, "r+") do f
        f["fields/Y/soilco2/SOC"][:, :, :, :, :] = full
    end
    damm = Dict(
        "soilCO2_reference_rate" => toml_dict["soilCO2_reference_rate"] * rate,
        "michaelis_constant" => toml_dict["michaelis_constant"] * s,
        "O2_michaelis_constant" => toml_dict["O2_michaelis_constant"] * a.kmo2,
        "oxygen_diffusion_coefficient" => toml_dict["oxygen_diffusion_coefficient"] * a.doa,
        "soilCO2_activation_energy" => isnan(a.ea) ? toml_dict["soilCO2_activation_energy"] : a.ea,
    )
    open(joinpath(out_dir, "respiration_override.toml"), "w") do io
        println(io, "# DAMM parameters for the land-sink study (PLAN decision 7), from damm_choice.jl")
        println(io, "# curve set $j $(a), V_ref ×$rate, kM ×$s; the carbon-pool (Ra) parameters keep their defaults")
        for (name, value) in sort(collect(damm))
            println(io, "\n[$name]\nvalue = $value\ntype = \"float\"")
        end
    end
    @printf("wrote %s: SOC 0-1 m %.0f, 0-2 m %.0f, full %.0f Pg C; capped litter %.2f + %.2f Pg C/yr\n",
        out_state, r.bands[1][2], r.soc2, integral(column(r.soc), all_land) * 1e-12, r.cap_top, r.cap_deep)
end

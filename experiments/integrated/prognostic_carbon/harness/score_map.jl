## Global skill of an equilibrium biomass map, and the payload for the figures.
##
## `global_equilibrium.jl` prints bias and spatial correlation as it runs. This
## adds RMSE, reports every product in one table, and writes the binned map
## payload the artifact plots, so that every published number comes from a
## command rather than from a previous figure adjusted by eye.
##
## RMSE is reported alongside bias because they answer different questions and
## the products disagree enough that neither alone is sufficient: bias is a
## statement about the global total, RMSE about whether individual cells are
## right, and a map can be unbiased and wrong everywhere.
##
## Usage:
##   julia --project=.buildkite score_map.jl <equilibrium_carbon.nc> [maps.json]

include(joinpath(@__DIR__, "global_equilibrium.jl"))

"Area-weighted skill of `M` against product `name`, on the model grid."
function score(name, path, lons, lats, M)
    olons, olats, O = obs_grid(name, path)
    m, o = FT[], FT[]
    for i in eachindex(lons), j in eachindex(lats)
        valid(M[i, j]) || continue
        x = nearest(olons, olats, O, lons[i], lats[j])
        valid(x) || continue
        push!(m, FT(M[i, j]))
        push!(o, FT(x))
    end
    n = length(m)
    n == 0 && return nothing
    mb, ob = sum(m) / n, sum(o) / n
    sm = sqrt(sum((m .- mb) .^ 2) / n)
    so = sqrt(sum((o .- ob) .^ 2) / n)
    r = sum((m .- mb) .* (o .- ob)) / n / max(sm * so, eps(FT))
    rmse = sqrt(sum((m .- o) .^ 2) / n)
    return (; n, model = mb, obs = ob, bias = mb - ob, rmse, r)
end

"""
    block_means(lons, lats, A, step)

Coarsen to `step`-degree blocks for plotting. The figure does not need 1-degree
detail and a full-resolution payload is megabytes of JSON.
"""
function block_means(lons, lats, A, step)
    nb_lon = ceil(Int, length(lons) / step)
    nb_lat = ceil(Int, length(lats) / step)
    out = fill(NaN, nb_lon, nb_lat)
    for bi in 1:nb_lon, bj in 1:nb_lat
        s, n = 0.0, 0
        for i in ((bi - 1) * step + 1):min(bi * step, length(lons)),
            j in ((bj - 1) * step + 1):min(bj * step, length(lats))

            valid(A[i, j]) && (s += A[i, j]; n += 1)
        end
        n > 0 && (out[bi, bj] = s / n)
    end
    return out
end

function main(ncpath, jsonpath = nothing)
    ds = NCDatasets.NCDataset(ncpath)
    lons = Array{FT}(Array(ds["lon"][:]))
    lats = Array{FT}(Array(ds["lat"][:]))
    # C_stem, not cVeg. Most of these products report woody or above-ground
    # biomass, and more importantly every bias and correlation reported for this
    # model so far is against C_stem — scoring cVeg here would produce a table
    # that silently disagrees with all of them.
    M = Array{Union{Missing, FT}}(Array(ds["C_stem"][:, :]))
    close(ds)

    println(
        rpad("product", 14),
        rpad("cells", 8),
        rpad("model", 9),
        rpad("obs", 9),
        rpad("bias", 9),
        rpad("RMSE", 9),
        "r",
    )
    rows = []
    for (name, path) in PRODUCTS
        isfile(path) || continue
        s = score(name, path, lons, lats, M)
        s === nothing && continue
        push!(rows, (name, s))
        println(
            rpad(name, 14),
            rpad(s.n, 8),
            rpad(round(s.model, digits = 2), 9),
            rpad(round(s.obs, digits = 2), 9),
            rpad(round(s.bias, digits = 2), 9),
            rpad(round(s.rmse, digits = 2), 9),
            round(s.r, digits = 3),
        )
    end
    k = length(rows)
    println(
        "\nacross products: mean bias ",
        round(sum(r[2].bias for r in rows) / k, digits = 2),
        ", mean RMSE ",
        round(sum(r[2].rmse for r in rows) / k, digits = 2),
        ", mean r ",
        round(sum(r[2].r for r in rows) / k, digits = 3),
        " (kg C m-2)",
    )

    jsonpath === nothing && return
    # Reference product for the observed and bias panels: XuSaatchi reports
    # carbon directly and covers the most cells, so it needs no dry-matter
    # conversion and no extrapolation.
    olons, olats, O = obs_grid("XuSaatchi", PRODUCTS[1][2])
    Og = fill(FT(NaN), length(lons), length(lats))
    for i in eachindex(lons), j in eachindex(lats)
        x = nearest(olons, olats, O, lons[i], lats[j])
        valid(x) && valid(M[i, j]) && (Og[i, j] = x)
    end
    B = fill(FT(NaN), length(lons), length(lats))
    for i in eachindex(lons), j in eachindex(lats)
        valid(M[i, j]) && valid(Og[i, j]) && (B[i, j] = M[i, j] - Og[i, j])
    end
    # Payload format is the artifact's, not a convenience: a flat lat-major
    # array indexed `j*nx + i`, values scaled by 10 and stored as integers, -9999
    # for missing. Emitting it directly means the figure is redrawn from a
    # command rather than hand-edited, which is how a previous set of published
    # ratios came to be adjusted by eye.
    step = 2
    scale, miss = 10, -9999
    mm = block_means(lons, lats, M, step)
    oo = block_means(lons, lats, Og, step)
    bb = block_means(lons, lats, B, step)
    nx, ny = size(mm)
    function flat(A)
        v = Vector{Int}(undef, nx * ny)
        for j in 1:ny, i in 1:nx
            x = A[i, j]
            v[(j - 1) * nx + i] = isnan(x) ? miss : round(Int, x * scale)
        end
        return "[" * join(v, ",") * "]"
    end
    open(jsonpath, "w") do io
        print(io, "{\"nx\":$nx,\"ny\":$ny,\"scale\":$scale,\"miss\":$miss,")
        print(io, "\"model\":$(flat(mm)),")
        print(io, "\"obs\":$(flat(oo)),")
        print(io, "\"bias\":$(flat(bb))}")
    end
    println("wrote $jsonpath")
end

isempty(ARGS) &&
    error("usage: julia score_map.jl <equilibrium_carbon.nc> [maps.json]")
main(ARGS[1], length(ARGS) > 1 ? ARGS[2] : nothing)

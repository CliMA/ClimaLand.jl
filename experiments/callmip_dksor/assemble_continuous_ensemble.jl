#=
Stack the continuous per-member files (output_callmip_sims/members/member_<k>.jld2,
each a full 1997-2014 daily series) into the ensemble_diagnostics.jld2 that
write_callmip_netcdf.jl reads for the posterior _unc fields:
  member_nee/lhf/shf : (n_members × n_days),  dates : Vector{Date}.
All members share the same continuous date axis (prior/posterior too).
=#
using JLD2, Dates

const BASE = joinpath(@__DIR__, "output_callmip_sims")
const MDIR = joinpath(BASE, "members")

# NOTE: wrapped in a function — at top-level script scope, assigning `dates`
# inside the `for` loop makes it loop-local and `dates === nothing` throws
# UndefVarError (Julia soft-scope). A function barrier fixes it.
function assemble()
    files = sort(filter(f -> startswith(f, "member_") && endswith(f, ".jld2"), readdir(MDIR)),
                 by = f -> parse(Int, match(r"member_(\d+)\.jld2", f).captures[1]))
    isempty(files) && error("No member files in $MDIR")

    dates = nothing
    nee = Vector{Float64}[]; lhf = Vector{Float64}[]; shf = Vector{Float64}[]; pn = nothing
    for f in files
        d = JLD2.load(joinpath(MDIR, f))
        if dates === nothing
            dates = Date.(d["dates"])
        else
            @assert Date.(d["dates"]) == dates "member $f has a different date axis"
        end
        push!(nee, Float64.(d["nee"])); push!(lhf, Float64.(d["lhf"])); push!(shf, Float64.(d["shf"]))
        pn = d["param_names"]
    end
    NEE = permutedims(reduce(hcat, nee))   # (n_members × n_days)
    LHF = permutedims(reduce(hcat, lhf))
    SHF = permutedims(reduce(hcat, shf))
    jldsave(joinpath(BASE, "ensemble_diagnostics.jld2");
        member_nee = NEE, member_lhf = LHF, member_shf = SHF, dates = dates, param_names = pn)
    @info "Assembled $(length(files)) members × $(length(dates)) days → ensemble_diagnostics.jld2"
end
assemble()

#=
Build the deliverable posterior — posterior_mean.jld2 + posterior_samples.jld2 —
from the EKI checkpoint: the UTKI (impose_prior=true) final-ensemble mean and the
final ensemble itself form the posterior used for the CalLMIP submission and its
uncertainty bands.

Why EKI and not CES for the deliverable: the CalLMIP Scenario-1 objective is
NEE + LHF + SHF only, which constrains total ecosystem respiration but NOT the
autotrophic/heterotrophic split (a flat likelihood ridge). The regularized UTKI
posterior stays physical on that ridge (AR/GPP ~80%). CES (emulate_sample.jl) runs
as an independent cross-check — its Sample step anchored on this EKI posterior
(see priors.jl::build_dk_sor_eki_prior) — and converges to EKI; it is not the
deliverable because, being anchored on EKI, it largely reproduces it. See the PR
description for the full EKI-vs-CES comparison.
=#
import ClimaComms; ClimaComms.@import_required_backends
using JLD2, Statistics
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
include(joinpath(@__DIR__, "priors.jl"))

const OUT = joinpath(@__DIR__, "output_calibration")
prior, pv = build_dk_sor_priors()
pnames = [only(PD.get_name(d)) for d in pv]

ekp = JLD2.load(joinpath(OUT, "ekp_checkpoint.jld2"), "ekp")
@info "EKI loaded: $(length(EKP.get_u(ekp))-1) iterations, ensemble size $(EKP.get_N_ens(ekp))"

eki_mean = Float64.(EKP.get_ϕ_mean_final(prior, ekp))     # constrained, length 11
ens      = Float64.(EKP.get_ϕ_final(prior, ekp))          # constrained, 11 × N_ens

@info "EKI posterior mean (constrained):"
for (n, v) in zip(pnames, eki_mean)
    @info "  $(rpad(n,34)) = $(round(v; sigdigits=5))"
end

jldsave(joinpath(OUT, "posterior_mean.jld2");
    posterior_mean = eki_mean, param_names = pnames,
    eki_mean = eki_mean, coverage = NaN, used_eki_fallback = true)
jldsave(joinpath(OUT, "posterior_samples.jld2");
    constrained_posterior = ens, param_names = pnames)
@info "Wrote posterior_mean.jld2 + posterior_samples.jld2 (EKI fallback; $(size(ens,2)) ensemble members)"

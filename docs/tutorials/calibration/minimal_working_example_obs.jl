# ## Calibration: single site (Ozark) latent heat flux observations.

# We calibrate observations of latent heat flux (LHF) retrieved from Ozark eddy-covariance data.
# We calibrate three parameters: g1, which is a core parameters of stomatal conductance, and
# two parameters controling plant sensitivity to moisture stress, sc and pc.

# We use EnsembleKalmanProcess.jl to perform this calibration, detailed step by step below:

# We will need:
# 1. A function returning our model LHF output given the parameters we want to calibrate.
# 2. The "truth" target data, to calibrate on.
# 3. The prior distribution of these parameters.

# ## Load packages
import ClimaLandSimulations.Fluxnet as CLS # to run the model
import EnsembleKalmanProcesses as EKP # to perform the calibration
import Random # to use the same seed each run in the tutorial, optional
import Logging
Logging.disable_logging(Logging.Warn); # hide julia warnings

# ## Write a function returning our model LHF output given the parameters to calibrate
function Ozark_LHF(params) # params is a 2 element Array
    g1, sc, pc = params
    sv = CLS.run_fluxnet(
        "US-MOz";
        params = CLS.ozark_default_params(;
            conductance = CLS.conductance_ozark(; g1 = g1),
            photosynthesis = CLS.photosynthesis_ozark(; sc = sc, pc = pc),
        ),
    )[1]

    inputs = CLS.make_inputs_df("US-MOz")[1]

    simulation_output = CLS.make_output_df("US-MOz", sv, inputs)

    LHF_soil =
        [parent(sv.saveval[k].soil.turbulent_fluxes.lhf)[1] for k in 1:1441]
    LHF_canopy = [
        parent(sv.saveval[k].canopy.energy.turbulent_fluxes.lhf)[1] for
        k in 1:1441
    ]
    LHF = LHF_soil + LHF_canopy

    return LHF
end;

# ## "Truth" target data to calibrate on
t0 = Float64(120 * 3600 * 24);
N_spinup_days = 30;
N_days_sim = 30;
N_days = N_spinup_days + N_days_sim;
index_t_start = Int(t0 / (3600 * 24) * 48);
index_t_end = Int(index_t_start + (N_days - N_spinup_days) * 48);
inputs = CLS.make_inputs_df("US-MOz")[1];
LHF_target = Float64.(inputs.LE)[index_t_start:index_t_end];

# ## Parameters prior
prior_g1 = EKP.constrained_gaussian("g1", 141, 100, 0, Inf); # mean of 221 sqrt(Pa) = 7 sqrt(kPa), std of 100 (3 kPa)
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 5e-4, 0, Inf);
prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, Inf);
prior = EKP.combine_distributions([prior_g1, prior_sc, prior_pc]);

# To use the same seed each run, optional
rng_seed = 2
rng = Random.MersenneTwister(rng_seed)

# ## Calibration

# generate the initial ensemble and set up the ensemble Kalman inversion
N_ensemble = 3 # should be 30, 10*3 parameters
N_iterations = 3 # should be 10
Γ = 5.0 * EKP.I

initial_ensemble = EKP.construct_initial_ensemble(rng, prior, N_ensemble);

ensemble_kalman_process = EKP.EnsembleKalmanProcess(
    initial_ensemble,
    LHF_target,
    Γ,
    EKP.Inversion();
    rng = rng,
);

# We are now ready to carry out the inversion. At each iteration, we get the ensemble from the last iteration, apply
# Ozark_LHF(params) to each ensemble member, and apply the Kalman update to the ensemble.

params_i = [] # initialize empty array
ClimaLand_out = []
for i in 1:N_iterations
    push!(params_i, EKP.get_ϕ_final(prior, ensemble_kalman_process))
    push!(ClimaLand_out, [Ozark_LHF(params_i[i][:, j]) for j in 1:N_ensemble])
    ClimaLand_ens = hcat(ClimaLand_out[i]...)
    EKP.update_ensemble!(ensemble_kalman_process, ClimaLand_ens)
end;

# Done! Here are the parameters:

final_ensemble = EKP.get_ϕ_final(prior, ensemble_kalman_process)


# ## Plot
using CairoMakie
CairoMakie.activate!()
fig = Figure(size = (2000, 1000))
ax = Axis(
    fig[1:3, 1],
    ylabel = "Latent heat flux (W m^-2)",
    xlabel = "Half-hour",
)
trange = 1:1:length(LHF_target)

l2 = [
    lines!(
        ax,
        trange,
        ClimaLand_out[1][i],
        color = :red,
    ) for i in 1:N_ensemble
][1]
l3 = [
      lines!(ax, trange, ClimaLand_out[3][i], color = :green) for
    i in 1:N_ensemble
][1]
l1 = lines!(ax, trange, LHF_target, color = :black)
axislegend(
    ax,
    [l1, l2, l3],
    ["Observations (target)", "Before calibration", "After calibration"],
)
ax2 = Axis(fig[4, 1], ylabel = "Shortwave in (W m^-2)")
lines!(ax2, trange, Float64.(inputs.SW_IN)[index_t_start:index_t_end])
ax3 = Axis(fig[5, 1], ylabel = "Air temperature (K)")
lines!(ax3, trange, Float64.(inputs.TA)[index_t_start:index_t_end])
ax4 = Axis(fig[6, 1], ylabel = "Soil moisture (m^3 m^-3)")
lines!(ax4, trange, Float64.(inputs.SWC)[index_t_start:index_t_end])

save("fig_lhf_obs.png", fig);
# ![](fig_lhf_obs.png)

# with `N_ensemble` = 30 and `N_iterations` = 10, you get:

# ![calibration_obs](src/assets/calibration.png)

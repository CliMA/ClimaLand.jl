using ClimaLandSimulations
using ClimaLandSimulations.Fluxnet

# default parameters, except for custom parameter hetero_resp.kM_o2
sv_test, sol_test, Y_test, p_test = run_fluxnet(
    "US-MOz";
    params = ozark_default_params(;
        hetero_resp = hetero_resp_ozark(; kM_o2 = 0.003),
    ),
)

# defaults, except start date
sv, sol, Y, p = run_fluxnet(
    "US-MOz";
    setup = make_setup(;
        ozark_default_args(; t0 = Float64(125 * 1800 * 24))...,
    ),
)

# default parameters
sv, sol, Y, p = run_fluxnet("US-MOz")

# DataFrame for input and output
inputs, inputs_SI, inputs_commonly_used = make_inputs_df("US-MOz")
simulation_output = make_output_df("US-MOz", sv, inputs)

# Will save figures in current repo /figures
make_plots(inputs, simulation_output)

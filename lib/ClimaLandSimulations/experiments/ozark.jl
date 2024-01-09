using ClimaLandSimulations

# Should inputs and outputs be separate functions?
sv, sol, Y, p, inputs, climalsm = fluxnet_simulation("US-MOz");

# Will save figures in current repo /figures
make_plots(inputs, climalsm)

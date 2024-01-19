using ClimaLandSimulations

sv, sol, Y, p, inputs, climalsm = fluxnet_simulation("US-MOz");

# Will save figures in current repo /figures
make_plots(inputs, climalsm)

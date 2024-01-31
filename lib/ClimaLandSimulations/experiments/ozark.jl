using ClimaLandSimulations

sv, sol, Y, p, inputs, climaland = fluxnet_simulation("US-MOz");

# Will save figures in current repo /figures
make_plots(inputs, climaland)

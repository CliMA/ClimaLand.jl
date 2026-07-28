# ERA5 energy-flux leaderboard for a calibrated-model diagnostics directory.
# Compares the model's turbulent + radiative fluxes (lhf, shf, lwu, swu) against
# the ERA5 monthly product (era5_monthly_averages_surface_single_level_1979_2024)
# via the standard ERA5DataLoader — the energy-side analogue of the inversion
# carbon leaderboard (experiments/calibration/plot_inversion_leaderboard.jl).
#
# Usage (julia 1.10, .buildkite env has the CairoMakie/GeoMakie/ClimaAnalysis stack):
#   julia +1.10 --startup-file=no --project=.buildkite \
#     experiments/calibration/plot_energy_leaderboard.jl \
#     /home/alexisr/GitHub/climaland_out_invcal_iter0 [/path/to/savedir]

using CairoMakie, GeoMakie, ClimaAnalysis
using ClimaLand

const ext = Base.get_extension(ClimaLand, :LandSimulationVisualizationExt)
ext === nothing && error(
    "LandSimulationVisualizationExt is not loaded. Load CairoMakie, GeoMakie, ClimaAnalysis first.",
)

const diag_dir =
    length(ARGS) >= 1 ? ARGS[1] :
    "/home/alexisr/GitHub/climaland_out_invcal_iter0"
const savedir =
    length(ARGS) >= 2 ? ARGS[2] : mktempdir(; prefix = "energy_leaderboard_")
isdir(savedir) || mkpath(savedir)

@info "Diagnostics dir: $diag_dir"
@info "Save dir:        $savedir"

ext.compute_seasonal_leaderboard(savedir, diag_dir, "ERA5")

@info "Done. ERA5 energy leaderboard PNGs:"
for f in sort(readdir(savedir))
    endswith(f, ".png") && println("  ", joinpath(savedir, f))
end

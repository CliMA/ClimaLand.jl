# Seasonal leaderboard comparing the calibrated model diagnostics against the
# inversion-derived carbon obs (inversion NEE [CT2022], GOSIF GPP, residual ER)
# from the `inversion_nee` artifact plus the MODIS LAI target (Yuan et al. 2017).
# LAI replaces the earlier Hashimoto Rh panel, so the four panels are NEE / GPP /
# ER / LAI — the calibrated carbon + LAI targets. This is the inversion analogue
# of the ILAMB (FLUXCOM) leaderboard (uses InversionDataLoader; see
# ext/.../leaderboard/data_sources.jl).
#
# Usage (julia 1.10, .buildkite env has the CairoMakie/GeoMakie/ClimaAnalysis stack):
#   julia +1.10 --startup-file=no --project=.buildkite \
#     experiments/calibration/plot_inversion_leaderboard.jl \
#     /home/alexisr/GitHub/climaland_out_invcal_2 [/path/to/savedir]

using CairoMakie, GeoMakie, ClimaAnalysis
using ClimaLand

const ext = Base.get_extension(ClimaLand, :LandSimulationVisualizationExt)
ext === nothing && error(
    "LandSimulationVisualizationExt is not loaded. Load CairoMakie, GeoMakie, ClimaAnalysis first.",
)

const diag_dir =
    length(ARGS) >= 1 ? ARGS[1] : "/home/alexisr/GitHub/climaland_out_invcal_2"
const savedir =
    length(ARGS) >= 2 ? ARGS[2] : mktempdir(; prefix = "inversion_leaderboard_")
isdir(savedir) || mkpath(savedir)

@info "Diagnostics dir: $diag_dir"
@info "Save dir:        $savedir"

ext.compute_seasonal_leaderboard(savedir, diag_dir, "FlagshipCarbonMetrics")

@info "Done. Inversion leaderboard PNGs:"
for f in sort(readdir(savedir))
    endswith(f, ".png") && println("  ", joinpath(savedir, f))
end

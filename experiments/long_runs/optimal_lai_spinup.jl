# # Spin-up of the optimal-LAI model

# Spins up the time-integrated variables of `ZhouOptimalLAIModel` for the long runs,
# which start on March 1 2008: the global land model with prognostic LAI runs for six
# years of ERA5 (March 2002 to March 2008) with a 1-year memory of the annual totals,
# so that their initial error (from the climatology of the `optimal_lai_inputs`
# artifact, e.g. the potential evaporation of humid regions) decays by e^-6. Its final
# state is written to `optimal_lai_spinup_<device>/optimal_lai_state.nc`, the state
# the long runs start from (`ClimaLand.Artifacts.optimal_lai_state_path`). It takes
# about two hours on a GPU; run it again when the model changes enough to matter.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates

include(joinpath(@__DIR__, "optimal_lai_state.jl"))

const FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "optimal_lai_spinup_$(device_suffix)"
mkpath(root_path)

start_date = DateTime("2002-03-01")
stop_date = DateTime("2008-03-01")
Δt = 900.0
domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
toml_dict = LP.create_toml_dict(
    FT;
    override_files = [joinpath(@__DIR__, "optimal_lai_spinup.toml")],
)
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    context,
)
model = LandModel{FT}(
    (; atmos, radiation),
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
)
simulation = LandSimulation(start_date, stop_date, Δt, model; diagnostics = ())
@info "Optimal-LAI spin-up" start_date stop_date Δt domain.nelements
solve!(simulation)
write_optimal_lai_state(joinpath(root_path, "optimal_lai_state.nc"), simulation)

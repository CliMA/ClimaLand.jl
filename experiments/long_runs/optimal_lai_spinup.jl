# # Spin-up of the optimal-LAI model

# Spins up the time-integrated variables of `ZhouOptimalLAIModel` for runs that start
# from them, such as the long runs on March 1 2008: the global land model with prognostic
# LAI runs from March 2002, starting from the climatology of the `optimal_lai_inputs`
# artifact, with a 1-year memory of the annual totals, so that the error of the
# climatology (e.g. the potential evaporation of humid regions) decays by e^-6 by
# March 2008. The states of the first day of each month of 2008 are written to
# `optimal_lai_spinup_<device>/optimal_lai_state.nc`; a run starts from the state of
# the day of the year of its start date (`set_canopy_component_initial_conditions!`).
# It takes about two and a half hours on a GPU. Its output is the state of the
# `optimal_lai_state_and_clm_maps` artifact, to recreate when the model changes.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, step!
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
snapshot_dates = [DateTime(2008, month, 1) for month in 1:12]
stop_date = last(snapshot_dates)
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
set_ic! = ClimaLand.Simulations.make_set_initial_state_from_file(
    ClimaLand.Artifacts.rosetta_spunup_ic_path(; context),
    model;
    optimal_lai_state_path = nothing,
)
simulation =
    LandSimulation(start_date, stop_date, Δt, model; set_ic!, diagnostics = ())
@info "Optimal-LAI spin-up" start_date stop_date Δt domain.nelements
current_date() = ClimaUtilities.TimeManager.date(simulation._integrator.t)
snapshots = map(snapshot_dates) do snapshot_date
    while current_date() < snapshot_date
        step!(simulation)
    end
    @info "Optimal-LAI state" current_date()
    optimal_lai_state_snapshot(simulation)
end
write_optimal_lai_state(joinpath(root_path, "optimal_lai_state.nc"), snapshots)

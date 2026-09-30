# # Pre-industrial spin-up

# Initial conditions for the canopy carbon pools and the optimal-LAI model under a
# proxy for the pre-industrial climate: the earliest ERA5 years (1979-1982) with
# pre-industrial CO2. The global land model with prognostic LAI runs for three years,
# ending on March 1 as the long runs start. The optimal-LAI state is its final state;
# the carbon pools are set to their equilibrium under the climate of the last two years
# (see `equilibrium_biomass.jl`). The carbon pools do not feed back on GPP, LAI or
# climate, so they are not needed in the simulation itself.

# Output, in `preindustrial_spinup_<device>/`:
# - `initial_conditions.nc`: the carbon pools C_leaf, C_stem and C_root, their mean
#   annual temperature T_annual and precipitation P_annual, and the time-integrated
#   variables of `ZhouOptimalLAIModel`, on the lon-lat grid of the diagnostics
# - `equilibrium_woody_carbon.png`: the equilibrium woody carbon against XuSaatchi

import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaCore
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates

include(joinpath(@__DIR__, "equilibrium_biomass.jl"))

const FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "preindustrial_spinup_$(device_suffix)"
outdir = ClimaUtilities.OutputPathGenerator.generate_output_path(
    joinpath(root_path, "global_diagnostics"),
)

start_date = DateTime("1979-03-01")
stop_date = DateTime("1982-03-01")
Δt = 900.0
# Pre-industrial atmospheric CO2 (mol mol^-1)
c_co2 = TimeVaryingInput((t) -> 2.8e-4)

domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
toml_dict = LP.create_toml_dict(FT)
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    context,
    c_co2,
)
model = LandModel{FT}(
    (; atmos, radiation),
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
)
diagnostics = ClimaLand.default_diagnostics(
    model,
    start_date,
    outdir;
    output_vars = ["gpp", "crd", "ct", "tair", "precip", "fc3"],
)
simulation =
    LandSimulation(start_date, stop_date, Δt, model; outdir, diagnostics)
@info "Pre-industrial spin-up" start_date stop_date Δt domain.nelements
ClimaLand.Simulations.solve!(simulation)

skip_months = 12
parameters = ClimaLand.Canopy.PrognosticCarbonParameters(toml_dict)
lon, lat, pools = equilibrium_pools(outdir, parameters; skip_months)

# Final state of the optimal-LAI model, on the lon-lat grid of the diagnostics
Y = simulation._integrator.u
remapper = ClimaCore.Remapping.Remapper(
    axes(Y.canopy.biomass.LAI),
    [ClimaCore.Geometry.LatLongPoint(φ, λ) for λ in lon, φ in lat],
)
lai_names = ClimaLand.prognostic_vars(model.canopy.biomass)
lai_state = NamedTuple{lai_names}(
    map(lai_names) do name
        field = getproperty(Y.canopy.biomass, name)
        Array(ClimaCore.Remapping.interpolate(remapper, field))
    end,
)
write_initial_conditions(
    joinpath(root_path, "initial_conditions.nc"),
    lon,
    lat,
    (; pools..., lai_state...),
)
(; bias, rmse) =
    plot_equilibrium_woody_carbon(lon, lat, pools.C_stem; savedir = root_path)
@info "Equilibrium woody carbon against XuSaatchi (kg m^-2)" bias rmse

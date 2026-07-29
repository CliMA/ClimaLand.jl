# Compare the two forcing-loading paths on the SAME site data:
#
#   main : FLUXNET CSV -> 0-D `TimeVaryingInput(times, values)` on a `Column`
#   new  : FLUXNET-style NetCDF -> `DataSource` -> `MultiColumnDataHandler` ->
#          `TimeVaryingInput` on a single-column `ColumnEnsemble`
#
# The 0-D path is space-agnostic: it interpolates to one scalar and `fill!`s the
# destination field, so it runs on any space but gives every column the same value. The
# handler path reads one file per column onto a `PointCloudSpace` and is the only way for
# columns to differ. At N = 1 they should agree, so any difference here is attributable to
# the loading machinery rather than to the data or the domain.
#
# The NetCDF is generated from the site's own CSV by `fluxnet_netcdf.jl`, so the two paths
# read identical physical values and the machinery is what is left over. It does mean this
# is not a test against real PLUMBER2 data.
#
# For the same comparison through the full integrated fluxnet model rather than a bare soil
# column, see `compare_fluxnet_forcing_paths.jl`.
#
# Run: julia --project=.buildkite experiments/comparison/compare_forcing_paths.jl [SITE_ID]
#      (SITE_ID defaults to US-MOz; must be a site with a CSV and a ClimaLand site file)

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Dates
using DelimitedFiles # loads the FluxnetSimulations extension
using NCDatasets
using Printf

using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
using ClimaLand.Soil
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, step!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaParams as CP

include(joinpath(@__DIR__, "field_rmse.jl"))
include(joinpath(@__DIR__, "fluxnet_netcdf.jl"))

const FT = Float64
const SITE_ID = isempty(ARGS) ? "US-MOz" : ARGS[1]
const SITE_VAL = Val(FluxnetSimulations.replace_hyphen(SITE_ID))

toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)

(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, SITE_VAL)
(; time_offset, lat, long) = FluxnetSimulations.get_location(FT, SITE_VAL)
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, SITE_VAL)
(start_date, stop_date) = FluxnetSimulations.get_data_dates(SITE_ID, time_offset)

met_nc_path = joinpath(mkpath(joinpath(@__DIR__, "output")), "$(SITE_ID)_Met.nc")
isfile(met_nc_path) && rm(met_nc_path)
(; first_precip_date) = write_fluxnet_netcdf(
    met_nc_path,
    SITE_ID,
    lat,
    long,
    time_offset,
    thermo_params,
    FT,
)
@info "Generated NetCDF from CSV" met_nc_path SITE_ID time_offset lat long first_precip_date

# ---------------------------------------------------------------------------
# Build both forcings.
# ---------------------------------------------------------------------------
csv_forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
    SITE_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
)

column_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
ensemble_domain = ColumnEnsemble(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)

nc_forcing = FluxnetSimulations.prescribed_forcing_netcdf(
    met_nc_path,
    ensemble_domain.space.surface,
    atmos_h,
    start_date,
    toml_dict,
    FT;
    hour_offset_from_UTC = time_offset,
)

# ---------------------------------------------------------------------------
# Driver-by-driver comparison. The 0-D input writes into a plain 1-element array; the
# handler input writes into a Field on the point-cloud surface.
# ---------------------------------------------------------------------------
scalar_dest = zeros(FT, 1)
field_dest = ClimaCore.Fields.zeros(ensemble_domain.space.surface)
csv_at(tvi, t) = (ClimaLand.evaluate!(scalar_dest, tvi, t); scalar_dest[1])
nc_at(tvi, t) = (ClimaLand.evaluate!(field_dest, tvi, t); parent(field_dest)[1])

drivers = (
    ("T (K)", csv_forcing.atmos.T, nc_forcing.atmos.T),
    ("u (m/s)", csv_forcing.atmos.u, nc_forcing.atmos.u),
    ("q (kg/kg)", csv_forcing.atmos.q, nc_forcing.atmos.q),
    ("P (Pa)", csv_forcing.atmos.P, nc_forcing.atmos.P),
    ("c_co2 (mol/mol)", csv_forcing.atmos.c_co2, nc_forcing.atmos.c_co2),
    (
        "liquid_precip (m/s)",
        csv_forcing.atmos.liquid_precip,
        nc_forcing.atmos.liquid_precip,
    ),
    (
        "snow_precip (m/s)",
        csv_forcing.atmos.snow_precip,
        nc_forcing.atmos.snow_precip,
    ),
    ("SW_d (W/m^2)", csv_forcing.radiation.SW_d, nc_forcing.radiation.SW_d),
    ("LW_d (W/m^2)", csv_forcing.radiation.LW_d, nc_forcing.radiation.LW_d),
)

# Sampled at 15 minutes, off the 30-minute file grid on every other sample, so time
# interpolation is exercised and not just snapshot lookup. The window has to reach past the
# site's first precipitation event, or the composed `Precip`/`VPD` path stays identically
# zero and its agreement is vacuous; at US-MOz the first rain is at t = 900000 s, and it
# falls below freezing, so both the rain and the snow branch of the split are covered.
sweep_seconds = FT.(0:900:(12 * 86400))

@info "Driver comparison over $(length(sweep_seconds)) times" start_date
for (name, csv_tvi, nc_tvi) in drivers
    max_abs = FT(0)
    max_rel = FT(0)
    worst_t = FT(0)
    sum_csv = FT(0)
    for t in sweep_seconds
        csv_val = csv_at(csv_tvi, t)
        nc_val = nc_at(nc_tvi, t)
        d = abs(csv_val - nc_val)
        sum_csv += abs(csv_val)
        if d > max_abs
            max_abs = d
            worst_t = t
        end
        max_rel = max(max_rel, d / max(abs(csv_val), eps(FT)))
    end
    @printf(
        "%-22s max|Δ| = %.6e  max rel = %.6e  at t = %9.1f s   mean|csv| = %.6e\n",
        name,
        max_abs,
        max_rel,
        worst_t,
        sum_csv / length(sweep_seconds)
    )
end

# ---------------------------------------------------------------------------
# Same forcing difference, now through a model: a standalone EnergyHydrology soil
# column driven by each path.
# ---------------------------------------------------------------------------
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_vg_n,
    soil_vg_α,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_α_PAR,
    soil_α_NIR,
) = FluxnetSimulations.get_parameters(FT, SITE_VAL)

function build_soil_sim(domain, forcing)
    soil = Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        toml_dict;
        albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
            PAR_albedo = soil_α_PAR,
            NIR_albedo = soil_α_NIR,
        ),
        runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
        retention_parameters = (;
            ν = soil_ν,
            θ_r,
            K_sat = soil_K_sat,
            hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
        ),
        composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel),
        S_s = soil_S_s,
        z_0m = z_0m_soil,
        z_0b = z_0b_soil,
        emissivity = soil_ϵ,
    )
    function init_soil!(Y, p, t0, model)
        Y.soil.ϑ_l .= FT(0.4)
        Y.soil.θ_i .= 0
        ρc_s = Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            FT(0),
            model.parameters.ρc_ds,
            model.parameters.earth_param_set,
        )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                FT(0),
                ρc_s,
                FT(290),
                model.parameters.earth_param_set,
            )
    end
    dt = FT(450)
    return LandSimulation(
        start_date,
        start_date + Day(1),
        dt,
        soil;
        set_ic! = init_soil!,
        updateat = dt,
        diagnostics = nothing,
    )
end

sim_csv = build_soil_sim(column_domain, csv_forcing)
sim_nc = build_soil_sim(ensemble_domain, nc_forcing)

report_state_diffs(sim_csv, sim_nc; col = 1, label = "t = 0 (no steps)")

const NSTEPS = 100
for _ in 1:NSTEPS
    step!(sim_csv)
    step!(sim_nc)
end
@assert sum(isnan.(sim_nc._integrator.u)) == 0

report_state_diffs(sim_csv, sim_nc; col = 1, label = "after $NSTEPS steps")

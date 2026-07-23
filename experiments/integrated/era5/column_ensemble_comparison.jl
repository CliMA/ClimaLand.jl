# Run a ERA5 single site simulation on a single `Column` and a
# `ColumnEnsemble` of duplicate columns. We test that
# 1. The ensemble columns are identical to each other,
# 2. Every ensemble column matches the single column run within some tolerance.
# This involves checking the prognostic state, the cache, and the halfhourly
# diagnostics.

# As of now, the main differences between a single `Column` and a
# `ColumnEnsemble` of duplicate columns is the geometry. A single column uses
# a `ClimaCore.Geometry.CartesianGlobalGeometry` while a `ColumnEnsemble` use
# `ClimaCore.Geometry.SphericalGlobalGeometry`. Analytically, you would get the
# same results regardless of the geometry. However, due to floating point error,
# this results in different simulation states.

import ClimaComms
ClimaComms.@import_required_backends
using Dates
using Test
import ClimaDiagnostics

using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

include("comparison_utils.jl")

const FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
) where {FT}
    surface_space = domain.space.surface
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
        use_lowres_forcing = true,
    )
    forcing = (; atmos, radiation)

    land = LandModel{FT}(forcing, toml_dict, domain, Δt)

    output_vars = [
        "shf",
        "lhf",
        "trans",
        "swu",
        "lwu",
        "sr",
        "ssr",
        "precip",
        "et",
        "lai",
    ]
    writer = ClimaDiagnostics.Writers.DictWriter()
    diagnostics = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = writer,
        reduction_period = :halfhourly,
        reduction_type = :average,
        output_vars,
    )

    simulation = LandSimulation(start_date, stop_date, Δt, land; diagnostics)
    return simulation, writer
end

start_date = DateTime("2008-03-01")
stop_date = start_date + Day(30)
Δt = Float64(900)
longlat = FT.((-77.0, 0.1))
zlim = FT.((-15, 0))
nelements = 15
dz_tuple = FT.((3, 0.05))
toml_dict = LP.create_toml_dict(FT)

column_domain = Column(; zlim, nelements, dz_tuple, longlat)
N_columns = 3
ensemble_domain = ColumnEnsemble(;
    zlim,
    nelements,
    dz_tuple,
    longlat = fill(longlat, N_columns),
)

@info "Run: ERA5 Soil-Canopy-Snow Model on a Column and a ColumnEnsemble"
@info "Columns in the ensemble: $N_columns"
@info "Resolution: $(column_domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"

sim_column, writer_column =
    setup_model(FT, start_date, stop_date, Δt, column_domain, toml_dict)
@time solve!(sim_column)
sim_ensemble, writer_ensemble =
    setup_model(FT, start_date, stop_date, Δt, ensemble_domain, toml_dict)
@time solve!(sim_ensemble)

@testset "ERA5 Column vs ColumnEnsemble ($longlat)" begin
    for sim in (sim_column, sim_ensemble)
        @test sim._integrator.t == sim._integrator.sol.prob.tspan[2]
    end

    Y_column, p_column = sim_column._integrator.u, sim_column._integrator.p
    Y_ensemble, p_ensemble =
        sim_ensemble._integrator.u, sim_ensemble._integrator.p

    @testset "Duplicate ensemble columns are identical" begin
        # Only check the first and second columns of the ensemble
        k = 2
        diffs_Y =
            field_diffs(Y_ensemble, Y_ensemble; col1 = 1, col2 = k, name = "Y")
        diffs_p =
            field_diffs(p_ensemble, p_ensemble; col1 = 1, col2 = k, name = "p")
        report_diffs(
            diffs_Y;
            label = "Y: ensemble column $k vs ensemble column 1 from $start_date to $stop_date",
        )
        report_diffs(
            diffs_p;
            label = "p: ensemble column $k vs ensemble column 1 from $start_date to $stop_date",
        )

        # For two duplicate columns, we should get the exactly same result
        ens_rtol = 0.0
        for diffs in (diffs_Y, diffs_p)
            @test isempty(filter(((_, d),) -> !(d.err <= ens_rtol), diffs))
        end
    end

    k = 1
    @testset "Ensemble column $k vs single column" begin
        rtol = 1e-9
        diffs_Y = field_diffs(Y_column, Y_ensemble; col2 = k, name = "Y")
        diffs_p = field_diffs(p_column, p_ensemble; col2 = k, name = "p")
        diffs_diag = diagnostic_diffs(writer_column, writer_ensemble; col2 = k)
        report_diffs(
            diffs_Y;
            label = "Y: single vs ensemble column $k from $start_date to $stop_date",
        )
        report_diffs(
            diffs_p;
            label = "p: single vs ensemble column $k from $start_date to $stop_date",
        )
        report_diffs(
            diffs_diag;
            label = "Diagnostics (worst over saved times): single vs ensemble column $k",
        )
        for diffs in (diffs_Y, diffs_p, diffs_diag)
            @test isempty(filter(((_, d),) -> !(d.err <= rtol), diffs))
        end
    end
end

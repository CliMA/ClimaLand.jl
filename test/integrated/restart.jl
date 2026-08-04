using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore
import ClimaUtilities
import ClimaUtilities.TimeManager: ITime, date
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, step!
using Dates

# Integrate the global land model (soil, snow, soil CO2, inland water, and a
# canopy which uses the P-model for photosynthesis; this is the configuration of
# experiments/long_runs/snowy_land_pmodel.jl) for `n_steps`, saving a checkpoint
# after `n_steps_to_checkpoint`. Then restart from that checkpoint and integrate
# to the same final time, and compare against the uninterrupted simulation.
FT = Float64
context = ClimaComms.context()
start_date = DateTime(2008, 3, 1)
Δt = 900.0
n_steps = 18
# 3 hours, the interval at which `LandSimulation` updates the drivers by default.
# The checkpoint is taken on one of those updates so that the restarted
# simulation, whose update schedule starts at the checkpoint rather than at
# `start_date`, sees the same forcing as the original run.
n_steps_to_checkpoint = 12
stop_date = start_date + Second(n_steps * Δt)
checkpoint_date = start_date + Second(n_steps_to_checkpoint * Δt)
t0 = ITime(0, Second(1), start_date)

root_path = "land_pmodel_restart"
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(root_path)

# A coarse version of the long run: the restart machinery does not depend on
# resolution, and this keeps the solves cheap.
domain = ClimaLand.Domains.global_domain(
    FT;
    nelements = (10, 5),
    mask_threshold = FT(0.99),
    context,
)
toml_dict = LP.create_toml_dict(FT)
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = true,
    context,
)
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
)
model = ClimaLand.LandModel{FT}(
    (; atmos, radiation),
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
)
@test model.canopy.photosynthesis isa ClimaLand.Canopy.PModel

set_ic! = ClimaLand.Simulations.make_set_initial_state_from_file(
    ClimaLand.Artifacts.soil_ic_2008_50m_path(; context),
    model;
    enforce_constraints = true,
)

simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    model;
    set_ic!,
    user_callbacks = (
        ClimaLand.CheckpointCallback(
            ITime(n_steps_to_checkpoint * Δt, epoch = start_date),
            output_dir,
            t0;
            model,
        ),
    ),
    diagnostics = nothing,
)
for _ in 1:n_steps_to_checkpoint
    step!(simulation)
end
Y_checkpoint = deepcopy(simulation._integrator.u)
snow_T_sfc_checkpoint = deepcopy(simulation._integrator.p.snow.T_sfc)
for _ in 1:(n_steps - n_steps_to_checkpoint)
    step!(simulation)
end
@test date(simulation._integrator.t) == stop_date

restart_file = ClimaLand.find_restart(root_path)
@test !isnothing(restart_file)
t_restart = ClimaLand.initial_time_from_checkpoint(restart_file; model)
@test date(t_restart) == checkpoint_date

restarted_simulation() = LandSimulation(
    date(t_restart),
    stop_date,
    Δt,
    model;
    set_ic! = (Y, p, t, model) ->
        ClimaLand.set_initial_conditions_from_checkpoint!(
            Y,
            restart_file;
            model,
        ),
    user_callbacks = (),
    diagnostics = nothing,
)

"""
    prognostic_fields(Y, prefix = "")

Return the fields of the state `Y` paired with their names, e.g.
`"soil.ϑ_l" => Y.soil.ϑ_l`.
"""
function prognostic_fields(Y, prefix = "")
    fields = Pair{String, Any}[]
    for property_name in propertynames(Y)
        field = getproperty(Y, property_name)
        name = isempty(prefix) ? "$property_name" : "$prefix.$property_name"
        if field isa ClimaCore.Fields.Field
            push!(fields, name => field)
        else
            append!(fields, prognostic_fields(field, name))
        end
    end
    return fields
end

"""
    compare_states(expected, actual; rtol = 0)

Test that every field of the state `actual` matches `expected` to a relative
tolerance of `rtol`, comparing the largest deviation of a field against the
largest value it takes. With the default `rtol` the states must be identical.
"""
function compare_states(expected, actual; rtol = 0)
    for ((name, expected_field), (_, actual_field)) in
        zip(prognostic_fields(expected), prognostic_fields(actual))
        @testset "$name" begin
            expected_values = Array(parent(expected_field))
            actual_values = Array(parent(actual_field))
            @test maximum(abs.(expected_values .- actual_values)) <=
                  rtol * maximum(abs.(expected_values))
        end
    end
end

restart = restarted_simulation()
@testset "State read from the checkpoint" begin
    compare_states(Y_checkpoint, restart._integrator.u)
end

for _ in 1:(n_steps - n_steps_to_checkpoint)
    step!(restart)
end
@test date(restart._integrator.t) == stop_date
# The state evolves after the checkpoint, so the comparisons below are not
# satisfied trivially.
@test Array(parent(restart._integrator.u.soil.ρe_int)) !=
      Array(parent(Y_checkpoint.soil.ρe_int))

# Only the state is checkpointed, so the restarted simulation rebuilds its cache
# from it; the small differences that remain come entirely from the snow surface
# temperature, see below.
@testset "Final state after restarting" begin
    compare_states(simulation._integrator.u, restart._integrator.u; rtol = 1e-3)
end

# `p.snow.T_sfc` is the initial guess for the snow surface temperature root find,
# which takes a single Newton step per timestep. It is therefore lagged state
# rather than a function of `Y` alone, and it is the only quantity in the cache
# which is: restoring it, and nothing else, makes the restarted simulation
# identical to the uninterrupted one.
restart_with_snow_T_sfc = restarted_simulation()
parent(restart_with_snow_T_sfc._integrator.p.snow.T_sfc) .=
    parent(snow_T_sfc_checkpoint)
for _ in 1:(n_steps - n_steps_to_checkpoint)
    step!(restart_with_snow_T_sfc)
end
@testset "Final state after restarting with the snow surface temperature" begin
    compare_states(
        simulation._integrator.u,
        restart_with_snow_T_sfc._integrator.u,
    )
end

rm(root_path; recursive = true)

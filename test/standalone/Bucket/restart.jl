using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaLand
import ClimaParams
import SciMLBase
import ClimaTimeSteppers
import Dates
import ClimaUtilities

# First, run a simulation for 3 hours, saving a checkpoint every hour.
# Then, restart the simulation after hour 2 and complete the simulation up to hour 3.
# Compare the final states, they should be identical.
FT = Float32
t0 = 0.0
Δt = 3600.0
tf = 3Δt
start_date = Dates.DateTime(2005)
root_path = "bucket_restart"
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(root_path)

bucket_domain = ClimaLand.SphericalShell(;
    radius = FT(100),
    depth = FT(3.5),
    nelements = (1, 10),
    npolynomial = 1,
)

bucket_atmos, bucket_rad = ClimaLand.prescribed_analytic_forcing(FT)
τc = FT(1.0)
α_bareground_func = (coordinate_point) -> 0.2
α_snow = FT(0.8)
z_0m = FT(1e-2)
z_0b = FT(1e-3)
albedo = ClimaLand.Bucket.PrescribedBaregroundAlbedo{FT}(
    α_snow,
    α_bareground_func,
    bucket_domain.space.surface,
)
bucket_parameters =
    ClimaLand.Bucket.BucketModelParameters(FT; albedo, z_0m, z_0b, τc)

model = ClimaLand.Bucket.BucketModel(
    parameters = bucket_parameters,
    domain = bucket_domain,
    atmosphere = bucket_atmos,
    radiation = bucket_rad,
)

Y, p, coords = ClimaLand.initialize(model)
Y.bucket.T .= 280.0
Y.bucket.W .= 0.05
Y.bucket.Ws .= 0.0
Y.bucket.σS .= 0.08

exp_tendency! = ClimaLand.make_exp_tendency(model)
set_initial_cache! = ClimaLand.make_set_initial_cache(model)
set_initial_cache!(p, Y, t0)

prob = SciMLBase.ODEProblem(
    ClimaTimeSteppers.ClimaODEFunction((T_exp!) = exp_tendency!),
    Y,
    (t0, tf),
    p,
)
updateat = collect(t0:Δt:tf)
checkpoint_frequency = 2Δt
drivers = ClimaLand.get_drivers(model)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
checkpoint_cb = ClimaLand.CheckpointCallback(
    checkpoint_frequency,
    output_dir,
    start_date,
    t0;
    model,
)
cb = SciMLBase.CallbackSet(driver_cb, checkpoint_cb)

timestepper = ClimaTimeSteppers.RK4()
ode_algo = ClimaTimeSteppers.ExplicitAlgorithm(timestepper)

sol = SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb)

# Now, let's restart from the checkpoint (we have the pass the root_path, not the output_dir)
restart_file = ClimaLand.find_restart(root_path)
Y_restart, p_restart, _ = ClimaLand.initialize(model)
ClimaLand.set_initial_conditions_from_checkpoint!(
    Y_restart,
    restart_file;
    model,
)
t_restart = ClimaLand.initial_time_from_checkpoint(restart_file; model)

set_initial_cache!(p_restart, Y_restart, t_restart)
prob_restart = SciMLBase.ODEProblem(
    ClimaTimeSteppers.ClimaODEFunction((T_exp!) = exp_tendency!),
    Y_restart,
    (t_restart, tf),
    p_restart,
)
updateat_restarted = collect(t_restart:Δt:tf)
driver_cb_restarted =
    ClimaLand.DriverUpdateCallback(updateat_restarted, updatefunc)
cb_restarted = SciMLBase.CallbackSet(driver_cb_restarted)
sol_restarted =
    SciMLBase.solve(prob_restart, ode_algo; dt = Δt, callback = cb_restarted)
for p in propertynames(Y.bucket)
    @test getproperty(Y.bucket, p) == getproperty(Y_restart.bucket, p)
end

rm(output_dir; recursive = true)

import ClimaCalibrate as CAL

# CAL.forward_model(iteration, member)
# --> just runs the model, doesn't return anything
# CAL.observation_map(iteration)
# --> returns the vector of obs

# @everywhere, see https://github.com/CliMA/ClimaAtmos.jl/blob/30ac26aafbb27c66bc65201dec5397116588a9af/calibration/test/e2e_test.jl

CAL.calibrate(
              CAL.WorkerBackend,
              ensemble_size,
              n_iterations,
              observations,
              noise,
              prior,
              output_dir
             )


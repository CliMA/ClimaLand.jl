import TOML as toml

# Explain things here - tutorial style
"""

"""
function CAL.forward_model(iteration, member)
    ensemble_member_path = path_to_ensemble_member(caldir, iteration, member)
    params_path = parameter_path(caldir, iteration, member)
    params = toml.parsefile(params_path) # should load a Dict, that needs to be converted to namedtuple

    @info ensemble_member_path
    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    diagdir = ClimaUtilities.OutputPathGenerator.generate_output_path(
        diagnostics_dir,
    )

    prob, cb = setup_prob(t0, tf, Δt, params, diagdir; nelements)

    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
    return nothing
end

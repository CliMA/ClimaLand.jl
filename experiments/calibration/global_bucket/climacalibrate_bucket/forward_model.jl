# Explain things here - tutorial style
# (Nat: why is it CAL.forward_model, and not CAL.parameter_path?
"""

"""
function CAL.forward_model(iteration, member)
    ensemble_member_path = parameter_path(caldir, iteration, member)
    params = toml.parsefile(ensemble_member_path) # should load a Dict, that needs to be converted to namedtuple

    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    diagdir = ClimaUtilities.OutputPathGenerator.generate_output_path(
        diagnostics_dir,
    )

    prob, cb = setup_prob(t0, tf, Δt, params, diagdir; nelements)

    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
end

export make_timestepper

"""
    make_timestepper(site_setup_out;
        N_spinup_days = 30,
        N_days_sim = 30,
        timestepper = CTS.RK4(),
        ode_algo = CTS.ExplicitAlgorith(timestepper)
        )

Define the setup for the simulation timestepper. 
the default timestepper is 4th order Runge Kutta method,
 other timesteppers from ClimaTimeSteppers.jl can be used. 
"""
function make_timestepper(
    site_setup_out;
    N_spinup_days = 30,
    N_days_sim = 30,
    timestepper = CTS.RK4(),
    ode_algo = CTS.ExplicitAlgorithm(timestepper),
)
    N_days = N_spinup_days + N_days_sim
    tf = Float64(site_setup_out.t0 + 3600 * 24 * N_days)
    t_spinup = Float64(site_setup_out.t0 + N_spinup_days * 3600 * 24)
    saveat = Array(t_spinup:(site_setup_out.n * site_setup_out.dt):tf)
    return (
        N_days = N_days,
        tf = tf,
        t_spinup = t_spinup,
        saveat = saveat,
        timestepper = timestepper,
        ode_algo = ode_algo,
    )
end

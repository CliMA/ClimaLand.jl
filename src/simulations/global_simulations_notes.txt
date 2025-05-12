function GlobalDomain(...)
    return SphericalShell()
end

context = ;
FT = ;
param_set = ...earth_param_set(FT)
parameter_maps = spatially_varying_maps(FT, ...)

domain = GlobalDomain(context, )
land = LandModel(FT, domain, params, component_choices; parameter_maps)
# if hires fails to download, use lowres (era5 forcing?)


#set_ic!(Y, land, t0)
timestepper = (alg::ImEx = ARS111, t0, dt, tf; start_date = nothing) # handled internally with itime, e.g. pass in seconds
timestepper = (alg::ImEx = ARS111, t0::DateTime, dt::Seconds, duration::Period) # handled internally with itime
user_callbacks # all SciML callbacks, but we provide helper functions for specific instances (nancheck, report, checkpoint)
diagnostics # Tuple of scheduled diagnostics, diagnostics handler is constructed internally
# driver is combine with user_callbacks and diagnostics callback internally to make the CB set
output_dir # will have default
#context  # will have default, this includes the device, or take the device and create the context
config = # specifies a default postprocessing/plots to make
    simulation = Simulation(land, set_ic!, domain, timestepper, user_callbacks, diagnostics, config)



   


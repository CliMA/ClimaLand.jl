function param_dashboard(model_parameters, # e.g., SoilCO2ModelParameters
                         model_functions::Vector{Function}, # e.g., [microbe_source, co2_diffusivity]
                         model_drivers, # e.g., soil temperature, soil moisture
                         drivers_limit # e.g., ([273, 303], [0.0, 0.8]) 
                         ) 

  # Figure, Axis3 (3D) and Axis (2D), Menu
  fig = Figure(resolution = (1200, 1200))
  labels_d = ["Soil temperature [K]", "Soil moisture [m³ m⁻³]"] # Can be retrieved from PrecribedSoil struct (would be a different struct for prognostic) 
  menu_opt = ["CO2 production", "CO2 diffusivity"] # Should get these name directly from AbtractModel param function names
  menu = Menu(fig[1,1], options = menu_opt); m = menu.selection
  ax3D = Axis3(fig[2,1], xlabel = labels_d[1], ylabel = labels_d[2])
  axT = Axis(fig[3,1], xlabel = labels_d[1])
  axM = Axis(fig[4,1], xlabel = labels_d[2])

  # SliderGrid for parameters
  FT = Float64 
  earth_param_set = create_lsm_parameters(FT)
  params = SoilCO2ModelParameters{FT}(; earth_param_set = earth_param_set)
  labels = ["$(s)" for s in fieldnames(SoilCO2ModelParameters)[1:end-1]] # without earth_param_set
  ranges = [(val/2 : val/4: val*2, val) for val in [getfield(params, i) for i in fieldnames(SoilCO2ModelParameters)[1:end-1]]]
  sliders = [(label = label, range = range, startvalue = startvalue) for (label, (range, startvalue)) in zip(labels, ranges)]
  sg = SliderGrid(fig[2, 2], sliders..., width = 350, tellheight = true)
  param_title = Label(fig[1, 2], "Parameters", tellheight = false)

  # SliderGrid for drivers
  startval_d = [293.13, 0.3] # Again, should this values taken from AbstractModel somehow, or given by user with an arg? 
  ranges_d = [(val/2 : val/4: val*2, val) for val in startval_d]
  sliders_d = [(label = label, range = range, startvalue = startvalue) for (label, (range, startvalue)) in zip(labels_d, ranges_d)]
  sg_d = SliderGrid(fig[4, 2], sliders_d..., width = 350, tellheight = false)
  param_title_d = Label(fig[3, 2], "Drivers", tellheight = false)

  # Get Observable and their values from SliderGrid
  sd = Dict(i => sg.sliders[i].value for i in 1:length(sliders))
  sd_v = Dict(i => sg.sliders[i].value[] for i in 1:length(sliders))
  sd_d = Dict(i => sg_d.sliders[i].value for i in 1:length(sliders_d))
  sd_v_d = Dict(i => sg_d.sliders[i].value[] for i in 1:length(sliders_d))

  # Create struct of AbtractModel (example SoilCO2ModelParameters for now, but should be an argument of function) from SliderGrid values
  keys = Symbol.(labels) 
  s = collect(values(sort(sd_v)))
  args = (; zip(keys, s)...)
  parameters = Observable(SoilCO2ModelParameters{FT}(; args..., earth_param_set = earth_param_set)) # this works

  # Menu select the model to display. x is driver #1 y is driver #2, p is parameters. 
  # This is a tricky part. Some parameterization will depend on 1 or 2 drivers, and sometimes on more stuff (e.g. below, Csom is not a param nor driver)
  # Need to find a way to generalize a solution to all parameterization (could be specified by users in args)
  fun = Dict(menu_opt[1] => (x, y, p) -> microbe_source(x, y, Csom, p), # need to give a value for Csom
             menu_opt[2] => co2_diffusivity) # Note, name of function should be retrieved from the AbstractModel (how? Kat)

  Csom = FT(5.0) # Should be either retrieved from AbstractModel or given by user as an arg
  x = @lift(mat([283, 313], [0.0, 0.5], 30, fun[$m], $parameters)[1]) 
  y = @lift(mat([283, 313], [0.0, 0.5], 30, fun[$m], $parameters)[2])
  z = @lift(mat([283, 313], [0.0, 0.5], 30, fun[$m], $parameters)[3])
  surface!(ax3D, x, y, z)

  # plot on 2D axis
  x_axT = collect(283:1:313) # should be ax1 and ax2, drivers can be other than temperature or moisture
  x_axM = collect(0.0:0.0163:0.5)

  # special stuff for driver 1 and driver 2 (e.g., temperature and moisture)
  # so need a separate on(sd[drivers[i]])
  y_axT = @lift(fvec(x_axT, $(sd_d[2]), fun[$m], $parameters)) # fun should be Observable from Menu
  y_axM = @lift(fvecM($(sd_d[1]), x_axM, fun[$m], $parameters))

  lines!(axT, x_axT, y_axT, color = :red, linewidth = 4)
  lines!(axM, x_axM, y_axM, color = :blue, linewidth = 4)

  cM = @lift(repeat([$(sd_d[2])], 31)) # Should output an Observable? See old code from on the on(sd[i]) loop
  cT = @lift(repeat([$(sd_d[1])], 31))
  lines!(ax3D, x_axT, cM, y_axT, color = :red, linewidth = 4) 
  lines!(ax3D, cT, x_axM, y_axM, color = :blue, linewidth = 4)

  for i in 1:length(sd)
    on(sd[i]) do val # update parameters and rescale x and y limits
      sd_v = Dict(i => sg.sliders[i].value[] for i in 1:length(sliders))
      s = collect(values(sort(sd_v)))
      args = (; zip(keys, s)...)
      parameters[] = SoilCO2ModelParameters{FT}(; args..., earth_param_set = earth_param_set)  # new args
      autolimits!(ax3D)
      autolimits!(axT)
      autolimits!(axM)
    end
  end

  for i in 1:2
    on(sd_d[i]) do val
      autolimits!(axT)
      autolimits!(axM)
    end
  end

  on(menu.selection) do val
    autolimits!(ax3D)
    autolimits!(axT)
    autolimits!(axM)
  end

  return fig
end


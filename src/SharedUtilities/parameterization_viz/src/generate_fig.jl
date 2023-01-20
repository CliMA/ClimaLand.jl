#=
using parameterization_viz # local
using ClimaLSM, ClimaLSM.Soil.Biogeochemistry

model_parameters = SoilCO2ModelParameters
model_functions = Dict("CO2 production" => (d1, d2, p) -> microbe_source(d1, d2, 5.0, p),
                       "CO2 diffusivity" => co2_diffusivity)
drivers_name = ["T_soil", "M_soil"]
drivers_limit = ([273, 303], [0.0, 0.5])

param_dashboard(model_parameters, model_functions, drivers_name, drivers_limit)
=#

"""
    param_dashboard(model_parameters, model_functions, drivers_name, drivers_limit)

Generates an interactive web dashboard of the parameterization functions of a model. 
"""
function param_dashboard(model_parameters, model_functions, drivers_name, drivers_limit) 
  # Create a Figure and its layout: a Menu, a 3D axis, and 2 2D axis 
  fig = Figure(resolution = (1200, 1200))
  menu_opt = collect(keys(model_functions)) 
  menu = Menu(fig[1,1:2], options = menu_opt); m = menu.selection
  ax3D = Axis3(fig[2,2], xlabel = drivers_name[1], ylabel = drivers_name[2])
  ax_d1 = Axis(fig[3,1], xlabel = drivers_name[1])
  ax_d2 = Axis(fig[3,2], xlabel = drivers_name[2])

  # Get SliderGrid args for parameters
  FT = Float64 
  earth_param_set = create_lsm_parameters(FT)
  params = model_parameters{FT}(; earth_param_set = earth_param_set)
  labels = ["$(s)" for s in fieldnames(model_parameters)[1:end-1]] # without earth_param_set
  ranges = [(val/2 : val/4: val*2, val) for val in [getfield(params, i) for i in fieldnames(model_parameters)[1:end-1]]]
  sliders = [(label = label, range = range, startvalue = startvalue) for (label, (range, startvalue)) in zip(labels, ranges)]

  # Get SliderGrid args for drivers
  startval_d = [mean(drivers_limit[1]), mean(drivers_limit[2])]  
  ranges_d = [(val/2 : val/4: val*2, val) for val in startval_d]
  sliders_d = [(label = label, range = range, startvalue = startvalue) for (label, (range, startvalue)) in zip(drivers_name, ranges_d)]
  
  # Layouting sliders
  s_layout = GridLayout()
  param_title = s_layout[1,1] = Label(fig, "Parameters")
  sg = s_layout[2,1] = SliderGrid(fig, sliders..., width = 250)
  param_title_d = s_layout[1,2] = Label(fig, "Drivers")
  sg_d = s_layout[2,2] = SliderGrid(fig, sliders_d..., width = 220)
  fig.layout[2,1] = s_layout

  # Get Observable and their values from SliderGrid
  sd = Dict(i => sg.sliders[i].value for i in 1:length(sliders))
  sd_v = Dict(i => sg.sliders[i].value[] for i in 1:length(sliders))
  sd_d = Dict(i => sg_d.sliders[i].value for i in 1:length(sliders_d))
  sd_v_d = Dict(i => sg_d.sliders[i].value[] for i in 1:length(sliders_d))

  # Create struct of parameters from SliderGrid values
  param_keys = Symbol.(labels) 
  s = collect(values(sort(sd_v)))
  args = (; zip(param_keys, s)...)
  parameters = Observable(model_parameters{FT}(; args..., earth_param_set = earth_param_set)) 

  # Plot 3D surface of model(drivers, params)
  x = @lift(mat(drivers_limit[1], drivers_limit[2], 30, model_functions[$m], $parameters)[1]) 
  y = @lift(mat(drivers_limit[1], drivers_limit[2], 30, model_functions[$m], $parameters)[2])
  z = @lift(mat(drivers_limit[1], drivers_limit[2], 30, model_functions[$m], $parameters)[3])
  surface!(ax3D, x, y, z, colormap = Reverse(:Spectral), transparency = true, alpha = 0.2, shading = false)

  # Plot 2D lines of model(drivers, params)
  x_d1 = collect(range(drivers_limit[1][1], drivers_limit[1][2], 31)) 
  x_d2 = collect(range(drivers_limit[2][1], drivers_limit[2][2], 31))
  y_d1 = @lift(d1_vec(x_d1, $(sd_d[2]), model_functions[$m], $parameters)) 
  y_d2 = @lift(d2_vec($(sd_d[1]), x_d2, model_functions[$m], $parameters))
  lines!(ax_d1, x_d1, y_d1, color = :red, linewidth = 4)
  lines!(ax_d2, x_d2, y_d2, color = :blue, linewidth = 4)

  # Plot 3D lines of model(drivers, params)
  c_d2 = @lift(repeat([$(sd_d[2])], 31)) 
  c_d1 = @lift(repeat([$(sd_d[1])], 31))
  lines!(ax3D, x_d1, c_d2, y_d1, color = :red, linewidth = 4) 
  lines!(ax3D, c_d1, x_d2, y_d2, color = :blue, linewidth = 4)

  # Update parameters and rescale x and y limits
  for i in 1:length(sd)
    on(sd[i]) do val 
      sd_v = Dict(i => sg.sliders[i].value[] for i in 1:length(sliders))
      s = collect(values(sort(sd_v)))
      args = (; zip(param_keys, s)...)
      parameters[] = model_parameters{FT}(; args..., earth_param_set = earth_param_set)  # new args
      autolimits!(ax3D)
      autolimits!(ax_d1)
      autolimits!(ax_d2)
    end
  end
  for i in 1:2
    on(sd_d[i]) do val
      autolimits!(ax_d1)
      autolimits!(ax_d2)
    end
  end
  on(menu.selection) do val
    autolimits!(ax3D)
    autolimits!(ax_d1)
    autolimits!(ax_d2)
  end

  DataInspector(fig)
  return fig
end


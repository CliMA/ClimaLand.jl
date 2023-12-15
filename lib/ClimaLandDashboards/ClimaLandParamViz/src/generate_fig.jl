"""
    param_dashboard(parameterisation::Function, inputs::Inputs, sliders)

Generates a dashboard of a parameterisation(drivers, parameters, constants) function,
where the user can interact with driver and parameter values via sliders. 
"""
function param_dashboard(parameterisation::Function, inputs::Inputs, drivers_sliders, parameters_sliders, output)
  fig = Figure(resolution = (800, 800))

  drivers_ranges_unitconverted = [ustrip.(uconvert.(inputs.drivers.units[i][2], (inputs.drivers.ranges[i])inputs.drivers.units[i][1])) for i = 1:2]
  parameters_ranges_unitconverted = [ustrip.(uconvert.(inputs.parameters.units[i][2], (inputs.parameters.ranges[i])inputs.parameters.units[i][1])) for i = 1:length(inputs.parameters.units)]
  output_range_unitconverted = ustrip.(uconvert.(output.unit[2], (output.range)output.unit[1]))

  # JSServe layout
  ax3D = Axis3(fig[1,1:2][1,1], xlabel = inputs.drivers.names[1], ylabel = inputs.drivers.names[2], zlabel = output.name); zlims!(ax3D, output_range_unitconverted);
  xlims!(ax3D, drivers_ranges_unitconverted[1][1], drivers_ranges_unitconverted[1][2]); ylims!(ax3D, drivers_ranges_unitconverted[2][1], drivers_ranges_unitconverted[2][2])
  ax_d1 = Axis(fig[2,1], xlabel = inputs.drivers.names[1], ylabel = output.name); ylims!(ax_d1, output_range_unitconverted); xlims!(ax_d1, drivers_ranges_unitconverted[1])
  ax_d2 = Axis(fig[2,2], xlabel = inputs.drivers.names[2], ylabel = output.name); ylims!(ax_d2, output_range_unitconverted); xlims!(ax_d2, drivers_ranges_unitconverted[2])

  n_drivers = 2  
  n_parameters = length(inputs.parameters.names)
  drivers_vals = [@lift(ustrip.(uconvert.(inputs.drivers.units[i][1], ($(drivers_sliders[i].value))*inputs.drivers.units[i][2]))) for i in 1:n_drivers] |> Tuple
  parameters_vals = [@lift(ustrip.(uconvert.(inputs.parameters.units[i][1], ($(parameters_sliders[i].value))*inputs.parameters.units[i][2]))) for i in 1:n_parameters] |> Tuple 

  steps = 30
  drivers = lift((args...,) -> args, drivers_vals...)
  parameters = lift((args...,) -> args, parameters_vals...)
  constants = inputs.constants.values

  x_d1 = collect(range(drivers_ranges_unitconverted[1][1], drivers_ranges_unitconverted[1][2], steps)) # min d1 to max d1, n steps
  x_d2 = collect(range(drivers_ranges_unitconverted[2][1], drivers_ranges_unitconverted[2][2], steps)) # min d2 to max d2, n steps
  c_d1 = @lift(repeat([ustrip.(uconvert.(inputs.drivers.units[2][2], ($(drivers_vals[2])*inputs.drivers.units[2][1])))], steps)) # constant d1 val, length n
  c_d2 = @lift(repeat([ustrip.(uconvert.(inputs.drivers.units[1][2], ($(drivers_vals[1])*inputs.drivers.units[1][1])))], steps)) # constant d2 val, length n
  y_d1 = @lift(ustrip.(uconvert.(output.unit[2], (d1_vec($(drivers_vals[2]), parameterisation, inputs, $parameters, steps))*output.unit[1]))) # function output at constant driver 1
  y_d2 = @lift(ustrip.(uconvert.(output.unit[2], (d2_vec($(drivers_vals[1]), parameterisation, inputs, $parameters, steps))*output.unit[1]))) # function output at constant driver 2
  val = @lift(ustrip.(uconvert.(output.unit[2], (parameterisation($drivers, $parameters, constants))*output.unit[1])))
  point3D = @lift(Vec3f.(ustrip.(uconvert.(inputs.drivers.units[1][2], ($(drivers_vals[1]))*inputs.drivers.units[1][1])), ustrip.(uconvert.(inputs.drivers.units[2][2], ($(drivers_vals[2]))*inputs.drivers.units[2][1])), $val))
  point2D_ax1 = @lift(Vec2f.(ustrip.(uconvert.(inputs.drivers.units[1][2], ($(drivers_vals[1]))*inputs.drivers.units[1][1])), $val))
  point2D_ax2 = @lift(Vec2f.(ustrip.(uconvert.(inputs.drivers.units[2][2], ($(drivers_vals[2]))*inputs.drivers.units[2][1])), $val))

  # Plot 3D surface of model(drivers, params)
  
  x = @lift(ustrip.(uconvert.(inputs.drivers.units[1][2], (mat(parameterisation, inputs, $parameters, steps)[1])*inputs.drivers.units[1][1])))
  y = @lift(ustrip.(uconvert.(inputs.drivers.units[2][2], (mat(parameterisation, inputs, $parameters, steps)[2])*inputs.drivers.units[2][1])))
  z = @lift(ustrip.(uconvert.(output.unit[2], (mat(parameterisation, inputs, $parameters, steps)[3])*output.unit[1])))
  surface!(ax3D, x, y, z, colormap = Reverse(:Spectral), transparency = true, alpha = 0.8, shading = false, colorrange = output_range_unitconverted)
  cb = Colorbar(fig[1, 1:2][1, 2], colormap = Reverse(:Spectral), limits = output_range_unitconverted, label = output.name)
  cb.alignmode = Mixed(right = 0)

  # Plot 2D lines of model(drivers, params)
  lines!(ax_d1, x_d1, y_d1, color = :red, linewidth = 4)
  lines!(ax_d2, x_d2, y_d2, color = :blue, linewidth = 4)
  scatter!(ax_d1, point2D_ax1, color = :black, markersize = 20)
  scatter!(ax_d2, point2D_ax2, color = :black, markersize = 20)

  # Plot 3D lines of model(drivers, params)
  lines!(ax3D, x_d1, c_d1, y_d1, color = :red, linewidth = 4) 
  lines!(ax3D, c_d2, x_d2, y_d2, color = :blue, linewidth = 4)
  scatter!(ax3D, point3D, color = :black, markersize = 20, colormap = Reverse(:Spectral), colorrange = output.range,
          strokewidth = 10, strokecolor = :black) # stroke not supported in WGLMakie?

  DataInspector(fig)

  return fig, val  
end

function webapp(parameterisation, inputs, output)
  Param_app = App() do
    drivers_ranges_unitconverted = [ustrip.(uconvert.(inputs.drivers.units[i][2], (inputs.drivers.ranges[i])inputs.drivers.units[i][1])) for i = 1:2]
    parameters_ranges_unitconverted = [ustrip.(uconvert.(inputs.parameters.units[i][2], (inputs.parameters.ranges[i])inputs.parameters.units[i][1])) for i = 1:length(inputs.parameters.units)]
    output_range_unitconverted = ustrip.(uconvert.(output.unit[2], (output.range)output.unit[1]))
    n_drivers = 2  
    n_parameters = length(inputs.parameters.names)
    drivers_range = [round.(range(drivers_ranges_unitconverted[i][1], drivers_ranges_unitconverted[i][2], 12), sigdigits = 2) for i in 1:n_drivers]
    parameters_range = [round.(range(parameters_ranges_unitconverted[i][1], parameters_ranges_unitconverted[i][2], 12), sigdigits = 2) for i in 1:n_parameters]
    drivers_sliders = [JSServe.TailwindDashboard.Slider(inputs.drivers.names[i], drivers_range[i], value = drivers_range[i][6]) for i in 1:n_drivers] |> Tuple
    parameters_sliders = [JSServe.TailwindDashboard.Slider(inputs.parameters.names[i], parameters_range[i], value = parameters_range[i][6]) for i in 1:n_parameters] |> Tuple
    fig, out = param_dashboard(parameterisation, inputs, drivers_sliders, parameters_sliders, output)
    output_value = DOM.div(output.name, " = ", @lift(round($(out), sigdigits = 2)); style="font-size: 20px; font-weight: bold")
    drivers_label = DOM.div("Drivers:"; style="font-size: 16px; font-weight: bold")
    parameters_label = DOM.div("Parameters:"; style="font-size: 16px; font-weight: bold")
    return DOM.div(
                   JSServe.TailwindDashboard.Card(
                   JSServe.TailwindDashboard.FlexCol(
                                                     JSServe.TailwindDashboard.Card(output_value; class="container mx-auto"),
                                                     JSServe.TailwindDashboard.FlexRow(
                                                                                       JSServe.TailwindDashboard.Card(JSServe.TailwindDashboard.FlexCol(parameters_label, parameters_sliders...)),
                                                                                       JSServe.TailwindDashboard.Card(JSServe.TailwindDashboard.FlexCol(drivers_label, drivers_sliders...))
                                                                                      ),
                                                     fig)      
                                                    )
                  )
  end
  return Param_app
end


#= with units
using ParamViz
using JSServe
using WGLMakie

using Unitful: R, L, mol, K, kJ, °C, m, g, cm, hr, mg, s, μmol
using UnitfulMoles: molC
using Unitful, UnitfulMoles
@compound CO₂

using ClimaLSM
using ClimaLSM.Canopy
FT = Float64

drivers = Drivers(("PAR (μmol m⁻² s⁻¹))", "LAI (m² m⁻²))"),
                  (FT.([0, 1500 * 1e-6]), FT.([0, 10])),
                  ((mol*m^-2*s^-1, μmol*m^-2*s^-1), (m^2*m^-2, m^2*m^-2))
                 )

parameters = Parameters(("canopy reflectance, ρ_leaf",
                         "extinction coefficient, K",
                         "clumping index, Ω"),
                        (FT.([0, 1]), FT.([0, 1]), FT.([0, 1])),
                        ((m, m), (m, m), (m, m)) # dummy units, no conversion
                       )

# need a method with no constant! 
# hack: useless constants
constants = Constants(("a", "b"), (FT(1), FT(2)))

inputs = Inputs(drivers, parameters, constants)

output = Output("APAR (μmol m⁻² s⁻¹)", [0, 1500 * 1e-6], (mol*m^-2*s^-1, μmol*m^-2*s^-1))

import ParamViz.parameterisation
function parameterisation(PAR, LAI, ρ_leaf, K, Ω, a, b)   
  APAR = plant_absorbed_ppfd(PAR, ρ_leaf, K, LAI, Ω) 
  return APAR
end

beer_app = webapp(parameterisation, inputs, output)

=#

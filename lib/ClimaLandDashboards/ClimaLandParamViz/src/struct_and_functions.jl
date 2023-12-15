"""
  Drivers(names, ranges, units)

A struct to store drivers names, ranges and units
"""
struct Drivers # Need 2 drivers always
  names 
  ranges
  units
end

"""
  Parameters(names, ranges, units)

A struct to store parameters names, ranges and units
"""
struct Parameters
  names
  ranges
  units
end

"""
  Constant(names, values)

A struct to store constants names and values
"""
struct Constants # should be able to be 0, or have a way to deal with 0
  names
  values
end

"""
  Inputs(drivers, parameters, constants)

A struct to store drivers, parameters and constants. 
"""
struct Inputs
  drivers::Drivers
  parameters::Parameters
  constants::Constants
end

"""
  Output(name, range, unit)

A struct to store output name, range and unit.
"""
struct Output
  name
  range
  unit
end

"""
    mat(parameterisation::Function, inputs::Inputs, parameters, steps::Real)

Evaluate the function parameterisation on grid values of drivers x and y. 
Returns vectors x and y, and a matrix of function output at those points, FMatrix.
"""
function mat(parameterisation::Function, inputs::Inputs, parameters, steps::Real)
  # range from min to max for n steps (size = steps) 
  x = collect(range(inputs.drivers.ranges[1][1], length=steps, stop=inputs.drivers.ranges[1][2]))   
  y = collect(range(inputs.drivers.ranges[2][1], length=steps, stop=inputs.drivers.ranges[2][2])) 
  # Grid (size = steps*steps)
  X = repeat(1:steps, inner=steps)  
  Y = repeat(1:steps, outer=steps) 
  # Grid (size = steps*steps)
  X2 = repeat(x, inner=steps)    
  Y2 = repeat(y, outer=steps)
  # args for parameterisation
  drivers = [(X2[i], Y2[i]) for i in 1:steps*steps]
  parameters = repeat([parameters], steps*steps)
  constants = repeat([inputs.constants.values], steps*steps)
  # parameterisation output on Grid
  FMatrix = Matrix(sparse(X, Y, parameterisation.(drivers, parameters, constants)))
  return x, y, FMatrix
end

"""
    d1_vec(y, parameterisation::Function, inputs::Inputs, parameters, steps::Real)

Evaluate the function parameterisation on a line of driver x, with constant y. 
"""
function d1_vec(y, parameterisation::Function, inputs::Inputs, parameters, steps::Real) 
  x = collect(range(inputs.drivers.ranges[1][1], inputs.drivers.ranges[1][2], steps)) # min d1 to max d1, n steps
  y = repeat([y], steps)
  drivers = [(x[i], y[i]) for i in 1:steps] 
  parameters = repeat([parameters], steps)
  constants = repeat([inputs.constants.values], steps)
  vec = parameterisation.(drivers, parameters, constants)
  return vec
end

"""
    d2_vec(x, parameterisation::Function, inputs::Inputs, parameters, steps::Real)

Evaluate the function parameterisation on a line of driver y, with constant x.
"""
function d2_vec(x, parameterisation::Function, inputs::Inputs, parameters, steps::Real)
  y = collect(range(inputs.drivers.ranges[2][1], inputs.drivers.ranges[2][2], steps)) # min d2 to max d2, n steps
  x = repeat([x], steps)
  drivers = [(x[i], y[i]) for i in 1:steps] 
  parameters = repeat([parameters], steps)
  constants = repeat([inputs.constants.values], steps)
  vecM = parameterisation.(drivers, parameters, constants)
  return vecM
end

function parameterisation(drivers, parameters, constants)
  return parameterisation((drivers..., parameters..., constants...)...)
end


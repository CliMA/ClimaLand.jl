#=
using ParamViz
using ClimaLSM
using ClimaLSM.Canopy
FT = Float64
=#

function Beer_app_f()
  drivers = Drivers(("PAR (μmol m⁻² s⁻¹))", "LAI (m² m⁻²))"),
                    FT.((500 * 1e-6, 5)),
                    (FT.([0, 1500 * 1e-6]), FT.([0, 10])),
                    (1e6, 1.0) # scalers
                   )

  parameters = Parameters(("canopy reflectance, ρ_leaf",
                           "extinction coefficient, K",
                           "clumping index, Ω"),
                          (FT(0.1), FT(0.6), FT(0.69)),
                          (FT.([0, 1]), FT.([0, 1]), FT.([0, 1])),
                          (1, 1, 1) # scalers
                         )

  # need a method with no constant! 
  # hack: useless constants
  constants = Constants(("a", "b"), (FT(1), FT(2)))

  inputs = Inputs(drivers, parameters, constants)

  output = Output("APAR (μmol m⁻² s⁻¹)", [0, 1500 * 1e-6], 1e6)

  import ParamViz.parameterisation
  function parameterisation(PAR, LAI, ρ_leaf, K, Ω, a, b)   
    APAR = plant_absorbed_ppfd(PAR, ρ_leaf, K, LAI, Ω) 
    return APAR
  end

  beer_app = webapp(parameterisation, inputs, output)
  return beer_app
end

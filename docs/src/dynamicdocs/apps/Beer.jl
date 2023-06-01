#=
using ParamViz
using ClimaLSM
using ClimaLSM.Canopy
FT = Float64
=#

function Beer_app_f()
  drivers = Drivers(("PAR (μmol m⁻² s⁻¹))", "LAI (m² m⁻²))"),
                    FT.((500 * 1e-6, 5)),
                    (FT.([0, 1500 * 1e-6]), FT.([0, 10]))
                   )

  parameters = Parameters(("canopy reflectance, ρ_leaf",
                           "extinction coefficient, K",
                           "clumping index, Ω"),
                          (FT(0.1), FT(0.6), FT(0.69)),
                          (FT.([0, 1]), FT.([0, 1]), FT.([0, 1]))
                         )

  constants = Constants(("nothing1", "nothing2"), (FT(1), FT(2))) # just need something... could add method without

  inputs = Inputs(drivers, parameters, constants)

  output = Output("APAR (μmol m⁻² s⁻¹)", [0, 1500])

  function beer_APAR(PAR, LAI, ρ_leaf, K, Ω)   
    APAR = plant_absorbed_ppfd(PAR, ρ_leaf, K, LAI, Ω) * 1e6 # from mol to μmol
    return APAR
  end

  function beer_APAR(inputs)
    PAR, LAI = inputs.drivers.values[1], inputs.drivers.values[2] 
    ρ_leaf, K, Ω = [inputs.parameters.values[i] for i in 1:length(inputs.parameters.values)]
    return beer_APAR(PAR, LAI, ρ_leaf, K, Ω)
  end

  function beer_APAR(drivers, parameters, constants)
    PAR, LAI = drivers[1], drivers[2]
    ρ_leaf, K, Ω = [parameters[i] for i in 1:3]
    return beer_APAR(PAR, LAI, ρ_leaf, K, Ω)
  end

  beer_app = webapp(beer_APAR, inputs, output)
  return beer_app
end

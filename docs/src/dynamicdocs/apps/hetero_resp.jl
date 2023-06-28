#=
using ParamViz 
using ClimaLSM
using ClimaLSM.Soil.Biogeochemistry
FT = Float64
=#

function Rh_app_f()
  drivers = Drivers( #Drivers(names, values, ranges)
                    ("Soil temperature", "Soil moisture"), #drivers.names
                    FT.((283, 0.4)), #drivers.values
                    (FT.([273, 323]),FT.([0.1, 0.8])) #drivers.ranges
                   )

  parameters = Parameters(
                          (# names   
                           "Soil porosity (m³ m⁻³)",
                           "Pre-exponential factor (kg C m-3 s-1)",                            
                           "Activation energy (J mol-1)",                            
                           "Michaelis constant (kg C m-3)",                          
                           "Michaelis constant for O2 (m3 m-3)",                            
                           "Volumetric fraction of O₂ in the soil air, dimensionless",                                                      
                           "Fraction of soil carbon that is considered soluble, dimensionless",),                       
                          (# values
                           FT(0.556), # ν   
                           FT(194e3), # α_sx
                           FT(61e3), # Ea_sx
                           FT(5e-3), # kM_sx
                           FT(0.004), # kM_o2
                           FT(0.209), # O2_a
                           FT(0.024)), # p_sx
                          (#ranges
                           FT.([0.1, 0.9]), # porosity
                           FT.([100e3, 300e3]), # α_sx
                           FT.([30e3, 90e3]), # Ea_sx
                           FT.([1e-3, 9e-3]), # kM_sx
                           FT.([0.001, 0.009]), # kM_o2
                           FT.([0.005, 0.5]), # O2_a
                           FT.([0.005, 0.5])) # p_sx
                         )

  constants = Constants(
                      (# names   
                       "Pressure at the surface of the soil (Pa)",
                       "Diffusivity of soil C substrate in liquid (unitless)",
                       "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)",
                       "Diffusion coefficient for CO₂ in air at standard temperature and pressure (m² s⁻¹)",
                       "Absolute value of the slope of the line relating log(ψ) versus log(θ) (unitless)",
                       "Diffusion coefficient of oxygen in air, dimensionless"
                       ),                        
                       (# values
                       FT(101e3), # P_sfc
                       FT(0.1816), # θ_a100
                       FT(1.39e-5), # D_ref
                       FT(4.547), # b
                       FT(3.17), # D_liq
                       FT(1.67) # D_oa
                       )                     
                      )

  inputs = Inputs(drivers, parameters, constants)

  output = Output(
                "Rh (μmol m⁻² s⁻¹)", #name
                [0, 20] #range
                )

  function hetero_resp()

    return Rh 
  end

  function hetero_resp(inputs)

    return Rh 
  end

  function hetero_resp(drivers, parameters, constants)

    return Rh 
  end

  Rh_app = webapp(hetero_resp, inputs, output)
  return Rh_app
end

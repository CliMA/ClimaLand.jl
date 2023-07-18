#=
using ParamViz 
using ClimaLSM
using ClimaLSM.Soil.Biogeochemistry
FT = Float64

# this whole mess is because of parameters.jl ...
using CLIMAParameters
using Thermodynamics
using Insolation
using SurfaceFluxes
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
earth_param_set = create_lsm_parameters(FT)
=#

function Rh_app_f()
  drivers = Drivers( # Drivers(names, values, ranges)
                    ("Tₛ (°C)", "θ (m³ m⁻³)"), # drivers.names
                    FT.((283, 0.15)), # drivers.values
                    (FT.([273, 323]), FT.([0.0, 1.0])) # drivers.ranges #MIND THE UNIT FOR TEMPERATURE
                   )

  parameters = Parameters(
                          (# names   
                           "Soil porosity, ν (m³ m⁻³)",
                           "Pre-exponential factor, α_sx (kg C m⁻³ s⁻¹)",                            
                           "Activation energy, Ea_sx (J mol⁻¹)",                            
                           "Michaelis constant for soil, kM_sx (kg C m⁻³)",                          
                           "Michaelis constant for O₂, kM_O₂ (m³ m⁻³)",                            
                           "Volumetric fraction of O₂ in the soil air, O₂_a (dimensionless)",
                           "Fraction of soil carbon that is considered soluble, p_sx (dimensionless)",
                           "Soil organic C, Csom (kg C m⁻³)"
                          ),                       
                          (# values
                           FT(0.556), # ν   
                           FT(194e3), # α_sx
                           FT(61e3), # Ea_sx
                           FT(5e-3), # kM_sx
                           FT(0.004), # kM_o2
                           FT(0.209), # O2_a
                           FT(0.024), # p_sx
                           FT(5.0) # Csom
                          ), 
                          (# ranges
                           FT.([0.0, 1.0]), # porosity
                           FT.([100e3, 300e3]), # α_sx
                           FT.([50e3, 70e3]), # Ea_sx
                           FT.([1e-10, 0.1]), # kM_sx
                           FT.([1e-10, 0.1]), # kM_o2
                           FT.([0.005, 0.5]), # O2_a
                           FT.([0.005, 0.5]), # p_sx
                           FT.([1.0, 10.0]) # Csom
                          ) 
                         )

  constants = Constants(
                      (# names
                       "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)",
                       "Diffusivity of soil C substrate in liquid (unitless)",
                       "Diffusion coefficient of oxygen in air, dimensionless"
                      ),    
    
                       (# values
                       FT(0.1816), # θ_a100
                       FT(3.17), # D_liq
                       FT(1.67) # D_oa
                       )                     
                       )

  inputs = Inputs(drivers, parameters, constants)

  output = Output(
                "Rh (μmol m⁻² s⁻¹)", # name
                [0, 20] # range
                )

  function hetero_resp(Tₛ, θ, # drivers
                       ν, α_sx, Ea_sx, kM_sx, kM_O₂, O₂_a, p_sx, Csom, # parameters
                       θ_a100, D_liq, D_oa) # constants
    params = SoilCO2ModelParameters{FT}(;
                                        ν = ν,
                                        θ_a100 = θ_a100,
                                        D_liq = D_liq,                                
                                        α_sx = α_sx,
                                        Ea_sx = Ea_sx,
                                        kM_sx = kM_sx,
                                        kM_o2 = kM_O₂,
                                        O2_a = O₂_a,
                                        D_oa = D_oa,
                                        p_sx = p_sx,
                                        earth_param_set = earth_param_set,
                                       )
    # θ has to be lower than porosity
    if θ > ν
        θ = ν
    end

    Rh = microbe_source(Tₛ, θ, Csom, params) * 1e6 # from mol to μmol
    return Rh 
  end

  function hetero_resp(inputs)
    Tₛ, θ = inputs.drivers.values
    ν, α_sx, Ea_sx, kM_sx, kM_O₂, O₂_a, p_sx, Csom = inputs.parameters.values
    θ_a100, D_liq, D_oa = inputs.constants.values    
    return hetero_resp(Tₛ, θ, ν, α_sx, Ea_sx, kM_sx, kM_O₂, O₂_a, p_sx, Csom, θ_a100, D_liq, D_oa)  
  end

  function hetero_resp(drivers, parameters, constants)
    Tₛ, θ = drivers
    ν, α_sx, Ea_sx, kM_sx, kM_O₂, O₂_a, p_sx, Csom = parameters
    θ_a100, D_liq, D_oa = constants    
    return hetero_resp(Tₛ, θ, ν, α_sx, Ea_sx, kM_sx, kM_O₂, O₂_a, p_sx, Csom, θ_a100, D_liq, D_oa) 
  end

  Rh_app = webapp(hetero_resp, inputs, output)
  return Rh_app
end


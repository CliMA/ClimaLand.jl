# Test with ClimaLSM Farquhar An
# using local ClimaLSM 
using ParamViz

drivers = Drivers(("PAR (W m-2)", "T (K)"), (1, 1), ([-5, 5], [-5, 5]))

parameters = Parameters(("Moisture stress, β (unit)",
                         "Leaf area index, LAI (m² m⁻²)",
                         "CO2 concentration, ca (ppm)", # ppm?
                         "Vapor pressure deficit, VPD (Pa)"),
                        (1.0, # β
                         1.0, # LAI
                         400, # ca
                         100, # VPD, Pa
                         ),
                        ([-5, 5], # β
                         [-5, 5], # LAI
                         [300, 500], # ca
                         [50, 200], # VPD
                         ))

constants = Constants(("θs", # Sun zenith angle
                       "ld", # Leaf angle distribution
                       "p_leaf", # PAR canopy reflectance
                       "Ω", # Clumping index
                       "Rstar25", # co2 compensation at 25c
                       "ΔHΓstar", # a constant energy of activation
                       "To", # a standard temperature 
                       "R", # the universal gas constant
                       "Vcmax25", # the maximum rate of carboxylation of Rubisco
                       "ΔHJmax", # a constant
                       "θj", # an empirical "curvature parameter"
                       "ϕ", # the quantum yield of photosystem II
                       "ΔHVcmax", # a constant
                       "g1", # a constant
                       "Kc25", # a constant
                       "ΔHkc", # a constant
                       "Ko25", # a constant
                       "ΔHko", # a constant
                       "oi", # an empirical parameter
                       "f", # an empirical factor
                        ),
                      (1.0, # θs
                       0.5, # ld
                       0.1, # p_leaf
                       1, # Ω
                       4.275e-5, # Rstar25
                       37830, # ΔHΓstar
                       298.15, # To
                       8.314, # R, J/mol
                       5e-5, # Vcmax25
                       43540, # ΔHJmax
                       0.9, # θj
                       0.6, # ϕ
                       58520, # ΔHVcmax
                       141, # g1
                       4.049e-4, # Kc25
                       79430, # ΔHkc
                       0.2874, # Ko25
                       36380, # ΔHko
                       0.209, # oi
                       0.015, # f
                      )) 

inputs = Inputs(drivers, parameters, constants)

output = Output("output", [-12, 12])

function leaf_photosynthesis(PAR, T, β, LAI, θs, ld, p_leaf, Ω, )   

  mechanism = Canopy.C3()

  K = extinction_coeff(ld, θs)
  APAR = plant_absorbed_ppfd(PAR, ρ_leaf, K, LAI, Ω)
  Jmax = max_electron_transport(Vcmax25, ΔHJmax, T, To, R)
  J = electron_transport(APAR, Jmax, θj, ϕ)
  Vcmax = compute_Vcmax(Vcmax25, T, To, R, ΔHVcmax)
  Γstar = co2_compensation(Γstar25, ΔHΓstar, T, To, R)
  medlynt = medlyn_term(g1, VPD)
  ci = intercellular_co2(ca, Γstar, medlynt)
  Aj = light_assimilation(Ref(mechanism), J, ci, Γstar)
  Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
  Ko = MM_Ko(Ko25, ΔHko, T, To, R)
  Ac = rubisco_assimilation(Ref(mechanism), Vcmax, ci, Γstar, Kc, Ko, oi)
  Rd = dark_respiration(Vcmax25, β, f, ΔHRd, T, To, R)
  An = net_photosynthesis(Ac, Aj, Rd, β)
  return An
end














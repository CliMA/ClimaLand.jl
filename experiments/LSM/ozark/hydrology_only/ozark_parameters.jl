# Soil parameters
soil_ν = FT(0.45) # m3/m3, Natan's value: 0.55, Wang et al value: 0.45
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_S_s = FT(1e-3) # 1/m, guess
θ_r = FT(0.067) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
# With these parameters, the soil water potential bottoms out at ~0.1 during the second half of the year.
#soil_vg_n = FT(2.6) # unitless. We cited Yujie's paper for these, but I can't find them now
#soil_vg_α = FT(1.3) # inverse meters.

# With these parameters, we get only a small dip still, soil water potential still doesnt look correct
#soil_vg_n = FT(1.5) # unitless, Natan: 1.5
#soil_vg_α = FT(0.1089) # inverse meters. Natan

# With these parameters, we get the dip in leaf water potential, but it's small amplitude.
soil_vg_α = FT(2.6) # unitless.
soil_vg_n = FT(1.3) # inverse meters.

# Beer Lambert model parameters
Ω = FT(0.69)
ld = FT(0.5)
ρ_leaf = FT(0.1)
λ_γ = FT(5e-7)

# Conductance Model
g1 = FT(141) # Natan used 300; in Yujie's paper: 141 sqrt(Pa) for Medlyn model
Drel = FT(1.6)
g0 = FT(1e-4) # Natan used zero

#Photosynthesis model
oi = FT(0.209)
ϕ = FT(0.6)
θj = FT(0.9)
f = FT(0.015)
sc = FT(5e-6) # Bonan's book: range of 2-5e-6
pc = FT(-2e6) # Bonan's book
Vcmax25 = FT(4.5e-5) # from Yujie's paper, Natan used 9e-5
Γstar25 = FT(4.275e-5)
Kc25 = FT(4.049e-4)
Ko25 = FT(0.2874)
To = FT(298.15)
ΔHkc = FT(79430)
ΔHko = FT(36380)
ΔHVcmax = FT(58520)
ΔHΓstar = FT(37830)
ΔHJmax = FT(43540)
ΔHRd = FT(46390)

# Plant Hydraulics and general plant parameters
SAI = FT(0.00242) # m2/m2, estimated from Wang et al. 
LAI = FT(4.2) # m2/m2, from Wang et al.
f_root_to_shoot = FT(3.5) # from google?
RAI = (SAI + LAI) * f_root_to_shoot # CLM
K_sat_plant = 1.8e-8 # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
plant_ν = FT(0.7) # guess, m3/m3
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(0.5) # Agrees with Natan.
z0_m = FT(2)
z0_b = FT(0.2)

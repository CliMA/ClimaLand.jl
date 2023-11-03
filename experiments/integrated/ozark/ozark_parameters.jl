## Ozark parameters taken from the following:
## Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
## Holtzman, Natan; private communication
## CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
## # Bonan, G. Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.

## LAI data from MODIS
# Heterotrophic respiration parameters
θ_a100 = FT(0.1816)
D_ref = FT(1.39e-5)
b = FT(4.547)
D_liq = FT(3.17)
# DAMM
α_sx = FT(194e3)
Ea_sx = FT(61e3)
kM_sx = FT(5e-3)
kM_o2 = FT(0.004)
O2_a = FT(0.209)
D_oa = FT(1.67)
p_sx = FT(0.024)

# Autotrophic respiration parameters
ne = FT(8 * 1e-4)
ηsl = FT(0.01)
σl = FT(0.05)
μr = FT(1.0)
μs = FT(0.1)
f1 = FT(0.012)
f2 = FT(0.25)

# Soil parameters
soil_ν = FT(0.5) # m3/m3
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
θ_r = FT(0.067) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
κ_quartz = FT(7.7) # W/m/K
κ_minerals = FT(2.5) # W/m/K
κ_om = FT(0.25) # W/m/K
κ_liq = FT(0.57) # W/m/K
κ_ice = FT(2.29) # W/m/K
κ_air = FT(0.025); #W/m/K
ρp = FT(2700); # kg/m^3
κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
κ_dry = Soil.κ_dry(ρp, soil_ν, κ_solid, κ_air)
κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, soil_ν, κ_ice)
κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, soil_ν, κ_liq);
ρc_ds = FT((1 - soil_ν) * 4e6); # J/m^3/K
z_0m_soil = FT(0.01)
z_0b_soil = FT(0.001)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.3)

# TwoStreamModel parameters
Ω = FT(0.89)
ld = FT(0.5)
α_PAR_leaf = FT(0.1)
λ_γ_PAR = FT(5e-7)
λ_γ_NIR = FT(1.65e-6)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.45)
τ_NIR_leaf = FT(0.25)
n_layers = UInt64(20)
ϵ_canopy = FT(0.97)

# Energy Balance model
ac_canopy = FT(2.5e3)

# Conductance Model
g1 = FT(70) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
Drel = FT(1.6)
g0 = FT(1e-4)

#Photosynthesis model
oi = FT(0.209)
ϕ = FT(0.6)
θj = FT(0.9)
f = FT(0.015)
sc = FT(2e-6) # Bonan's book: range of 2-5e-6
pc = FT(-2e6) # Bonan's book: -2e6
Vcmax25 = FT(9e-5) # from Yujie's paper 4.5e-5 , Natan used 9e-5
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
SAI = FT(1.0) # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
maxLAI = FT(4.2) # m2/m2, from Wang et al.
f_root_to_shoot = FT(3.5)
RAI = (SAI + maxLAI) * f_root_to_shoot # CLM
K_sat_plant = 5e-9 # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.08 * 0.0098) # Natan used 0.05
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
capacity = FT(10) # kg/m^2
plant_ν = capacity / (maxLAI / 2 * h_leaf + SAI * h_stem) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(0.5) # from Natan
z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m

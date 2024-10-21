"""Site-specific model parameters for running Clima Land on the Niwot Ridge
fluxtower site."""

# Data download link
data_link = "https://caltech.box.com/shared/static/r6gvldgabk3mvtx53gevnlaq1ztsk41i.csv"

# Timezone (offset from UTC in hrs)
time_offset = 7

# Site latitude and longitude
lat = FT(40.0329) # degree
long = FT(-105.5464) # degree

# Height of sensor
atmos_h = FT(21.5)
# Metzger, Stefan & Burba, George & Burns, Sean & Blanken, Peter & Li, 
# Jiahong & Luo, Hongyan & Zulueta, Rommel. (2016). Optimization of an enclosed
# gas analyzer sampling system for measuring eddy covariance fluxes of H2O and
# CO2. Atmospheric Measurement Techniques. 9. 1341-1359. 10.5194/amt-9-1341-2016. 

# Soil parameters
soil_ν = FT(0.45) # m3/m3
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
θ_r = FT(0.067) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021

# Soil makeup
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.2)

# TwoStreamModel parameters
Ω = FT(0.71)
ld = FT(0.5)
G_Function = ConstantGFunction(ld)
α_PAR_leaf = FT(0.1)
λ_γ_PAR = FT(5e-7)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.35)
τ_NIR_leaf = FT(0.25)
ϵ_canopy = FT(0.97)

# Energy Balance model
ac_canopy = FT(3e3)

# Conductance Model
g1 = FT(141) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
Drel = FT(1.6)
g0 = FT(1e-4)

#Photosynthesis model
Vcmax25 = FT(9e-5) # from Yujie's paper 4.5e-5 , Natan used 9e-5

# Plant Hydraulics and general plant parameters
SAI = FT(1.0) # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
f_root_to_shoot = FT(3.5)
K_sat_plant = 5e-9 # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
plant_ν = FT(8.06e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(1.0)
z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m

## Some site parameters have been taken from
## Ma, S., Baldocchi, D. D., Xu, L., Hehn, T. (2007)
## Inter-Annual Variability In Carbon Dioxide Exchange Of An
## Oak/Grass Savanna And Open Grassland In California, Agricultural
## And Forest Meteorology, 147(3-4), 157-171. https://doi.org/10.1016/j.agrformet.2007.07.008 
## CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
# Bonan, G. Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.

# Autotrophic respiration parameters
ne = FT(8 * 1e-4)
ηsl = FT(0.01)
σl = FT(0.05)
μr = FT(1.0)
μs = FT(0.1)
f1 = FT(0.012)
f2 = FT(0.25)

# Soil parameters
soil_ν = FT(0.45) # m3/m3
soil_K_sat = FT(0.45 / 3600 / 100) # m/s,
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(2.0) # unitless
soil_vg_α = FT(2.0) # inverse meters
θ_r = FT(0.067) # m3/m3, 

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.38)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.1);
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
soil_α_PAR = FT(0.3)
soil_α_NIR = FT(0.4)

# TwoStreamModel parameters
Ω = FT(1.0)
ld = FT(0.5)
α_PAR_leaf = FT(0.11)
λ_γ_PAR = FT(5e-7)
λ_γ_NIR = FT(1.65e-6)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.35)
τ_NIR_leaf = FT(0.34)
n_layers = UInt64(20)
ϵ_canopy = FT(0.97)

# Conductance Model
g1 = FT(166) # CLM C3 grass
Drel = FT(1.6)
g0 = FT(1e-4)

#Photosynthesis model
oi = FT(0.209)
ϕ = FT(0.6)
θj = FT(0.9)
f = FT(0.015)
sc = FT(2e-6) # Bonan's book: range of 2-5e-6
pc = FT(-2e6) # Bonan's book: -2e6
Vcmax25 = FT(4.225e-5) # CLM C3 grass, Slevin et al. 2015
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

# Energy Balance model
ac_canopy = FT(745)

# Plant Hydraulics and general plant parameters
maxLAI = FT(maximum(LAI_timeseries))
SAI = FT(0)
f_root_to_shoot = FT(1.0)
RAI = maxLAI * f_root_to_shoot # CLM
K_sat_plant = 5e-9 # m/s # seems much too small?
ψ63 = FT(-2.7 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
capacity = FT(2.0) # kg/m^2
plant_ν = capacity / (maxLAI * h_leaf) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(2.6) # from Bonan Table 2.3
z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m

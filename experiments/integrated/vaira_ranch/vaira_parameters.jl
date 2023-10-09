# Autotrophic respiration parameters
ne = FT(8 * 1e-4)
ηsl = FT(0.01)
σl = FT(0.05)
μr = FT(1.0)
μs = FT(0.1)
f1 = FT(0.012)
f2 = FT(0.25)

# Soil parameters - Bonan silt loam
soil_ν = FT(0.485) # m3/m3
soil_K_sat = FT(2.59 / 3600 / 100) # m/s,
soil_S_s = FT(1e-3) # 1/m, guess
soil_bc_c = FT(0.188) # unitless
soil_bc_ψb = FT(-0.786) # inverse meters
θ_r = FT(0.0) # m3/m3, 

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.295)
ν_ss_om = FT(0.0)
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
z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.2)

# TwoStreamModel parameters
Ω = FT(0.8)
ld = FT(0.5)
α_PAR_leaf = FT(0.1)
λ_γ_PAR = FT(5e-7)
λ_γ_NIR = FT(1.65e-6)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.45)
τ_NIR_leaf = FT(0.25)
n_layers = UInt64(20)
diff_perc = FT(0.2)
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
Vcmax25 = FT(5e-5) # ????
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
capacity = FT(22.0) # kg/m^2
plant_ν = capacity / (maxLAI * h_leaf) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(2.6) # from Bonan Table 2.3
z0_m = FT(0.25)
z0_b = FT(0.25)

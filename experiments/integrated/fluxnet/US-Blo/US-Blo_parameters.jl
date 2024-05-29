"""Site-specific model parameters for running Clima Land on the Vaira Ranch
fluxtower site."""

## Some site parameters have been taken from
## Ma, S., Baldocchi, D. D., Xu, L., Hehn, T. (2007)
## Inter-Annual Variability In Carbon Dioxide Exchange Of An
## Oak/Grass Savanna And Open Grassland In California, Agricultural
## And Forest Meteorology, 147(3-4), 157-171. https://doi.org/10.1016/j.agrformet.2007.07.008 
## CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
# Bonan, G. Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.

# Data download link
#data_link = "https://caltech.box.com/shared/static/54huhy74kxnn23i6w1s54vddo2j4hl71.csv"
#LAI_link = "/Users/mitraasadollahi/Projects/CliMA/2_calib_physical_sims/Data/LAI/$(site_ID)_LAI.csv"#"https://caltech.box.com/shared/static/y5vf8s9qkoogglc1bc2eyu1k95sbjsc3.csv"

# Height of sensor
atmos_h = FT(21.5)

# Timezone (offset from UTC in hrs)
time_offset = 8

# Site latitude and longitude
lat = FT(38.88) # degree
long = FT(-120.63) # degree

# Soil parameters
soil_ν = FT(0.4) # m3/m3
soil_K_sat = FT(2.5/ 3600 / 100) # m/s,X unit is cm/hr where X/3600/100 = m/s
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(1.5) # unitless
soil_vg_α = FT(0.2) # inverse meters
θ_r = FT(0.07) # m3/m3, 

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.2)


# TwoStreamModel parameters
Ω = FT(0.74)
ld = FT(0.5)
G_Function = ConstantGFunction(ld)
α_PAR_leaf = FT(0.7)
λ_γ_PAR = FT(5e-7)
λ_γ_NIR = FT(1.65e-6)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.35)
τ_NIR_leaf = FT(0.1)
ϵ_canopy = FT(0.97)



# Conductance Model
g1 = FT(80) # CLM C3 grass
Drel = FT(1.6)
g0 = FT(1e-4)

#Photosynthesis model
Vcmax25 = FT(4.5e-5) # CLM C3 grass, Slevin et al. 2015

# Energy Balance model
ac_canopy = FT(2792)

# Plant Hydraulics and general plant parameters
SAI = FT(1.0)
f_root_to_shoot = FT(3.5)
K_sat_plant = 5e-9 # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
plant_ν = FT(2e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(10) # from Bonan Table 2.3
z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m

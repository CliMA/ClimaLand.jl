## Some site parameters have been taken from
## Ma, S., Baldocchi, D. D., Xu, L., Hehn, T. (2007)
## Inter-Annual Variability In Carbon Dioxide Exchange Of An
## Oak/Grass Savanna And Open Grassland In California, Agricultural
## And Forest Meteorology, 147(3-4), 157-171. https://doi.org/10.1016/j.agrformet.2007.07.008 
## CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
# Bonan, G. Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.


# Soil parameters
soil_ν = FT(0.45) # m3/m3
soil_K_sat = FT(0.45 / 3600 / 100) # m/s,
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(2.0) # unitless
soil_vg_α = FT(2.0) # inverse meters
θ_r = FT(0.067) # m3/m3, 

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.2)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.4);
κ_quartz = FT(7.7) # W/m/K
κ_minerals = FT(2.5) # W/m/K
κ_om = FT(0.25) # W/m/K
κ_liq = FT(0.57) # W/m/K
κ_ice = FT(2.29) # W/m/K
κ_air = FT(0.025); #W/m/K
ρp = FT(2700); # kg/m^3
κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
ρc_ds = FT((1 - soil_ν) * 4e6); # J/m^3/K
z_0m_soil = FT(0.01)
z_0b_soil = FT(0.001)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.3)
soil_α_NIR = FT(0.4)

(; ν_ss_om, ν_ss_quartz, ν_ss_gravel) =
    ClimaLand.Soil.soil_composition_parameters(subsurface_space, FT)
(; ν, hydrology_cm, K_sat, θ_r) =
    ClimaLand.Soil.soil_vangenuchten_parameters(subsurface_space, FT)
soil_albedo = Soil.CLMTwoBandSoilAlbedo{FT}(;
    ClimaLand.Soil.clm_soil_albedo_parameters(surface_space)...,
)
S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)
soil_params = Soil.EnergyHydrologyParameters(
    FT;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm,
    K_sat,
    S_s,
    θ_r,
    albedo = soil_albedo,
)
f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4 / 1000) # m/s
runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
    f_over = f_over,
    f_max = ClimaLand.Soil.topmodel_fmax(surface_space, FT),
    R_sb = R_sb,
)


# Spatially varying canopy parameters from CLM
g1 = ClimaLand.Canopy.clm_medlyn_g1(surface_space)
rooting_depth = ClimaLand.Canopy.clm_rooting_depth(surface_space)
(; is_c3, Vcmax25) =
    ClimaLand.Canopy.clm_photosynthesis_parameters(surface_space)
(; Ω, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf) =
    ClimaLand.Canopy.clm_canopy_radiation_parameters(surface_space)
# Energy Balance model
ac_canopy = FT(2.5e3)
# Plant Hydraulics and general plant parameters
SAI = FT(0.0) # m2/m2
f_root_to_shoot = FT(3.5)
RAI = FT(1.0)
K_sat_plant = FT(5e-9) # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
n_stem = 0
n_leaf = 1
h_stem = FT(0.0)
h_leaf = FT(1.0)
zmax = FT(0.0)
h_canopy = h_stem + h_leaf
compartment_midpoints =
    n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
compartment_surfaces = n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m

soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)

module FluxnetSimulations

using ClimaTimeSteppers
using ClimaComms
import ClimaComms: context, device
using SciMLBase
using Dates
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics
using ClimaLand
const FT = Float64

include("initial_conditions.jl")

import ..Diagnostics: close_output_writers

###################################
#     SIMULATION MODULE CODE      #
###################################

"""
    get_parameters(::Val{:US_Ha1}; kwargs...)

Gets parameters for the Fluxnet site US-Ha1 (Massachusetts Harvard Forest)
and returns them as a Named Tuple. Default parameters are provided and
can be overriden using keyword arguments.

Data sources:

Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
Xu and Baldocchi, 2003
https://atmos.seas.harvard.edu/research-harvard_forest-instrumentation

"""
function get_parameters(::Val{:US_Ha1};
    time_offset = 5,
    lat = FT(42.5378), # degree
    long = FT(-72.1715), # degree
    atmos_h = FT(30),
    soil_ν = FT(0.5), # m3/m3
    soil_K_sat = FT(4e-7), # m/s, matches Natan
    soil_S_s = FT(1e-3), # 1/m, guess
    soil_vg_n = FT(2.05), # unitless
    soil_vg_α = FT(0.04), # inverse meters
    θ_r = FT(0.067), # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.001),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),
    Ω = FT(0.69),
    ld = FT(0.5),
    G_Function = ConstantGFunction(ld),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),
    ac_canopy = FT(2.5e3),
    g1 = FT(141), # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
    Drel = FT(1.6),
    g0 = FT(1e-4),
    Vcmax25 = FT(9e-5), # from Yujie's paper 4.5e-5 , Natan used 9e-5
    SAI = FT(1.0), # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 5e-9, # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098), # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4), # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098), # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(2.46e-4),
    plant_S_s = FT(1e-2 * 0.0098), # m3/m3/MPa to m3/m3/m
    rooting_depth = FT(0.5), # from Natan
    n_stem = Int64(0),
    n_leaf = Int64(1),
    h_leaf = FT(0.5), # m, Xu and Baldocchi, 2003
    h_stem = FT(0), # m
    h_canopy = h_leaf + h_stem,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m
)

    return (; time_offset, lat, long, atmos_h, soil_ν, soil_K_sat,
        soil_S_s, soil_vg_n, soil_vg_α, θ_r, ν_ss_quartz, ν_ss_om,
        ν_ss_gravel, z_0m_soil, z_0b_soil, soil_ϵ, soil_α_PAR, soil_α_NIR,
        Ω, ld, G_Function, α_PAR_leaf, λ_γ_PAR, τ_PAR_leaf,
        α_NIR_leaf, τ_NIR_leaf, ϵ_canopy, ac_canopy, g1, Drel,
        g0, Vcmax25, SAI, f_root_to_shoot, K_sat_plant, ψ63,
        Weibull_param, a, conductivity_model, retention_model,
        plant_ν, plant_S_s, rooting_depth, n_stem, n_leaf,
        h_leaf, h_stem, h_canopy, z0_m, z0_b)

end


"""
    get_parameters(::Val{:US_Ha1}; kwargs...)

Gets parameters for the Fluxnet site US-Ha1 (Massachusetts Harvard Forest)
and returns them as a Named Tuple.

Data sources:

Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
Xu and Baldocchi, 2003
https://atmos.seas.harvard.edu/research-harvard_forest-instrumentation

"""
function get_parameters(::Val{:US_MOz};
    # Timezone (offset from UTC in hrs)
    time_offset = 7,

    # Height of sensor on flux tower
    atmos_h = FT(32),

    # Site latitude and longitude
    lat = FT(38.7441), # degree
    long = FT(-92.2000), # degree

    # Soil parameters
    soil_ν = FT(0.55), # m3/m3
    soil_K_sat = FT(4e-7), # m/s
    soil_S_s = FT(1e-2), # 1/m, guess
    soil_vg_n = FT(2.0), # unitless
    soil_vg_α = FT(0.05), # inverse meters
    θ_r = FT(0.04), # m3/m3, 

    # Soil makeup
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),

    # TwoStreamModel parameters
    Ω = FT(0.69),
    χl = FT(0.1),
    G_Function = CLMGFunction(χl),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),

    # Energy Balance model
    ac_canopy = FT(5e2),

    # Conductance Model
    g1 = FT(141), # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
    Drel = FT(1.6),
    g0 = FT(1e-4),

    #Photosynthesis model
    Vcmax25 = FT(6e-5), # from Yujie's paper 4.5e-5

    # Plant Hydraulics and general plant parameters
    pc = FT(-2.0e6),
    sc = FT(5e-6),
    SAI = FT(1.0), # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 7e-8, # m/s 
    ψ63 = FT(-4 / 0.0098), # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4), # unitless, Holtzman's original c param value
    a = FT(0.1 * 0.0098), # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(1.44e-4),
    plant_S_s = FT(1e-2 * 0.0098), # m3/m3/MPa to m3/m3/m
    rooting_depth = FT(0.5), # from Natan
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_stem = FT(9), # m
    h_leaf = FT(9.5), # m
    h_canopy = h_stem + h_leaf,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m
)
    return (; time_offset, lat, long, atmos_h, soil_ν, soil_K_sat,
        soil_S_s, soil_vg_n, soil_vg_α, θ_r, ν_ss_quartz, ν_ss_om,
        ν_ss_gravel, z_0m_soil, z_0b_soil, soil_ϵ, soil_α_PAR, soil_α_NIR,
        Ω, ld, G_Function, α_PAR_leaf, λ_γ_PAR, τ_PAR_leaf,
        α_NIR_leaf, τ_NIR_leaf, ϵ_canopy, ac_canopy, g1, Drel,
        g0, Vcmax25, SAI, f_root_to_shoot, K_sat_plant, ψ63,
        Weibull_param, a, conductivity_model, retention_model,
        plant_ν, plant_S_s, rooting_depth, n_stem, n_leaf,
        h_leaf, h_stem, h_canopy, z0_m, z0_b)
end

function get_parameters(::Val{:US_NR1};
    # Timezone (offset from UTC in hrs)
    time_offset = 7,

    # Site latitude and longitude
    lat = FT(40.0329), # degree
    long = FT(-105.5464), # degree

    # Height of sensor
    atmos_h = FT(21.5),
    # Metzger, Stefan & Burba, George & Burns, Sean & Blanken, Peter & Li, 
    # Jiahong & Luo, Hongyan & Zulueta, Rommel. (2016). Optimization of an enclosed
    # gas analyzer sampling system for measuring eddy covariance fluxes of H2O and
    # CO2. Atmospheric Measurement Techniques. 9. 1341-1359. 10.5194/amt-9-1341-2016. 

    # Soil parameters
    soil_ν = FT(0.45), # m3/m3
    soil_K_sat = FT(4e-7), # m/s,
    soil_S_s = FT(1e-3), # 1/m, guess
    soil_vg_n = FT(2.05), # unitless
    soil_vg_α = FT(0.04), # inverse meters
    θ_r = FT(0.0), # m3/m3, 

    # Soil makeup
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.1),
    z_0b_soil = FT(0.1),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),

    # TwoStreamModel parameters
    Ω = FT(0.71),
    ld = FT(0.5),
    G_Function = ConstantGFunction(ld),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.35),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),

    # Energy Balance model
    ac_canopy = FT(3e3),

    # Conductance Model
    g1 = FT(141), # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
    Drel = FT(1.6),
    g0 = FT(1e-4),

    #Photosynthesis model
    Vcmax25 = FT(9e-5), # from Yujie's paper 4.5e-5 , Natan used 9e-5

    # Plant Hydraulics and general plant parameters
    SAI = FT(1.0), # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 5e-9, # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098), # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4), # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098), # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(8.06e-4),
    plant_S_s = FT(1e-2 * 0.0098), # m3/m3/MPa to m3/m3/m
    rooting_depth = FT(1.0),
    n_stem = Int64(1),
    h_leaf = FT(6.5), # m
    h_stem = FT(7.5), # m
    h_canopy = h_leaf + h_stem,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m
)

    return (; time_offset, lat, long, atmos_h, soil_ν, soil_K_sat,
        soil_S_s, soil_vg_n, soil_vg_α, θ_r, ν_ss_quartz, ν_ss_om,
        ν_ss_gravel, z_0m_soil, z_0b_soil, soil_ϵ, soil_α_PAR, soil_α_NIR,
        Ω, ld, G_Function, α_PAR_leaf, λ_γ_PAR, τ_PAR_leaf,
        α_NIR_leaf, τ_NIR_leaf, ϵ_canopy, ac_canopy, g1, Drel,
        g0, Vcmax25, SAI, f_root_to_shoot, K_sat_plant, ψ63,
        Weibull_param, a, conductivity_model, retention_model,
        plant_ν, plant_S_s, rooting_depth, n_stem, n_leaf,
        h_leaf, h_stem, h_canopy, z0_m, z0_b)

end


"""


Some site parameters have been taken from
Ma, S., Baldocchi, D. D., Xu, L., Hehn, T. (2007)
Inter-Annual Variability In Carbon Dioxide Exchange Of An
Oak/Grass Savanna And Open Grassland In California, Agricultural
And Forest Meteorology, 147(3-4), 157-171. https://doi.org/10.1016/j.agrformet.2007.07.008 
CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
Bonan, G. Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.
"""
function get_parameters(::Val{:US_Var};
    # Timezone (offset from UTC in hrs)
    time_offset = 8,

    # Site latitude and longitude
    lat = FT(38.4133), # degree
    long = FT(-120.9508), # degree

    # Height of sensor
    atmos_h = FT(2), # from BADM

    # Soil parameters
    soil_ν = FT(0.45), # m3/m3, Bonan Table 8.3
    soil_K_sat = FT(0.45 / 3600 / 100), # m/s, Bonan Table 8.3
    soil_S_s = FT(1e-3), # 1/m, guess
    soil_vg_n = FT(2.0), # unitless, near Bonan Table 8.3
    soil_vg_α = FT(2.0), # inverse meters, Bonan Table 8.3
    θ_r = FT(0.067), # m3/m3, near Bonan Table 8.3

    # Soil heat transfer parameters; not needed for hydrology only test
    ν_ss_quartz = FT(0.3), # Xu and Baldocchi (2003)
    ν_ss_om = FT(0.02), # Xu and Baldocchi (2003)
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),

    # TwoStreamModel parameters
    Ω = FT(1.0),
    ld = FT(0.5),
    G_Function = ConstantGFunction(ld),
    α_PAR_leaf = FT(0.11),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.35),
    τ_NIR_leaf = FT(0.34),
    ϵ_canopy = FT(0.97),

    # Conductance Model
    g1 = FT(166), # CLM C3 grass
    Drel = FT(1.6),
    g0 = FT(1e-4),

    #Photosynthesis model
    Vcmax25 = FT(2 * 4.225e-5), # 2x CLM C3 grass, Slevin et al. 2015

    # Energy Balance model
    ac_canopy = FT(745),

    # Plant Hydraulics and general plant parameters
    pc = FT(-3e5),
    sc = FT(1e-3),
    SAI = FT(0),
    f_root_to_shoot = FT(1.0),
    K_sat_plant = 2e-8, # m/s
    ψ63 = FT(-2.7 / 0.0098), # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4), # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098), # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(8.93e-3),
    plant_S_s = FT(1e-2 * 0.0098), # m3/m3/MPa to m3/m3/m
    rooting_depth = FT(0.3), # based off of soil depth, Xu and Baldocchi
    n_stem = Int64(0),
    n_leaf = Int64(1),
    h_leaf = FT(0.5), # m, Xu and Baldocchi, 2003
    h_stem = FT(0), # m
    h_canopy = h_leaf + h_stem,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m
)
    return (; time_offset, lat, long, atmos_h, soil_ν, soil_K_sat,
        soil_S_s, soil_vg_n, soil_vg_α, θ_r, ν_ss_quartz, ν_ss_om,
        ν_ss_gravel, z_0m_soil, z_0b_soil, soil_ϵ, soil_α_PAR, soil_α_NIR,
        Ω, ld, G_Function, α_PAR_leaf, λ_γ_PAR, τ_PAR_leaf,
        α_NIR_leaf, τ_NIR_leaf, ϵ_canopy, ac_canopy, g1, Drel,
        g0, Vcmax25, SAI, f_root_to_shoot, K_sat_plant, ψ63,
        Weibull_param, a, conductivity_model, retention_model,
        plant_ν, plant_S_s, rooting_depth, n_stem, n_leaf,
        h_leaf, h_stem, h_canopy, z0_m, z0_b)
end

# fallback function, get generic values (decide later)
function get_parameters(site_ID;
    # Timezone (offset from UTC in hrs)
    time_offset = 7,

    # Height of sensor on flux tower
    atmos_h = FT(32),

    # Site latitude and longitude
    lat = FT(38.7441), # degree
    long = FT(-92.2000), # degree

    # Soil parameters
    soil_ν = FT(0.55), # m3/m3
    soil_K_sat = FT(4e-7), # m/s
    soil_S_s = FT(1e-2), # 1/m, guess
    soil_vg_n = FT(2.0), # unitless
    soil_vg_α = FT(0.05), # inverse meters
    θ_r = FT(0.04), # m3/m3, 

    # Soil makeup
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    ν_ss_gravel = FT(0.0),
    z_0m_soil = FT(0.01),
    z_0b_soil = FT(0.01),
    soil_ϵ = FT(0.98),
    soil_α_PAR = FT(0.2),
    soil_α_NIR = FT(0.2),

    # TwoStreamModel parameters
    Ω = FT(0.69),
    χl = FT(0.1),
    G_Function = CLMGFunction(χl),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.97),

    # Energy Balance model
    ac_canopy = FT(5e2),

    # Conductance Model
    g1 = FT(141), # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
    Drel = FT(1.6),
    g0 = FT(1e-4),

    #Photosynthesis model
    Vcmax25 = FT(6e-5), # from Yujie's paper 4.5e-5

    # Plant Hydraulics and general plant parameters
    pc = FT(-2.0e6),
    sc = FT(5e-6),
    SAI = FT(1.0), # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 7e-8, # m/s 
    ψ63 = FT(-4 / 0.0098), # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4), # unitless, Holtzman's original c param value
    a = FT(0.1 * 0.0098), # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param),
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a),
    plant_ν = FT(1.44e-4),
    plant_S_s = FT(1e-2 * 0.0098), # m3/m3/MPa to m3/m3/m
    rooting_depth = FT(0.5), # from Natan
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_stem = FT(9), # m
    h_leaf = FT(9.5), # m
    h_canopy = h_stem + h_leaf,
    z0_m = FT(0.13) * h_canopy,
    z0_b = FT(0.1) * z0_m
)
    return (; time_offset, lat, long, atmos_h, soil_ν, soil_K_sat,
        soil_S_s, soil_vg_n, soil_vg_α, θ_r, ν_ss_quartz, ν_ss_om,
        ν_ss_gravel, z_0m_soil, z_0b_soil, soil_ϵ, soil_α_PAR, soil_α_NIR,
        Ω, ld, G_Function, α_PAR_leaf, λ_γ_PAR, τ_PAR_leaf,
        α_NIR_leaf, τ_NIR_leaf, ϵ_canopy, ac_canopy, g1, Drel,
        g0, Vcmax25, SAI, f_root_to_shoot, K_sat_plant, ψ63,
        Weibull_param, a, conductivity_model, retention_model,
        plant_ν, plant_S_s, rooting_depth, n_stem, n_leaf,
        h_leaf, h_stem, h_canopy, z0_m, z0_b)
end

"""
    get_domain_info(::Val{:US_Ha1}; dz_bottom = FT(1.5), dz_top = FT(0.025),
        nelements = 20, zmin = FT(-10), zmax = FT(0)

Gets and returns primary domain information for the US-Ha1 (Massachusetts
Harvard Forest) Fluxnet site.

The data source comes from: Unknown.
"""
# for US-Ha1
function get_domain_info(::Val{:US_Ha1};
    dz_bottom = FT(1.5),
    dz_top = FT(0.025),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0)
)

    dz_tuple = (dz_bottom, dz_top),

    return (dz_tuple=dz_tuple,
        nelements=nelements,
        zmin=zmin,
        zmax=zmax)
end

"""
    get_domain_info(::Val{:US_MOz}; dz_bottom = FT(1.5), dz_top = FT(0.1),
        nelements = 20, zmin = FT(-10), zmax = FT(0))

Gets and returns primary domain information for the US-MOz (Missouri Ozark)
Fluxnet site.

The data source comes from: Unknown.
"""
function get_domain_info(::Val{:US_MOz};
    dz_bottom = FT(1.5),
    dz_top = FT(0.1),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0)
)
    
    dz_tuple = (dz_bottom, dz_top)

    return (dz_tuple=dz_tuple,
        nelements=nelements,
        zmin=zmin,
        zmax=zmax)
end

"""
    get_domain_info(::Val{:US_NR1}; dz_bottom = FT(1.25),
        dz_top = FT(0.05), nelements = 20, zmin = FT(-10), zmax = FT(0))

Gets and returns primary domain information for the US-NR1 (Colorado Niwot Ridge)
Fluxnet site.

The data source comes from: Unknown.
"""
function get_domain_info(::Val{:US_NR1};
    dz_bottom = FT(1.25),
    dz_top = FT(0.05),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0)
)

    dz_tuple = (dz_bottom, dz_top)

    return (dz_tuple=dz_tuple,
        nelements=nelements,
        zmin=zmin,
        zmax=zmax)
end

"""
    function get_domain_info(::Val{:US_Var}; dz_tuple = nothing,
        nelements = 14, zmin = FT(-0.5), zmax = FT(0)
Gets and returns primary domain information for the US-Var (California Vaira
Ranch Ione) Fluxnet site.

The data source comes from: Xu and Baldocchi, 2003.
"""
function get_domain_info(::Val{:US_Var};
    dz_tuple = nothing,
    nelements = 14,
    zmin = FT(-0.5),
    zmax = FT(0)
), where Val{:US_Var}

    return (dz_tuple=dz_tuple,
        nelements=nelements,
        zmin=zmin,
        zmax=zmax)
end

"""
    get_domain_info(site_ID::Symbol; dz_bottom = FT(1.5), dz_top = FT(0.1),
        nelements = 20, zmin = FT(-10),  zmax = FT(0))

Gets and returns primary domain information for a generic Fluxnet site,
using autofilled values from US-MOz (Missouri Ozark) site.
"""
function get_domain_info(site_ID::Symbol;
    dz_bottom = FT(1.5),
    dz_top = FT(0.1),
    nelements = 20,
    zmin = FT(-10),
    zmax = FT(0)
)
    
    dz_tuple = (dz_bottom, dz_top)

    return (dz_tuple=dz_tuple,
        nelements=nelements,
        zmin=zmin,
        zmax=zmax)
end

###################################
#            UTILITIES            #
###################################

"""
    replace_hyphen(old_site_ID::String)

Replaces all instances of hyphens in a given site ID string with underscores
and returns a Symbol of the reformatted site ID to be used as a Val{} type.
"""
function replace_hyphen(old_site_ID::String)
    new_site_ID = replace(old_site_ID, "-" => "_")

    return Symbol(new_site_ID)
end

# Things to consider:
# - fluxnet_domain.jl + fluxnet_simulation.jl - IGNORE
# - any imports? - LATER ISSUE
# - how to test - SEE TO DO
# - syntax issue for site ID - SEE TO DO

# TO DO
# make hypen replacing function in this file - DONE

# test parameters in a new file in test/fluxnet_sim.jl
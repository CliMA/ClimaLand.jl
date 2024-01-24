"""
This module is used to define plant functional types and how they can be used to
build up the parameter sets for the model. A user can specify a percentage of
land cover for each PFT within the domain over which they're running the model,
and the parameters fed to the model are then the linear combination of the PFT 
parameters weighted by the land cover percentages.
"""

export params_from_pfts, Pft

using LinearAlgebra
using DataFrames


# List of parameters that PFTs need to define. Will become the column names in
# the PFT DataFrame.
param_list = [

    # TwoStreamModel parameters
    "Ω",          # Clumping Index
    "α_PAR_leaf", # PAR Leaf Reflectance
    "α_NIR_leaf", # Near-infrared Leaf Reflectance
    "τ_PAR_leaf", # PAR Leaf Transmittance
    "τ_NIR_leaf", # Near-infrared Leaf Transmittance

    # Energy Balance model
    "ac_canopy", # Specific Heat per Emmiting Area

    # Conductance Model
    "g1", # Medlyn CO2 Sensitivity

    # Photosynthesis model
    "Vcmax25", # Max rate of carboxylation of Rubisco

    # Plant Hydraulics and general plant parameters
    "SAI",             # Avg Stem Area Index
    "f_root_to_shoot", # Root to Shoot Ratio
    "K_sat_plant",     # Hydraulic Conductivity at Saturation
    "ψ63",             # Xylem water potential at 63% loss of conductivity
    "Weibull_param",   # Conductivity Weibull Parameter
    "a",               # Retention Curve Curvature Parameter
    "capacity",        # Water Storage Capacity
    "plant_S_s",       # Specific Storativity
    "rooting_depth"    # Plant Rooting Depth
]


"""
Define a PFT type that can be used to store the parameters for each PFT.
Each PFT has a name and a list of parameters. The parameters are stored as 
NamedTuples with the parameter name as the first element and the parameter value
as the second element.
"""
struct Pft
    name::String
    parameters::NamedTuple{String, FT}
end


"""
A constructor for the PFT type which ensures that all required parameters listed
in the param_list are defined for the PFT.
"""
function Pft(name::String, parameters::List{NamedTuple{String, FT}})
    # Check that all required parameters are defined
    defined_param_names = [p[1] for p in parameters]
    for param in param_list
        if param ∉ defined_param_names
            error("PFT $name is missing parameter $param")
        end
    end
    # Return the PFT
    return Pft(name, parameters)
end


# The default PFTs used in the model use parameter values drawn from defferent 
# sources in the literature. The sources are listed here with abbreviations that
# are used in the comments in the code below.
#
#    Chen, J. M., G. Mo, J. Pisek, J. Liu, F. Deng, M. Ishizawa, 
#    and D. Chan (2012), Effects of foliage clumping on the estimation of global
#    terrestrial gross primary productivity, Global Biogeochem. Cycles, 26, 
#    GB1019, doi:10.1029/2010GB003996.
#    abb. as (Chen et al. 2012)
#
#    CLM5.0 Technical Note
#    abb. as (CLM5.0)
#
#    Best, M. J., Pryor, M., Clark, D. B., Rooney, G. G., Essery, R. L. H., 
#    Ménard, C. B., Edwards, J. M., Hendry, M. A., Porson, A., Gedney, N.,
#    Mercado, L. M., Sitch, S., Blyth, E., Boucher, O., Cox, P. M., Grimmond, 
#    C. S. B., and Harding, R. J.: The Joint UK Land Environment Simulator 
#    (JULES), model description – Part 1: Energy and water fluxes, Geosci. Model
#    Dev., 4, 677–699, https://doi.org/10.5194/gmd-4-677-2011, 2011.
#    abb as (JULES Model Description) 
#
#    Oliver, Rebecca & Mercado, Lina & Clark, Doug & Huntingford,
#    Chris & Taylor, Christopher & Vidale, P.L. & McGuire, Patrick & Todt, 
#    Markus & Folwell, Sonja & Semeena, Valiyaveetil & Medlyn, Belinda. (2022).
#    Improved representation of plant physiology in the JULES-vn5.6 land surface
#    model: Photosynthesis, stomatal conductance and thermal acclimation. 
#    10.5194/gmd-2022-11. 
#    abb. as (Oliver et al. 2022)

# Define the default PFTs that are used in the model.

NET_Temp = Pft(
    "Temperate Needleleaf Evergreen Tree",
    (
        # Radiative Transfer Model
        Ω = 0.74,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.07, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.35, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.10, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 2.35, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.080e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 3.90,      # Bonan Table 2.3
    )
)


NET_Bor = Pft(
    "Boreal Needleleaf Evergreen Tree",
    (
        # Radiative Transfer Model
        Ω = 0.74,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.07, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.35, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.10, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 2.35, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.080e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 2.0,       # Bonan Table 2.3
    )
)


NDT_Bor = Pft(
    "Boreal Needleleaf Deciduous Tree",
    (
        # Radiative Transfer Model
        Ω = 0.78,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.07, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.35, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.10, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 2.35, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.080e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 2.0,       # Bonan Table 2.3
    )
)

BET_Trop = Pft(
    "Tropical Broadleaf Evergreen Tree",
    (
        # Radiative Transfer Model
        Ω = 0.66,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.10, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.45, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.25, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.12, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 3.950e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 7.30,      # Bonan Table 2.3
    )
)


BET_Temp = Pft(
    "Temperate Broadleaf Evergreen Tree",
    (
        # Radiative Transfer Model
        Ω = 0.66,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.10, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.45, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.25, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.12, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 6.895e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 2.90,      # Bonan Table 2.3
    )
)


BDT_Trop = Pft(
    "Tropical Broadleaf Deciduous Tree",
    (
        # Radiative Transfer Model
        Ω = 0.70,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.10, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.45, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.25, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.45, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.524e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 3.70,      # Bonan Table 2.3
    )
)

BDT_Temp = Pft(
    "Temperate Broadleaf Deciduous Tree",
    (
        # Radiative Transfer Model
        Ω = 0.70,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.10, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.45, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.25, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.45, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.524e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 2.90,      # Bonan Table 2.3
    )
)


BDT_Bor = Pft(
    "Boreal Broadleaf Deciduous Tree",
    (
        # Radiative Transfer Model
        Ω = 0.70,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.10, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.45, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.25, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.45, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.524e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 2.00,      # Bonan Table 2.3
    )
)


BES_Temp = Pft(
    "Temperate Broadleaf Evergreen Shrub",
    (
        # Radiative Transfer Model
        Ω = 0.75,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.07, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.35, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.10, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.70, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 6.896e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 5.20,      # Bonan Table 2.3
    )
)


BDS_Temp = Pft(
    "Temperate Broadleaf Deciduous Shrub",
    (
        # Radiative Transfer Model
        Ω = 0.75,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.10, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.45, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.25, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.70, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.524e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 5.20,      # Bonan Table 2.3
    )
)


BDS_Bor = Pft(
    "Boreal Broadleaf Deciduous Shrub",
    (
        # Radiative Transfer Model
        Ω = 0.75,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.10, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.45, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.25, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 4.70, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 5.524e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 5.20,      # Bonan Table 2.3
    )
)


C3G_A = Pft(
    "C3 Arctic Grass",
    (
        # Radiative Transfer Model
        Ω = 0.76,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.11, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.35, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.34, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 2.22, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 4.383e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 2.60,      # Bonan Table 2.3
    )
)

C3G_NA = Pft(
    "Non-arctic C3 Grass",
    (
        # Radiative Transfer Model
        Ω = 0.76,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.11, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.35, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.34, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 5.25, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 4.383e-5, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 2.60,      # Bonan Table 2.3
    )
)


C4G = Pft(
    "C4 Grass",
    (
        # Radiative Transfer Model
        Ω = 0.75,          # Chen et al. 2012 Table 1
        α_PAR_leaf = 0.11, # CLM5.0 Table 2.3.1
        α_NIR_leaf = 0.35, # CLM5.0 Table 2.3.1
        τ_PAR_leaf = 0.05  # CLM5.0 Table 2.3.1
        τ_NIR_leaf = 0.34, # CLM5.0 Table 2.3.1

        # Canopy Energy Model
        ac_canopy = 2.5e3, # Taken from site params

        # Stomatal Conductance Model
        g1 = 1.62, # CLM5.0 Table 2.9.1

        # Photosynthesis Model,
        Vcmax25 = 4.383, # Oliver et al. 2022 Table 2

        # Hydraulics Model and Plant Parameters
        SAI = 2.50,                # Guess
        f_root_to_shoot = 3.50,    # Guess
        K_sat_plant = 2.78e-6,     # Guess
        ψ63 = -4.0 / 0.0098,       # Guess
        Weibull_param = 4.00,      # Taken from site params
        a = 0.05 * 0.0098,         # Taken from site params
        capacity = 15.00,          # Guess
        plant_S_s = 1e-2 * 0.0098, # Taken from site params
        rooting_depth = 15.0,      # Bonan Table 2.3
    )
)

# Lists of defaults to pass to parameter constructing function
default_pfts = [
    NET_Temp,
    NET_Bor,
    NDT_Bor,
    BET_Trop,
    BET_Temp,
    BDT_Trop,
    BDT_Temp,
    BDT_Bor,
    BES_Temp,
    BDS_Temp,
    BDT_Bor,
    C3G_A,
    C3G_NA,
    C4G
]


function params_from_pfts(pft_cover::Vector{Float64}, 
                            pfts::Vector{Pft}=default_pfts)
    """
    This function takes in a vector of PFT cover percentages and creates
    variables for each parameter in the model that are the linear combination of
    the PFT parameters weighted by the land cover percentages. For example,
    after running this function the variable `Ω` will be the linear combination 
    of the clumping indices of each PFT.

    May optionally take in a vector of PFTs to use instead of the default PFTs.
    In this case the pfts_cover vector must be the same length and in the
    corresponding order to the PFTs in the pfts vector. This allows a user to 
    define their own PFT scheme and plug it into the model.
    """
    # Construct the PFT DataFrame from the specified PFTs
    pft_df = DataFrame()
    for pft in pfts
        # Push the internal NamedTuple of parameters to a row in the df. Since
        # the PFT is stored as a NamedTuple conversion to a DF row is ensured 
        # to work and ensure the ordering of the parameters.
        push!(pft_df, pft.parameters)
    end
    
    # Instantiate a variable with the correct name for the weighted average of 
    # each parameter.
    for param in param_list
        # each param value is the dot product of the PFT cover percentages and
        # the parameter values
        param_value = pft_cover ⋅ pft_df[:, param]
        # create a variable with the same name as the parameter and assign it
        eval(Meta.parse("$param = $param_value"))
    end

    # Return the PFT DataFrame for reference for the user
    return pft_df
end

"""
Adds tests for the PFTs module, ensuring that PFTs cannot be created without
populating all required parameters, and ensuring that loading parameters based 
on percentage PFT cover works as expected.
"""

using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaLand
using ClimaLand.Canopy

# Test that PFTs cannot be created without all required parameters

# Create a valid parameter set for a PFT
pft_params = (
    # Radiative Transfer Model
    Ω = 0.74,          # Unitless ⋅ Chen et al. 2012 Table 1
    α_PAR_leaf = 0.07, # Unitless ⋅ CLM5.0 Table 2.3.1
    α_NIR_leaf = 0.35, # Unitless ⋅ CLM5.0 Table 2.3.1
    τ_PAR_leaf = 0.05, # Unitless ⋅ CLM5.0 Table 2.3.1
    τ_NIR_leaf = 0.10, # Unitless ⋅ CLM5.0 Table 2.3.1

    # Stomatal Conductance Model
    g1 = 74.31, # kPa^(1/2) ⋅ CLM5.0 Table 2.9.1

    # Photosynthesis Model
    Vcmax25 = 5.1e-5, # mol CO2/m^2/s ⋅ CLM2 PFT Data

    # Plant Hydraulics and general plant parameters
    SAI = 1.0,             # m^2/m^2
    f_root_to_shoot = 1.0, # Unitless ⋅ CLM2 PFT Data
    K_sat_plant = 5e-9,    # m/s
    ψ63 = -4 / 0.0098,     # / MPa to m
    capacity = 10.0,       # kg/m^2r
    rooting_depth = 3.90,  # m ⋅ Bonan Table 2.3
)

# For each parameter, remove it from the parameter set and check that the PFT
# constructor will throw an error
for param in keys(pft_params)
    # Create a copy of the parameter set minus the parameter
    pft_params_copy = Dict{Symbol, Any}()
    for (key, value) in zip(keys(pft_params), pft_params)
        if key != param
            pft_params_copy[key] = value
        end
    end
    pft_params_copy = NamedTuple{Tuple(keys(pft_params_copy))}(pft_params_copy)
    # Check that the PFT constructor will throw an error
    @test_throws ErrorException Pft("Test PFT", pft_params_copy)
end

# Test that loading parameters based on percentage PFT cover works as expected

# Create several PFTs
pft_1 = Pft("PFT1", pft_params)

pft_param_2 = (
    # Radiative Transfer Model
    Ω = 0.75,          # Unitless ⋅ Chen et al. 2012 Table 1
    α_PAR_leaf = 0.11, # Unitless ⋅ CLM5.0 Table 2.3.1
    α_NIR_leaf = 0.35, # Unitless ⋅ CLM5.0 Table 2.3.1
    τ_PAR_leaf = 0.05, # Unitless ⋅ CLM5.0 Table 2.3.1
    τ_NIR_leaf = 0.34, # Unitless ⋅ CLM5.0 Table 2.3.1

    # Stomatal Conductance Model
    g1 = 51.22, # kPa^(1/2) ⋅ CLM5.0 Table 2.9.1

    # Photosynthesis Model
    Vcmax25 = 2.4e-5, # mol CO2/m^2/s ⋅ CLM2 PFT Data

    # Plant Hydraulics and general plant parameters
    SAI = 0.0,               # m^2/m^2
    f_root_to_shoot = 3.0,   # Unitless ⋅ CLM2 PFT Data
    K_sat_plant = 5e-9,      # m/s
    ψ63 = -2.7 / 0.0098,     # / MPa to m
    capacity = 2.0,          # kg/m^2
    rooting_depth = 2.60,    # m ⋅ Bonan Table 2.3
)

pft_2 = Pft("PFT2", pft_param_2)

pft_param_3 = (
    # Radiative Transfer Model
    Ω = 0.75,          # Unitless ⋅ Chen et al. 2012 Table 1
    α_PAR_leaf = 0.10, # Unitless ⋅ CLM5.0 Table 2.3.1
    α_NIR_leaf = 0.45, # Unitless ⋅ CLM5.0 Table 2.3.1
    τ_PAR_leaf = 0.05, # Unitless ⋅ CLM5.0 Table 2.3.1
    τ_NIR_leaf = 0.25, # Unitless ⋅ CLM5.0 Table 2.3.1

    # Stomatal Conductance Model
    g1 = 148.63, # kPa^(1/2) ⋅ CLM5.0 Table 2.9.1

    # Photosynthesis Model
    Vcmax25 = 1.7e-5, # mol CO2/m^2/s ⋅ CLM2 PFT Data

    # Plant Hydraulics and general plant parameters
    SAI = 1.0,             # m^2/m^2
    f_root_to_shoot = 1.0, # Unitless ⋅ CLM2 PFT Data
    K_sat_plant = 5e-9,    # m/s
    ψ63 = -4 / 0.0098,     # / MPa to m
    capacity = 10.0,       # kg/m^2
    rooting_depth = 5.20,  # m ⋅ Bonan Table 2.3
)

pft_3 = Pft("PFT3", pft_param_3)
pfts = [pft_1, pft_2, pft_3]

# Define cover fractions with PFT1 dominant
cover_fractions = [0.5, 0.3, 0.2]

# Load the parameters based on the cover fractions
param_set = params_from_pfts(cover_fractions, pfts)

# Check that the loaded parameters are as expected - the same as param set 1
i = 1
for param in keys(pft_params)
    @test param_set[i] == pft_params[param]
    global i += 1
end

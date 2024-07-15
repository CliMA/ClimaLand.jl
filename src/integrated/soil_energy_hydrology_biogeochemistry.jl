export LandSoilBiogeochemistry, PrognosticMet

"""
    struct LandSoilBiogeochemistry{
        FT,
        SEH <: Soil.EnergyHydrology{FT},
        SB <: Soil.Biogeochemistry.SoilCO2Model{FT},
    } <: AbstractLandModel{FT}

A concrete type of land model used for simulating systems with a
soil energy, hydrology, and biogeochemistry component.
$(DocStringExtensions.FIELDS)"""
struct LandSoilBiogeochemistry{
    FT,
    SEH <: Soil.Soil.EnergyHydrology{FT},
    SB <: Soil.Biogeochemistry.SoilCO2Model{FT},
} <: AbstractLandModel{FT}
    "The soil model"
    soil::SEH
    "The biochemistry model"
    soilco2::SB
end

"""
    LandSoilBiogeochemistry{FT}(;
        soil_args::NamedTuple = (;),
        biogeochemistry_args::NamedTuple = (;),
    ) where {FT}
A constructor for the `LandSoilBiogeochemistry` model, which takes in
the required arguments for each component, constructs those models,
and constructs the `LandSoilBiogeochemistry` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.

Additional arguments, like parameters and driving atmospheric data, can be passed
in as needed.
"""
function LandSoilBiogeochemistry{FT}(;
    land_args::NamedTuple,
    soil_args::NamedTuple = (;),
    soilco2_args::NamedTuple = (;),
) where {FT}

    (; atmos, soil_organic_carbon) = land_args
    soil = Soil.EnergyHydrology{FT}(;
        soil_args..., # soil_args must have sources, boundary_conditions, domain, parameters
    )
    prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    soil_co2_drivers = Soil.Biogeochemistry.SoilDrivers(
        prognostic_soil,
        soil_organic_carbon,
        atmos,
    )
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(;
        soilco2_args...,
        drivers = soil_co2_drivers,
    )
    args = (soil, soilco2)
    return LandSoilBiogeochemistry{FT, typeof.(args)...}(args...)
end

struct PrognosticMet{FT, F <: Union{AbstractFloat, ClimaCore.Fields.Field}} <:
       Soil.Biogeochemistry.AbstractSoilDriver
    "Soil porosity (m³ m⁻³)"
    ν::F
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)"
    θ_a100::F
    "Absolute value of the slope of the line relating log(ψ) versus log(S) (unitless)"
    b::F
    function PrognosticMet(
        soil_params::Soil.EnergyHydrologyParameters{FT},
    ) where {FT}
        ν = soil_params.ν
        θ_r = soil_params.θ_r
        hcm = soil_params.hydrology_cm
        F = typeof(ν)
        if F <: AbstractFloat
            θ_a100 =
                Soil.inverse_matric_potential(hcm, -FT(1)) * (ν - θ_r) + θ_r
            b = Soil.approximate_ψ_S_slope(hcm)
        else
            θ_a100 =
                @. Soil.inverse_matric_potential(hcm, -FT(1)) * (ν - θ_r) + θ_r
            b = @. Soil.approximate_ψ_S_slope(hcm)
        end
        return new{FT, F}(ν, θ_a100, b)
    end
end

"""
    soil_temperature(driver::PrognosticSoil, p, Y, t, z)
Returns the prognostic soil temperature.
"""
function soil_temperature(driver::PrognosticMet, p, Y, t, z)
    return p.soil.T
end

"""
    soil_moisture(driver::PrognosticSoil, p, Y, t, z)

Returns the volumetric liquid fraction, computed by the soil
model from the prognostic liquid and ice fractions.
"""
function soil_moisture(driver::PrognosticMet, p, Y, t, z)
    return p.soil.θ_l
end


function ClimaLand.get_drivers(model::LandSoilBiogeochemistry)
    bc = model.soil.boundary_conditions.top
    if typeof(bc) <: AtmosDrivenFluxBC{
        <:PrescribedAtmosphere,
        <:AbstractRadiativeDrivers,
        <:Soil.AbstractRunoffModel,
    }
        return (bc.atmos, bc.radiation, model.soilco2.drivers.soc)
    else
        return (model.soilco2.drivers.atmos, model.soilco2.drivers.soc)
    end
end

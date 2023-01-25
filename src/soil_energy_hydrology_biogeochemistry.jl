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
    soil_args::NamedTuple = (;),
    soilco2_args::NamedTuple = (;),
) where {FT}

    soil = Soil.EnergyHydrology{FT}(;
        soil_args..., # soil_args must have sources, boundary_conditions, domain, parameters
    )

    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(; soilco2_args...)
    if typeof(soilco2_args.drivers.met) != PrognosticMet
        throw(
            AssertionError(
                "To run with a prognostic soil energy and hydrology model, the met driver must be of type PrognosticSoil.",
            ),
        )
    end
    if soil_args.parameters.ν != soilco2_args.parameters.ν
        throw(
            AssertionError(
                "The provided porosity must be the same in both model parameter structures.",
            ),
        )
    end

    @assert soil_args.parameters.ν == soilco2_args.parameters.ν
    args = (soil, soilco2)
    return LandSoilBiogeochemistry{FT, typeof.(args)...}(args...)
end

struct PrognosticMet <: Soil.Biogeochemistry.AbstractSoilDriver end

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

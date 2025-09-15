export LandSoilBiogeochemistry, PrognosticMet

"""
    struct LandSoilBiogeochemistry{
        FT,
        SEH <: Soil.EnergyHydrology{FT},
        SB <: Soil.Biogeochemistry.SoilCO2Model{FT},
    } <: AbstractLandModel{FT}

A concrete type of land model used for simulating systems with a
soil energy, hydrology, and biogeochemistry component.

ClimaLand v1: SoilCO2 is still under testing, but errors in soilco2
do not propagate into the other components.

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
    function LandSoilBiogeochemistry{FT}(
        soil::SEH,
        soilco2::SB,
    ) where {
        FT,
        SEH <: Soil.Soil.EnergyHydrology{FT},
        SB <: Soil.Biogeochemistry.SoilCO2Model{FT},
    }
        @assert soil.domain == soilco2.domain

        @assert soil.parameters.earth_param_set ==
                soilco2.parameters.earth_param_set

        @assert soilco2.drivers.met isa PrognosticMet
        comparison = PrognosticMet(soil.parameters)
        # check_land_equality allocates, and should only be used in initialization
        for property in propertynames(soilco2.drivers.met)
            check_land_equality(
                getproperty(soilco2.drivers.met, property),
                getproperty(comparison, property),
            )
        end
        if soil.boundary_conditions.top isa Soil.AtmosDrivenFluxBC
            @assert soil.boundary_conditions.top.atmos == soilco2.drivers.atmos
        end

        return new{FT, typeof(soil), typeof(soilco2)}(soil, soilco2)
    end
end

"""
    LandSoilBiogeochemistry{FT}(
        forcing,
        toml_dict::CP.ParamDict,
        domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell};
        soil = Soil.EnergyHydrology{FT}(
            domain,
            forcing,
            toml_dict;
        ),
        soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
            domain,
            Soil.Biogeochemistry.SoilDrivers(
               Soil.Biogeochemistry.PrognosticMet(soil.parameters),
                PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5)),
                forcing.atmos,
            ),
        ),
    ) where {FT}

A convenience constructor for setting up the default `LandSoilBiogeochemistry`,
where all the parameterizations and parameter values are set to default values
or passed in via the `toml_dict`. The boundary conditions of all models
correspond to `forcing` with the atmosphere, as specified by `forcing`, a NamedTuple
of the form `(;atmos, radiation)`, with `atmos` an `AbstractAtmosphericDriver` and `radiation`
an `AbstractRadiativeDriver`. The domain must be a ClimaLand domain with a vertical extent.
"""
function LandSoilBiogeochemistry{FT}(
    forcing,
    toml_dict::CP.ParamDict,
    domain::Union{
        ClimaLand.Domains.Column,
        ClimaLand.Domains.SphericalShell,
        ClimaLand.Domains.HybridBox,
    };
    soil = Soil.EnergyHydrology{FT}(domain, forcing, toml_dict;),
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        domain,
        Soil.Biogeochemistry.SoilDrivers(
            Soil.Biogeochemistry.PrognosticMet(soil.parameters),
            PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5)),
            forcing.atmos,
        ),
    ),
) where {FT}
    return LandSoilBiogeochemistry{FT}(soil, soilco2)
end


"""
    PrognosticMet <: AbstractSoilDriver

A container which holds the soil parameters needed for running biogeochemistry model with the
soil model.

$(DocStringExtensions.FIELDS)
"""
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

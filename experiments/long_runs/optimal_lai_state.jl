# Writes the state of `ZhouOptimalLAIModel` at the end of a simulation, to start later
# runs from it rather than from the climatology of the `optimal_lai_inputs` artifact.

import ClimaCore
import ClimaUtilities
import NCDatasets
using ClimaLand

const OPTIMAL_LAI_STATE_UNITS = Dict(
    :A0_daily => "mol CO2 m^-2 day^-1",
    :A0_annual => "mol CO2 m^-2 yr^-1",
    :precip_annual => "mol H2O m^-2 yr^-1",
    :PET_annual => "mol H2O m^-2 yr^-1",
    :VPDgs_annual => "Pa s yr^-1",
    :growing_days => "days",
    :A0c3_annual => "mol CO2 m^-2 yr^-1",
    :A0c4_annual => "mol CO2 m^-2 yr^-1",
    :GPPc3_annual => "mol CO2 m^-2 yr^-1",
    :LAI => "m^2 m^-2",
    :precip_30d => "mol H2O m^-2 (30 days)^-1",
    :PET_30d => "mol H2O m^-2 (30 days)^-1",
    :VPD_moist_annual => "Pa s",
    :moist_days => "days",
    :degree_days => "K days",
    :warm_days => "days",
    :age => "s",
)

"""
    write_optimal_lai_state(path, simulation)

Write the prognostic variables of the `ZhouOptimalLAIModel` of `simulation` at its
current time to the NetCDF file `path`, on a 1° lon-lat grid, NaN over ocean.
"""
function write_optimal_lai_state(path, simulation)
    model = simulation.model
    Y = simulation._integrator.u
    lon = collect(-180.0:1.0:179.0)
    lat = collect(-90.0:1.0:89.0)
    remapper = ClimaCore.Remapping.Remapper(
        axes(Y.canopy.biomass.LAI),
        [ClimaCore.Geometry.LatLongPoint(φ, λ) for λ in lon, φ in lat],
    )
    remap(field) = Array(ClimaCore.Remapping.interpolate(remapper, field))
    mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model))
    land =
        isnothing(mask) ? trues(length(lon), length(lat)) : remap(mask) .> 0.5
    date = ClimaUtilities.TimeManager.date(simulation._integrator.t)
    NCDatasets.NCDataset(path, "c") do ds
        ds.attrib["title"] = "State of ZhouOptimalLAIModel at the end of a simulation"
        ds.attrib["date"] = string(date)
        NCDatasets.defVar(
            ds,
            "lon",
            lon,
            ("lon",);
            attrib = ["units" => "degrees_east"],
        )
        NCDatasets.defVar(
            ds,
            "lat",
            lat,
            ("lat",);
            attrib = ["units" => "degrees_north"],
        )
        for name in ClimaLand.prognostic_vars(model.canopy.biomass)
            data = remap(getproperty(Y.canopy.biomass, name))
            NCDatasets.defVar(
                ds,
                String(name),
                Float32.(ifelse.(land, data, NaN)),
                ("lon", "lat");
                deflatelevel = 5,
                attrib = ["units" => OPTIMAL_LAI_STATE_UNITS[name]],
            )
        end
    end
    return nothing
end

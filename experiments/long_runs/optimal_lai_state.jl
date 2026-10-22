# Writes states of `ZhouOptimalLAIModel` during a simulation, to start later runs from
# them rather than from the climatology of the `optimal_lai_inputs` artifact.

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
    :snow_store => "mol H2O m^-2",
    :water_30d => "mol H2O m^-2 (30 days)^-1",
    :PET_30d => "mol H2O m^-2 (30 days)^-1",
    :VPD_moist_annual => "Pa s",
    :moist_days => "days",
    :degree_days => "K days",
    :warm_days => "days",
    :maintenance_days => "days",
    :age => "s",
)
const OPTIMAL_LAI_STATE_LON = collect(-180.0:1.0:179.0)
const OPTIMAL_LAI_STATE_LAT = collect(-90.0:1.0:89.0)

"""
    optimal_lai_state_snapshot(simulation)

The prognostic variables of the `ZhouOptimalLAIModel` of `simulation` at its current
time, on a 1° lon-lat grid, NaN over ocean, with their date: `(; date, fields)`.
"""
function optimal_lai_state_snapshot(simulation)
    model = simulation.model
    Y = simulation._integrator.u
    space = axes(Y.canopy.biomass.LAI)
    remapper = ClimaCore.Remapping.Remapper(
        space,
        [
            ClimaCore.Geometry.LatLongPoint(φ, λ) for
            λ in OPTIMAL_LAI_STATE_LON, φ in OPTIMAL_LAI_STATE_LAT
        ],
    )
    remap(field) = Array(ClimaCore.Remapping.interpolate(remapper, field))
    mask = ClimaCore.Spaces.get_mask(space)
    # A broadcast over a masked space leaves its ocean nodes unset, so the land cells
    # come from the active nodes, copied into every node of a field.
    land = if mask isa ClimaCore.DataLayouts.NoMask
        trues(length(OPTIMAL_LAI_STATE_LON), length(OPTIMAL_LAI_STATE_LAT))
    else
        active = ClimaCore.Fields.zeros(space)
        parent(active) .= parent(mask.is_active)
        remap(active) .> 0.5
    end
    fields = NamedTuple(
        name => Float32.(
            ifelse.(land, remap(getproperty(Y.canopy.biomass, name)), NaN),
        ) for name in ClimaLand.prognostic_vars(model.canopy.biomass)
    )
    date = ClimaUtilities.TimeManager.date(simulation._integrator.t)
    return (; date, fields)
end

"""
    write_optimal_lai_state(path, snapshots)
    write_optimal_lai_state(path, simulation)

Write the states of `ZhouOptimalLAIModel` in `snapshots` (see
`optimal_lai_state_snapshot`) to the NetCDF file `path`, on a (lon, lat) grid: one
state, or several along a `time` dimension, e.g. the first day of each month of a
calendar year, which runs start from at the day of the year of their start date. With
a `simulation`, write its current state.
"""
function write_optimal_lai_state(path, snapshots::AbstractVector)
    names = keys(first(snapshots).fields)
    NCDatasets.NCDataset(path, "c") do ds
        ds.attrib["title"] = "State of ZhouOptimalLAIModel"
        NCDatasets.defVar(
            ds,
            "lon",
            OPTIMAL_LAI_STATE_LON,
            ("lon",);
            attrib = ["units" => "degrees_east"],
        )
        NCDatasets.defVar(
            ds,
            "lat",
            OPTIMAL_LAI_STATE_LAT,
            ("lat",);
            attrib = ["units" => "degrees_north"],
        )
        if length(snapshots) == 1
            ds.attrib["date"] = string(only(snapshots).date)
            dims = ("lon", "lat")
            data = name -> only(snapshots).fields[name]
        else
            NCDatasets.defVar(
                ds,
                "time",
                [snapshot.date for snapshot in snapshots],
                ("time",);
                attrib = ["units" => "days since 2000-01-01 00:00:00"],
            )
            dims = ("lon", "lat", "time")
            data =
                name -> cat(
                    (snapshot.fields[name] for snapshot in snapshots)...;
                    dims = 3,
                )
        end
        for name in names
            NCDatasets.defVar(
                ds,
                String(name),
                data(name),
                dims;
                deflatelevel = 5,
                attrib = ["units" => OPTIMAL_LAI_STATE_UNITS[name]],
            )
        end
    end
    return nothing
end

write_optimal_lai_state(path, simulation) =
    write_optimal_lai_state(path, [optimal_lai_state_snapshot(simulation)])

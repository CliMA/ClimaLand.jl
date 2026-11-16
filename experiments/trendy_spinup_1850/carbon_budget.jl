import ClimaCore
import ClimaCore.Operators: column_integral_definite!
import ClimaLand
import ClimaLand.Canopy

"""
    CarbonBudget(simulation)

Global carbon fluxes and stocks of a `LandModel` with carbon pools, over land north of
60°S.

The fluxes are summed over the steps between two calls of `budget_row!`, in kg C m⁻²:
GPP, autotrophic respiration Ra, the soil surface CO2 efflux Rh, the microbial
respiration Sm and the litter input, both column-integrated. The stocks (kg C m⁻²)
are the four vegetation pools, SOC over the full column and the top metre, and the
soil CO2. `budget_row!` integrates all of them over land (Pg C) and checks the closure
ΔC_total = ∫(GPP − Ra − Rh).
"""
mutable struct CarbonBudget{F}
    sums::NamedTuple{(:GPP, :Ra, :Rh, :Sm, :litter), NTuple{5, F}}
    scratch::F
    nodes::Vector{Int}
    area::Vector{Float64}
    land_area::Float64
    seconds::Float64
    last_stocks::Union{Nothing, NamedTuple}
end

function CarbonBudget(simulation)
    Y = simulation._integrator.u
    space = axes(Y.canopy.biomass.C_leaf)
    zero_field() = ClimaCore.Fields.zeros(space)
    sums = (;
        GPP = zero_field(),
        Ra = zero_field(),
        Rh = zero_field(),
        Sm = zero_field(),
        litter = zero_field(),
    )
    nodes, area = land_nodes_and_areas(space)
    budget = CarbonBudget(sums, zero_field(), nodes, area, sum(area), 0.0, nothing)
    budget.last_stocks = global_stocks(budget, simulation)
    return budget
end

"""
    land_nodes_and_areas(space; radius = 6.371e6)

The linear indices, in the parent array of a field on the surface `space` (a regular
longitude-latitude grid), of the active (land) nodes north of 60°S, and their areas
(m²). Antarctica is left out, as in TRENDY; CRU-JRA has no data there. The nodes
masked out of the space hold no data, and the grid's Cartesian geometry does not give
areas, so integrals are taken over these nodes with these areas.
"""
function land_nodes_and_areas(space; radius = 6.371e6)
    coords = ClimaCore.Fields.coordinate_field(space)
    lat = vec(Array(parent(coords.lat)))
    lon = vec(Array(parent(coords.long)))
    mask = ClimaCore.Spaces.get_mask(space)
    active =
        mask isa ClimaCore.DataLayouts.NoMask ? trues(length(lat)) :
        vec(Array(parent(mask.is_active)))
    nodes = findall(active .& (lat .>= -60))
    dlat = π / length(unique(round.(lat; digits = 6)))
    dlon = 2π / length(unique(round.(lon; digits = 6)))
    area = @. radius^2 * cosd(lat[nodes]) * dlat * dlon
    return nodes, area
end

"""
    land_values(budget, field)

The values of the surface `field` at the land nodes of `budget`, on the host.
"""
land_values(budget::CarbonBudget, field) = Float64.(vec(Array(parent(field)))[budget.nodes])

"""
    land_integral(budget, field)

The integral over the land nodes of `budget` of the surface `field` (the field's units
times m²), computed on the host. Non-finite values are left out (counted by
`global_stocks`).
"""
land_integral(budget::CarbonBudget, field) =
    sum(ifelse(isfinite(x), x * a, 0.0) for (x, a) in zip(land_values(budget, field), budget.area))

"""
    accumulate_budget!(budget, simulation, Δt)

Adds the fluxes of the current step, times `Δt` (s), to the sums of `budget`.
"""
function accumulate_budget!(budget::CarbonBudget, simulation, Δt)
    land = simulation.model
    p = simulation._integrator.p
    M_C = land.canopy.biomass.parameters.M_C
    GPP = Canopy.get_GPP(p, land.canopy.photosynthesis)
    s = budget.sums
    @. s.GPP += GPP * M_C * Δt
    @. s.Ra += p.canopy.autotrophic_respiration.Ra * M_C * Δt
    @. s.Rh += p.soilco2.top_bc * Δt
    column_integral_definite!(budget.scratch, p.soilco2.Sm)
    @. s.Sm += budget.scratch * Δt
    column_integral_definite!(budget.scratch, p.soil_litter_input)
    @. s.litter += budget.scratch * Δt
    budget.seconds += Δt
    return nothing
end

"""
    global_stocks(budget, simulation)

Land-integrated carbon stocks (Pg C): the four vegetation pools, SOC over the full
column and the top metre, the soil CO2, and their total, over the nodes where all are
finite; and the number of land nodes where one is not.
"""
function global_stocks(budget::CarbonBudget, simulation)
    land = simulation.model
    Y = simulation._integrator.u
    b = Y.canopy.biomass
    to_Pg = 1e-12
    bad = falses(length(budget.nodes))
    function pool(x)
        bad .|= .!isfinite.(land_values(budget, x))
        return land_integral(budget, x) * to_Pg
    end
    column(x) = (column_integral_definite!(budget.scratch, x); pool(budget.scratch))
    z = land.soil.domain.fields.z
    top_metre_SOC = @. ifelse(z >= -1, Y.soilco2.SOC, zero(Y.soilco2.SOC))
    stocks = (;
        C_sugar = pool(b.C_sugar),
        C_leaf = pool(b.C_leaf),
        C_stem = pool(b.C_stem),
        C_root = pool(b.C_root),
        SOC = column(Y.soilco2.SOC),
        SOC_1m = column(top_metre_SOC),
        soil_CO2 = column(Y.soilco2.CO2),
    )
    total =
        stocks.C_sugar +
        stocks.C_leaf +
        stocks.C_stem +
        stocks.C_root +
        stocks.SOC +
        stocks.soil_CO2
    return merge(stocks, (; C_total = total, nonfinite_nodes = count(bad)))
end

"""
    budget_row!(budget, simulation)

The land-integrated fluxes (Pg C over the period since the last call) and the stocks
(Pg C) as a `NamedTuple`, with the closure residual ΔC_total − ∫(GPP − Ra − Rh). Resets
the sums.
"""
function budget_row!(budget::CarbonBudget, simulation)
    to_Pg = 1e-12
    flux(x) = land_integral(budget, x) * to_Pg
    s = budget.sums
    fluxes = (;
        GPP = flux(s.GPP),
        Ra = flux(s.Ra),
        Rh = flux(s.Rh),
        Sm = flux(s.Sm),
        litter = flux(s.litter),
    )
    NBP = fluxes.GPP - fluxes.Ra - fluxes.Rh
    stocks = global_stocks(budget, simulation)
    ΔC = stocks.C_total - budget.last_stocks.C_total
    LAI = land_integral(budget, simulation._integrator.p.canopy.biomass.area_index.leaf) / budget.land_area
    row = merge(
        (; days = budget.seconds / 86400),
        fluxes,
        (; NBP, closure = ΔC - NBP, LAI_end = LAI),
        stocks,
        (; land_area_Mkm2 = budget.land_area * 1e-12),
    )
    for x in s
        x .= 0
    end
    budget.seconds = 0.0
    budget.last_stocks = stocks
    return row
end

"""
    append_row(path, row; prefix = (;))

Appends `merge(prefix, row)` to the CSV at `path`, writing the header if the file is new.
"""
function append_row(path, row; prefix = (;))
    full = merge(prefix, row)
    isnew = !isfile(path)
    open(path, "a") do io
        isnew && println(io, join(string.(keys(full)), ","))
        println(io, join(string.(values(full)), ","))
    end
end

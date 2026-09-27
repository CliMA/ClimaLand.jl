#=
Surface layer between the mineral soil and the atmosphere.

`NoLitter` keeps the soil skin temperature scheme of `soil_surface_temperature.jl`:
a zero-heat-capacity skin connected to the top soil cell by the half-cell
conduction resistance plus a constant litter resistance.

`SlabLitter` adds a litter slab with heat capacity C_l = ρc_l d_l and
temperature T_l (`Y.soil.T_litter`) between the skin and the soil:

    atmosphere      F_atm(T_sfc) = SW_n + LW_n + H + L        (positive upward)
    ── skin T_sfc   zero heat capacity: F_atm = (T_l - T_sfc) / r_top
       r_top = d_l / (2 κ_l)
    ── litter T_l   C_l dT_l/dt = -F_atm + (T_top - T_l) / r_bot
       r_bot = d_l / (2 κ_l) + Δz_top / κ_top
    ── soil T_top   top energy flux (upward) = (T_top - T_l) / r_bot + F_soil

`F_soil` collects fluxes that enter the soil directly: infiltration energy,
snow excess fluxes, and the conduction from an overlying snowpack, which
exchanges with the top soil cell so that the litter under snow follows the
soil.
Water vapor from the soil diffuses through the slab, whose resistance
c_vap (d_l - d_min) / D_vapor adds in series to the dry-soil-layer resistance
in `soil_surface_vapor_conductance!`.

The litter thickness follows the plant area index of the canopy above with a
memory of about a year, d_l = max(d_PAI ⟨PAI⟩, d_min), where ⟨PAI⟩
(`Y.soil.PAI_mean`) is an exponentially weighted trailing mean of LAI + SAI:
litter stock is litterfall times residence time, and both average over
phenology. Without a canopy, ⟨PAI⟩ decays to zero and the slab reduces to the
skin scheme.

T_l relaxes on C_l / ∂F_atm∂T ~ 10² s, shorter than the time step, so it is
treated implicitly. The skin balance is solved once per step at T_l^n and the
atmospheric flux is linearized about it, F_atm(T_l) ≈ F_atm^n + Λ (T_l - T_l^n)
with Λ = ∂F_atm/∂T_l from the skin solve. The backward Euler litter equation is
then linear in T_l given T_top and is solved for T_l in closed form in every
Newton iteration of the soil solve; the soil top flux becomes a function of
T_top alone, and its derivative enters the soil energy Jacobian at the top
face. Y.soil.T_litter is advanced to this solution, and the linearized
atmospheric flux is used for the energy bookkeeping `∫F_e_dt`, together with
the energy of the litter mass gained or lost as ⟨PAI⟩ changes.
=#

"""
    AbstractSoilSurfaceLayer{FT}

Layer between the mineral soil and the atmosphere of an `EnergyHydrology` model.

Subtypes:
- [`NoLitter`](@ref): the soil skin exchanges with the top soil cell through a
  thermal resistance.
- [`SlabLitter`](@ref): a litter slab with heat capacity and a prognostic
  temperature between the skin and the soil.
"""
abstract type AbstractSoilSurfaceLayer{FT <: AbstractFloat} end

"""
    NoLitter{FT} <: AbstractSoilSurfaceLayer{FT}

Surface layer without heat capacity: the soil skin exchanges with the top soil
cell through the half-cell conduction resistance plus the litter resistance
`r_litter` passed to the skin solve (the canopy's `biomass.r_litter` in land
models with a canopy, zero for bare soil).
"""
struct NoLitter{FT} <: AbstractSoilSurfaceLayer{FT} end

"""
    SlabLitter{FT} <: AbstractSoilSurfaceLayer{FT}

Litter slab between the soil skin and the top soil cell, with temperature
`Y.soil.T_litter` and thickness `d_l = max(d_PAI * PAI_mean, d_min)`, where
`Y.soil.PAI_mean` is a trailing mean of the plant area index of the canopy
with memory `τ_PAI`.

# Fields
- `d_PAI`: Litter thickness per unit plant area index [m].
- `d_min`: Minimum litter thickness [m].
- `τ_PAI`: Memory timescale of the trailing mean plant area index [s].
- `κ_l`: Litter thermal conductivity [W/m/K].
- `ρc_l`: Litter volumetric heat capacity [J/m³/K].
- `c_vap`: Factor scaling the molecular vapor diffusion resistance
  `(d_l - d_min) / D_vapor` of the slab [-]; values below one represent
  ventilation of the loose litter by air motion.
- `Δt`: Simulation time step [s]; the litter temperature is advanced by a
  backward Euler step of this length within the soil solve.

# Constructor
    SlabLitter{FT}(toml_dict, Δt; d_PAI, d_min, τ_PAI, κ_l, ρc_l, c_vap)

Keyword arguments default to the TOML parameters `litter_thickness_per_pai`,
`litter_minimum_thickness`, `litter_memory_timescale`,
`litter_thermal_conductivity`, `litter_volumetric_heat_capacity`, and
`litter_vapor_resistance_factor`.

# Examples
```julia
soil = EnergyHydrology{FT}(
    domain,
    forcing,
    toml_dict;
    surface_layer = SlabLitter{FT}(toml_dict, 450),
)
```
"""
struct SlabLitter{FT} <: AbstractSoilSurfaceLayer{FT}
    d_PAI::FT
    d_min::FT
    τ_PAI::FT
    κ_l::FT
    ρc_l::FT
    c_vap::FT
    Δt::FT
end

function SlabLitter{FT}(
    toml_dict::CP.ParamDict,
    Δt;
    d_PAI = toml_dict["litter_thickness_per_pai"],
    d_min = toml_dict["litter_minimum_thickness"],
    τ_PAI = toml_dict["litter_memory_timescale"],
    κ_l = toml_dict["litter_thermal_conductivity"],
    ρc_l = toml_dict["litter_volumetric_heat_capacity"],
    c_vap = toml_dict["litter_vapor_resistance_factor"],
) where {FT}
    return SlabLitter{FT}(
        FT(d_PAI),
        FT(d_min),
        FT(τ_PAI),
        FT(κ_l),
        FT(ρc_l),
        FT(c_vap),
        FT(float(Δt)),
    )
end

# Prognostic and auxiliary variables contributed to the soil model
surface_layer_prognostic_vars(::NoLitter) = ()
surface_layer_prognostic_vars(::SlabLitter) = (:T_litter, :PAI_mean)
surface_layer_prognostic_types(::NoLitter) = ()
surface_layer_prognostic_types(::SlabLitter{FT}) where {FT} = (FT, FT)
surface_layer_prognostic_domain_names(::NoLitter) = ()
surface_layer_prognostic_domain_names(::SlabLitter) = (:surface, :surface)

surface_layer_aux_vars(::NoLitter) = ()
surface_layer_aux_vars(::SlabLitter) =
    (:litter, :skin_solve, :dfluxBCdY_heat, :topBC_heat_scratch)
surface_layer_aux_types(::NoLitter, soil) = ()
surface_layer_aux_types(::SlabLitter{FT}, soil) where {FT} = (
    NamedTuple{
        (:F_atm, :∂F_atm∂T, :F_soil, :T_n, :T),
        Tuple{FT, FT, FT, FT, FT},
    },
    skin_solve_type(
        first(
            boundary_var_types(
                soil,
                soil.boundary_conditions.top,
                ClimaLand.TopBoundary(),
            ),
        ),
        FT,
    ),
    Geometry.Covariant3Vector{FT},
    Geometry.Covariant3Vector{FT},
)

# The skin solve of a `SlabLitter` stores the turbulent fluxes, the skin
# temperature, and the flux sensitivity `∂F∂T` together
skin_solve_type(::Type{NamedTuple{names, T}}, FT) where {names, T} =
    NamedTuple{(names..., :∂F∂T), Tuple{fieldtypes(T)..., FT}}
surface_layer_aux_domain_names(::NoLitter) = ()
surface_layer_aux_domain_names(::SlabLitter) =
    (:surface, :surface, :surface, :subsurface_face)

"""
    litter_thickness(sl::SlabLitter, Y)

Return the litter thickness `max(d_PAI * PAI_mean, d_min)` [m] as a lazy
broadcast over `Y.soil.PAI_mean`.
"""
litter_thickness(sl::SlabLitter, Y) =
    @. lazy(max(sl.d_PAI * Y.soil.PAI_mean, sl.d_min))

"""
    litter_half_resistance(sl::SlabLitter, Y)

Return the thermal resistance between the litter node and either face of the
slab, `d_l / (2 κ_l)` [m² K/W].
"""
function litter_half_resistance(sl::SlabLitter, Y)
    d_l = litter_thickness(sl, Y)
    return @. lazy(d_l / (2 * sl.κ_l))
end

"""
    litter_heat_capacity(sl::SlabLitter, Y)

Return the heat capacity per unit area of the litter slab, `ρc_l d_l` [J/m²/K].
"""
function litter_heat_capacity(sl::SlabLitter, Y)
    d_l = litter_thickness(sl, Y)
    return @. lazy(sl.ρc_l * d_l)
end

"""
    litter_soil_resistance(sl::SlabLitter, model, Y, p)

Return the thermal resistance between the litter node and the center of the
top soil cell, `d_l / (2 κ_l) + Δz_top / κ_top` [m² K/W].
"""
function litter_soil_resistance(sl::SlabLitter, model, Y, p)
    κ_top = ClimaLand.Domains.top_center_to_surface(p.soil.κ)
    Δz_top = model.domain.fields.Δz_top
    r_half = litter_half_resistance(sl, Y)
    ε = eps(eltype(Δz_top))
    return @. lazy(r_half + Δz_top / max(κ_top, ε))
end

"""
    litter_vapor_resistance(sl::AbstractSoilSurfaceLayer, Y, _D_vapor)

Return the diffusive resistance [s/m] of the litter layer to water vapor:
`c_vap (d_l - d_min) / _D_vapor` for [`SlabLitter`](@ref) and zero for
[`NoLitter`](@ref). At the bare-soil floor `d_l = d_min` the resistance
vanishes.
"""
litter_vapor_resistance(::NoLitter, Y, _D_vapor::FT) where {FT} = FT(0)
function litter_vapor_resistance(sl::SlabLitter, Y, _D_vapor)
    d_l = litter_thickness(sl, Y)
    return @. lazy(sl.c_vap * (d_l - sl.d_min) / _D_vapor)
end

"""
    add_litter_vapor_resistance!(g_soil_sfc, sl::AbstractSoilSurfaceLayer, Y, _D_vapor)

Combine the mineral soil vapor conductance `g_soil_sfc` [m/s] in series with
the vapor diffusion resistance of the surface layer `sl` in place; return
`nothing`. For [`NoLitter`](@ref) this is a no-op.
"""
add_litter_vapor_resistance!(g_soil_sfc, ::NoLitter, Y, _D_vapor) = nothing
function add_litter_vapor_resistance!(g_soil_sfc, sl::SlabLitter, Y, _D_vapor)
    r_l = litter_vapor_resistance(sl, Y, _D_vapor)
    @. g_soil_sfc = 1 / (1 / g_soil_sfc + r_l)
    return nothing
end

"""
    plant_area_index_above_soil(p, prognostic_land_components, FT)

Return the plant area index `LAI + SAI` [-] of the canopy above the soil from
`p.canopy.biomass.area_index`, or zero when the land model has no canopy
component.
"""
function plant_area_index_above_soil(
    p,
    ::Val{components},
    FT,
) where {components}
    if :canopy in components
        area_index = p.canopy.biomass.area_index
        return @. lazy(area_index.leaf + area_index.stem)
    else
        return FT(0)
    end
end

"""
    skin_lower_node(sl, model, Y, p, r_litter)

Return the temperature of the node below the soil skin and the thermal
resistance [m² K/W] between them: the top soil cell through
`Δz_top/κ_top + r_litter` for [`NoLitter`](@ref), the litter node through half
the slab for [`SlabLitter`](@ref).

Called from [`update_soil_surface_temperature!`](@ref).
"""
function skin_lower_node(::NoLitter, model, Y, p, r_litter)
    T_top = ClimaLand.Domains.top_center_to_surface(p.soil.T)
    κ_top = ClimaLand.Domains.top_center_to_surface(p.soil.κ)
    Δz_top = model.domain.fields.Δz_top
    r = @. lazy(soil_surface_thermal_resistance(Δz_top, κ_top, r_litter))
    return T_top, r
end
skin_lower_node(sl::SlabLitter, model, Y, p, r_litter) =
    Y.soil.T_litter, litter_half_resistance(sl, Y)

without_sensitivity(x::NamedTuple{names}) where {names} =
    NamedTuple{Base.front(names)}(Base.front(Tuple(x)))

"""
    store_skin_solution!(sl, p, skin)

Store the result of the skin temperature solve, a lazy broadcast of the
turbulent fluxes and the skin temperature (see
[`solve_soil_surface_temperature_at_a_point`](@ref)), in
`p.soil.turbulent_fluxes`. For [`SlabLitter`](@ref), the solve also returns the
sensitivity `∂F∂T` of the atmospheric flux to the litter temperature, and the
full result is stored in `p.soil.skin_solve` first, so that the solve runs
once.

Called from [`update_soil_surface_temperature!`](@ref).
"""
function store_skin_solution!(::NoLitter, p, skin)
    p.soil.turbulent_fluxes .= skin
    return nothing
end
function store_skin_solution!(::SlabLitter, p, skin)
    p.soil.skin_solve .= skin
    p.soil.turbulent_fluxes .= without_sensitivity.(p.soil.skin_solve)
    return nothing
end

"""
    skin_flux_sensitivity(sl, p, FT)

Return the derivative of the atmospheric energy flux at the skin with respect
to the temperature of the node below it [W/m²/K]: zero for [`NoLitter`](@ref),
which does not use it, and the value stored by the skin solve for
[`SlabLitter`](@ref).
"""
skin_flux_sensitivity(::NoLitter, p, FT) = FT(0)
skin_flux_sensitivity(::SlabLitter, p, FT) = p.soil.skin_solve.∂F∂T

"""
    set_soil_top_heat_flux!(sl, model, Y, p; F_atm, ∂F_atm∂T, F_soil)

Set the top energy flux boundary condition of the soil, `p.soil.top_bc.heat`
(positive upward), from the fluxes at the surface. Return `nothing`.

For [`NoLitter`](@ref) the soil receives `F_atm + F_soil`. For
[`SlabLitter`](@ref) the atmospheric flux acts on the litter: the arguments and
the start-of-step litter temperature are stored in `p.soil.litter` for the
implicit solve, and the soil receives the litter-soil conduction plus `F_soil`
(see [`update_litter_soil_heat_flux!`](@ref)).

# Keyword Arguments
- `F_atm`: Net atmospheric energy flux at the skin [W/m²].
- `∂F_atm∂T`: Derivative of `F_atm` with respect to the temperature below the
  skin [W/m²/K].
- `F_soil`: Fluxes entering the soil directly: the energy of infiltrating
  water, snow excess fluxes, and the conduction from an overlying snowpack
  [W/m²].

Fluxes over partially snow-covered ground are passed weighted by area fraction.
Called from the `soil_boundary_fluxes!` methods.
"""
function set_soil_top_heat_flux!(
    ::NoLitter,
    model,
    Y,
    p;
    F_atm,
    ∂F_atm∂T,
    F_soil,
)
    @. p.soil.top_bc.heat = F_atm + F_soil
    return nothing
end
function set_soil_top_heat_flux!(
    sl::SlabLitter,
    model,
    Y,
    p;
    F_atm,
    ∂F_atm∂T,
    F_soil,
)
    litter = p.soil.litter
    @. litter.F_atm = F_atm
    @. litter.∂F_atm∂T = ∂F_atm∂T
    @. litter.F_soil = F_soil
    litter.T_n .= Y.soil.T_litter
    update_litter_soil_heat_flux!(sl, model, Y, p)
    return nothing
end

"""
    add_soil_top_heat_flux!(sl, model, Y, p, F)

Add the energy flux `F` [W/m²] (positive upward) to the fluxes that enter the
soil directly, after [`set_soil_top_heat_flux!`](@ref) has run; return
`nothing`. For [`SlabLitter`](@ref) the flux is stored in `p.soil.litter.F_soil`
so that the implicit recomputation of the soil top flux keeps it.

Called from `update_soil_heat_flux_with_lake_sediment_flux!`.
"""
function add_soil_top_heat_flux!(::NoLitter, model, Y, p, F)
    @. p.soil.top_bc.heat += F
    return nothing
end
function add_soil_top_heat_flux!(sl::SlabLitter, model, Y, p, F)
    @. p.soil.litter.F_soil += F
    update_litter_soil_heat_flux!(sl, model, Y, p)
    return nothing
end

"""
    litter_temperature(T_top, T_n, F_atm, Λ, C_l, r_bot, Δt)

Return the litter temperature [K] that solves the backward Euler energy balance
`C_l (T_l - T_n)/Δt = -F_atm - Λ (T_l - T_n) + (T_top - T_l)/r_bot`
for a given top soil cell temperature `T_top`.

# Arguments
- `T_top`: Temperature of the top soil cell [K].
- `T_n`: Litter temperature at the start of the step [K].
- `F_atm`: Atmospheric energy flux at the start of the step [W/m²].
- `Λ`: Derivative of the atmospheric flux with respect to the litter
  temperature [W/m²/K].
- `C_l`: Heat capacity of the litter per unit area [J/m²/K].
- `r_bot`: Thermal resistance between litter and top soil cell [m² K/W].
- `Δt`: Time step [s].
"""
function litter_temperature(T_top, T_n, F_atm, Λ, C_l, r_bot, Δt)
    a = C_l / Δt + Λ + 1 / r_bot
    return (C_l / Δt * T_n - F_atm + Λ * T_n + T_top / r_bot) / a
end

"""
    litter_flux_derivative(Λ, C_l, r_bot, Δt)

Return the derivative [W/m²/K] with respect to `T_top` of the litter-soil
conduction `(T_top - T_l(T_top))/r_bot`, with `T_l` from
[`litter_temperature`](@ref).
"""
function litter_flux_derivative(Λ, C_l, r_bot, Δt)
    a = C_l / Δt + Λ + 1 / r_bot
    return (1 - 1 / (r_bot * a)) / r_bot
end

"""
    update_litter_soil_heat_flux!(sl, model, Y, p)

Solve the litter energy balance for the litter temperature given the current
top soil temperature `p.soil.T` and update the cache in place; return
`nothing`. For [`NoLitter`](@ref) this is a no-op.

Reads `p.soil.T`, `p.soil.κ`, `p.soil.θ_l`, `Y.soil.θ_i`, and the fluxes stored
in `p.soil.litter`. Mutates `p.soil.litter.T` (the litter temperature after the
step), `p.soil.top_bc.heat` (the litter-soil conduction plus `F_soil`), and
`p.soil.dfluxBCdY_heat` (the derivative of that flux with respect to the top
cell internal energy, as a covariant vector on the top face, for the
Jacobian).

Called in the implicit stage of every Newton iteration, so that the coupling
is implicit in both temperatures. See also
[`add_surface_layer_heat_flux_jacobian!`](@ref).
"""
update_litter_soil_heat_flux!(::NoLitter, model, Y, p) = nothing
function update_litter_soil_heat_flux!(sl::SlabLitter, model, Y, p)
    (; ρc_ds, earth_param_set) = model.parameters
    litter = p.soil.litter
    T_top = ClimaLand.Domains.top_center_to_surface(p.soil.T)
    r_bot = litter_soil_resistance(sl, model, Y, p)
    C_l = litter_heat_capacity(sl, Y)
    @. litter.T = litter_temperature(
        T_top,
        litter.T_n,
        litter.F_atm,
        litter.∂F_atm∂T,
        C_l,
        r_bot,
        sl.Δt,
    )
    @. p.soil.top_bc.heat = (T_top - litter.T) / r_bot + litter.F_soil

    # ∂F_top/∂ρe_top = ∂F_top/∂T_top / ρc_top on the top face
    θ_l_top = ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)
    θ_i_top = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    ρc_ds_top = ClimaLand.Domains.top_center_to_surface(ρc_ds)
    local_geometry_faceN = ClimaLand.Domains.top_face_to_surface(
        ClimaCore.Fields.local_geometry_field(
            ClimaCore.Spaces.face_space(axes(p.soil.κ)),
        ),
        axes(p.soil.dfluxBCdY_heat),
    )
    @. p.soil.dfluxBCdY_heat =
        covariant3_unit_vector(local_geometry_faceN) * (
            litter_flux_derivative(litter.∂F_atm∂T, C_l, r_bot, sl.Δt) /
            volumetric_heat_capacity(
                θ_l_top,
                θ_i_top,
                ρc_ds_top,
                earth_param_set,
            )
        )
    return nothing
end

"""
    linearized_atmos_flux(sl::SlabLitter, p)

Return the atmospheric energy flux at the skin [W/m²] linearized about the
start of the step, `F_atm^n + ∂F_atm∂T (T_l - T_l^n)`, at the current litter
temperature `p.soil.litter.T`.
"""
function linearized_atmos_flux(::SlabLitter, p)
    litter = p.soil.litter
    return @. lazy(litter.F_atm + litter.∂F_atm∂T * (litter.T - litter.T_n))
end

"""
    soil_column_top_energy_flux(sl, Y, p)

Return the energy flux [W/m²] (positive upward) leaving the soil and its
surface layer at the top, which `∫F_e_dt` accumulates: `p.soil.top_bc.heat` for
[`NoLitter`](@ref), and the linearized atmospheric flux plus the direct soil
fluxes for [`SlabLitter`](@ref).
"""
soil_column_top_energy_flux(::NoLitter, Y, p) = p.soil.top_bc.heat
function soil_column_top_energy_flux(sl::SlabLitter, Y, p)
    litter = p.soil.litter
    F_atm = linearized_atmos_flux(sl, p)
    return @. lazy(F_atm + litter.F_soil)
end

"""
    surface_layer_exp_tendency!(sl, dY, Y, p, model)

Set the explicit tendencies of the surface layer state in `dY`; return
`nothing`. For [`NoLitter`](@ref) this is a no-op.

For [`SlabLitter`](@ref), `dY.soil.T_litter` is zero (the litter temperature is
advanced in the implicit stage), `dY.soil.PAI_mean` relaxes the trailing mean
toward the canopy plant area index with memory `τ_PAI`, and the energy of the
litter mass gained or lost as the thickness changes, `(T_l - T_0) dC_l/dt`, is
added to `dY.soil.∫F_e_dt`.

Called from `make_compute_exp_tendency(::EnergyHydrology)`.
"""
surface_layer_exp_tendency!(::NoLitter, dY, Y, p, model) = nothing
function surface_layer_exp_tendency!(sl::SlabLitter, dY, Y, p, model)
    FT = eltype(Y.soil.T_litter)
    dY.soil.T_litter .= 0
    PAI = plant_area_index_above_soil(
        p,
        Val(model.boundary_conditions.top.prognostic_land_components),
        FT,
    )
    reduction = ClimaLand.RunningMean(sl.τ_PAI)
    @. dY.soil.PAI_mean =
        ClimaLand.apply_time_reduction(PAI, Y.soil.PAI_mean, reduction)
    _T_ref = FT(LP.T_0(model.parameters.earth_param_set))
    # dC_l/dt = ρc_l d_PAI dPAI_mean/dt while the thickness is above its floor
    @. dY.soil.∫F_e_dt +=
        (Y.soil.T_litter - _T_ref) *
        sl.ρc_l *
        sl.d_PAI *
        dY.soil.PAI_mean *
        (sl.d_PAI * Y.soil.PAI_mean > sl.d_min)
    return nothing
end

"""
    surface_layer_imp_tendency!(sl, dY, Y, p, model)

Set the implicit tendencies of the surface layer state in `dY`; return
`nothing`. For [`SlabLitter`](@ref), `dY.soil.T_litter` advances the litter
temperature to the solution `p.soil.litter.T` of the implicit stage over the
time step. For [`NoLitter`](@ref) this is a no-op.

Called from `make_compute_imp_tendency(::EnergyHydrology)`.
"""
surface_layer_imp_tendency!(::NoLitter, dY, Y, p, model) = nothing
function surface_layer_imp_tendency!(sl::SlabLitter, dY, Y, p, model)
    litter = p.soil.litter
    @. dY.soil.T_litter = (litter.T - litter.T_n) / sl.Δt
    dY.soil.PAI_mean .= 0
    return nothing
end

"""
    add_surface_layer_heat_flux_jacobian!(sl, model, p, interpc2f_op, FT)

Add the derivative of the soil top energy flux with respect to the top cell
internal energy, `p.soil.dfluxBCdY_heat`, to the top face of the heat flux
matrix `p.soil.full_bidiag_matrix_scratch`; return `nothing`. For
[`NoLitter`](@ref) this is a no-op.

Called from `make_compute_jacobian(::EnergyHydrology)` before the flux matrix
is differenced into the Jacobian block of `ρe_int`.
"""
add_surface_layer_heat_flux_jacobian!(::NoLitter, model, p, interpc2f_op, FT) =
    nothing
function add_surface_layer_heat_flux_jacobian!(
    ::SlabLitter,
    model,
    p,
    interpc2f_op,
    FT,
)
    topBC_op = Operators.SetBoundaryOperator(
        top = Operators.SetValue(p.soil.dfluxBCdY_heat),
        bottom = Operators.SetValue(Geometry.Covariant3Vector(zero(FT))),
    )
    @. p.soil.topBC_heat_scratch =
        topBC_op(Geometry.Covariant3Vector(zero(interpc2f_op(p.soil.κ))))
    @. p.soil.full_bidiag_matrix_scratch +=
        MatrixFields.LowerDiagonalMatrixRow(p.soil.topBC_heat_scratch)
    return nothing
end

"""
    add_surface_layer_energy!(surface_field, sl, Y, model)

Add the energy per unit area of the surface layer [J/m²], `C_l (T_l - T_0)`
with the reference temperature `T_0` of the soil internal energy, to
`surface_field`; return `nothing`. For [`NoLitter`](@ref) this is a no-op.

Called from `total_energy_per_area!(::EnergyHydrology)`.
"""
add_surface_layer_energy!(surface_field, ::NoLitter, Y, model) = nothing
function add_surface_layer_energy!(surface_field, sl::SlabLitter, Y, model)
    FT = eltype(surface_field)
    _T_ref = FT(LP.T_0(model.parameters.earth_param_set))
    C_l = litter_heat_capacity(sl, Y)
    @. surface_field += C_l * (Y.soil.T_litter - _T_ref)
    return nothing
end

"""
    initialize_litter_temperature!(Y, model)

Set the litter temperature `Y.soil.T_litter` to the temperature of the top soil
cell computed from the soil state; return `nothing`. A no-op for models
without a litter layer.
"""
initialize_litter_temperature!(Y, model) =
    initialize_litter_temperature!(model.surface_layer, Y, model)
initialize_litter_temperature!(::NoLitter, Y, model) = nothing
function initialize_litter_temperature!(::SlabLitter, Y, model)
    (; ν, θ_r, ρc_ds, earth_param_set) = model.parameters
    θ_l = @. lazy(volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i, θ_r))
    ρc_s = @. lazy(
        volumetric_heat_capacity(θ_l, Y.soil.θ_i, ρc_ds, earth_param_set),
    )
    # Allocates a column field; only called when setting initial conditions
    T = @. temperature_from_ρe_int(
        Y.soil.ρe_int,
        Y.soil.θ_i,
        ρc_s,
        earth_param_set,
    )
    Y.soil.T_litter .= ClimaLand.Domains.top_center_to_surface(T)
    return nothing
end

"""
    initialize_litter_area_index!(Y, model, PAI)

Set the trailing mean plant area index `Y.soil.PAI_mean` to `PAI` [-], a scalar
or a field on the surface space; return `nothing`. A no-op for models without
a litter layer.
"""
initialize_litter_area_index!(Y, model, PAI) =
    initialize_litter_area_index!(model.surface_layer, Y, PAI)
initialize_litter_area_index!(::NoLitter, Y, PAI) = nothing
initialize_litter_area_index!(::SlabLitter, Y, PAI) =
    (Y.soil.PAI_mean .= PAI; nothing)

"""
    check_time_step(sl::AbstractSoilSurfaceLayer, Δt)

Throw an `ArgumentError` if the time step stored in a [`SlabLitter`](@ref)
differs from the simulation time step `Δt` [s]; return `nothing` otherwise.
The litter temperature is advanced by a backward Euler step of the stored
length, so a mismatch would break the energy budget.
"""
check_time_step(::AbstractSoilSurfaceLayer, Δt) = nothing
function check_time_step(sl::SlabLitter, Δt)
    isapprox(sl.Δt, float(Δt)) || throw(
        ArgumentError(
            "SlabLitter was created for a time step of $(sl.Δt) s but the simulation time step is $(float(Δt)) s",
        ),
    )
    return nothing
end

"""
    check_timestepper(sl::AbstractSoilSurfaceLayer, timestepper)

Throw an `ArgumentError` unless `timestepper` is an IMEX algorithm with the
`ARS111` tableau and at least two Newton iterations per step when `sl` is a
[`SlabLitter`](@ref); return `nothing` otherwise. The litter temperature is
advanced by a backward Euler substep of the full time step, which matches only
`ARS111`, and its update lags the soil by one Newton iteration, so a single
iteration would open the energy budget.
"""
check_timestepper(::AbstractSoilSurfaceLayer, timestepper) = nothing
function check_timestepper(::SlabLitter, timestepper)
    (
        timestepper isa ClimaTimeSteppers.IMEXAlgorithm &&
        timestepper.name isa ClimaTimeSteppers.ARS111
    ) || throw(
        ArgumentError(
            "SlabLitter requires the ARS111 (backward Euler) IMEX time stepper",
        ),
    )
    max_iters = timestepper.newtons_method.max_iters
    max_iters >= 2 || throw(
        ArgumentError(
            "SlabLitter requires at least two Newton iterations per step; max_iters = $max_iters",
        ),
    )
    return nothing
end

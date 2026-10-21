import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
using Dates
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import Interpolations: Constant
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
using DocStringExtensions

export prescribed_lai_era5,
    prescribed_lai_modis,
    prescribed_climatological_lai_modis,
    PrescribedAreaIndices,
    update_biomass!,
    ZhouOptimalLAIModel,
    mask_biomass!

"""
     prescribed_lai_era5(era5_lai_ncdata_path,
                         era5_lai_cover_ncdata_path,
                         surface_space,
                         start_date,
                         earth_param_set;
                         time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                         regridder_type = :InterpolationsRegridder,
                         interpolation_method = Interpolations.Constant(),)

A helper function which constructs the TimeVaryingInput object for Leaf Area Index, from a
file path pointing to the ERA5 LAI data in a netcdf file, a file path pointing to the ERA5
LAI cover data in a netcdf file, the surface_space, the start date, and the earth_param_set.

This currently one works when a single file is passed for both the era5 lai and era5 lai cover data.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported
by passing `interpolation_method = Interpolations.Linear()`.
"""
function prescribed_lai_era5(
    era5_lai_ncdata_path,
    era5_lai_cover_ncdata_path,
    surface_space,
    start_date;
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
)
    hvc_ds = NCFileReader(era5_lai_cover_ncdata_path, "cvh")
    lvc_ds = NCFileReader(era5_lai_cover_ncdata_path, "cvl")
    hv_cover = read(hvc_ds)
    lv_cover = read(lvc_ds)
    close(hvc_ds)
    close(lvc_ds)
    compose_function = let hv_cover = hv_cover, lv_cover = lv_cover
        (lai_hv, lai_lv) -> lai_hv .* hv_cover .+ lai_lv .* lv_cover
    end
    return TimeVaryingInput(
        era5_lai_ncdata_path,
        ["lai_hv", "lai_lv"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
        compose_function = compose_function,
    )
end

"""
     prescribed_lai_modis(surface_space,
                          start_date
                          stop_date,
                          earth_param_set;
                          time_interpolation_method = LinearInterpolation(),
                          regridder_type = :InterpolationsRegridder,
                          interpolation_method = Interpolations.Constant(),
                          modis_lai_ncdata_path = nothing,
                          context = ClimaComms.context(surface_space))

A helper function which constructs the TimeVaryingInput object for Leaf Area
Index using MODIS LAI data; requires the
surface_space, the start and stop dates, and the earth_param_set.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported
by passing `interpolation_method = Interpolations.Linear()`.

If `modis_lai_ncdata_path` is provided, it will be used directly.
Otherwise, the path will be inferred from the start and stop dates.
"""
function prescribed_lai_modis(
    surface_space,
    start_date,
    stop_date;
    time_interpolation_method = LinearInterpolation(),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
    modis_lai_ncdata_path = nothing,
    context = ClimaComms.context(surface_space),
)
    modis_lai_ncdata_path =
        isnothing(modis_lai_ncdata_path) ?
        ClimaLand.Artifacts.modis_lai_multiyear_paths(;
            context,
            start_date,
            stop_date,
        ) : modis_lai_ncdata_path
    return TimeVaryingInput(
        modis_lai_ncdata_path,
        ["lai"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
end

"""
     prescribed_climatological_lai_modis(surface_space,
                                         time_interpolation_method = LinearInterpolation(PeriodicCalendar(Year(1), DateTime(2000))),
                                         regridder_type = :InterpolationsRegridder,
                                         interpolation_method = Interpolations.Constant(),
                                         context = ClimaComms.context(surface_space))

A helper function which constructs the TimeVaryingInput object for Leaf Area
Index using MODIS climatological LAI data; requires the
surface_space.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported
by passing interpolation_method = Interpolations.Linear().
"""
function prescribed_climatological_lai_modis(
    surface_space,
    time_interpolation_method = LinearInterpolation(
        PeriodicCalendar(Year(1), DateTime(2000)),
    ),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
    context = ClimaComms.context(surface_space),
)
    modis_lai_ncdata_path =
        ClimaLand.Artifacts.modis_lai_climatology_data_path(; context)
    return TimeVaryingInput(
        modis_lai_ncdata_path,
        ["lai"],
        surface_space;
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
end

"""
     modis_max_lai(surface_space,
                   regridder_type = :InterpolationsRegridder,
                   interpolation_method = Interpolations.Constant(),
                   context = ClimaComms.context(surface_space))

A helper function which constructs the SpaceVaryingInput object for the maximum
Leaf Area Index using MODIS LAI data; requires the surface_space.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported by passing
`interpolation_method = Interpolations.Linear()`.
"""
function modis_max_lai(
    surface_space,
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
    context = ClimaComms.context(surface_space),
)
    modis_max_lai_ncdata_path =
        ClimaLand.Artifacts.modis_max_lai_data_path(; context)
    return SpaceVaryingInput(
        modis_max_lai_ncdata_path,
        ["lai"],
        surface_space;
        regridder_type,
        regridder_kwargs = (; interpolation_method),
    )
end

"""
    AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for modeling the biomass (above ground - LAI, SAI, canopy
height) and below ground (rooting depth, RAI).
"""
abstract type AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(::AbstractBiomassModel) = :biomass

abstract type AbstractAreaIndexModel end

"""
   PrescribedAreaIndices{FS <: Union{AbstractFloat, ClimaCore.Fields.Field}, F <: AbstractTimeVaryingInput}

A struct containing the area indices of the plants at a specific site;
LAI varies in time, while SAI and RAI are fixed in time and can either be a
scalar (spatially uniform) or a ClimaCore Field (spatially varying).

$(DocStringExtensions.FIELDS)
"""
struct PrescribedAreaIndices{
    FS <: Union{AbstractFloat, ClimaCore.Fields.Field},
    F <: AbstractTimeVaryingInput,
} <: AbstractAreaIndexModel
    "A function of simulation time `t` giving the leaf area index (LAI; m2/m2)"
    LAI::F
    "The constant-in-time stem area index (SAI; m2/m2), scalar or Field"
    SAI::FS
    "The constant-in-time root area index (RAI; m2/m2), scalar or Field"
    RAI::FS
end

"""
    PrescribedAreaIndices(
        LAI::AbstractTimeVaryingInput,
        SAI,
        RAI,
    )

An outer constructor for setting the PrescribedAreaIndices given LAI, SAI, and
RAI. SAI and RAI may be scalars or ClimaCore Fields.
"""
function PrescribedAreaIndices(LAI::AbstractTimeVaryingInput, SAI, RAI)
    PrescribedAreaIndices{typeof(SAI), typeof(LAI)}(LAI, SAI, RAI)
end

"""
    struct PrescribedBiomassModel{FT, PSAI <: PrescribedAreaIndices, RDTH <: Union{FT, ClimaCore.Fields.Field}} <: AbstractBiomassModel{FT}

A prescribed biomass model where LAI, SAI, RAI, rooting depth, and height are prescribed.

In  global run with patches
of bare soil, you can "turn off" the canopy model (to get zero root extraction, zero absorption and
emission, zero transpiration and sensible heat flux from the canopy), by setting:
- LAI = SAI = RAI = 0.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedBiomassModel{
    FT,
    PSAI <: PrescribedAreaIndices,
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractBiomassModel{FT}
    "The plant area index model for LAI, SAI, RAI"
    plant_area_index::PSAI
    "Rooting depth parameter (m) - a characteristic depth below which 1/e of the root mass lies"
    rooting_depth::RDTH
    "Canopy height (m) - can be scalar (uniform) or spatially-varying Field"
    height::HTH
    function PrescribedBiomassModel{FT, PSAI, RDTH, HTH}(
        plant_area_index,
        rooting_depth,
        height,
    ) where {FT, PSAI, RDTH, HTH}
        new{FT, PSAI, RDTH, HTH}(plant_area_index, rooting_depth, height)
    end
end

"""
    PrescribedBiomassModel{FT}(;LAI::AbstractTimeVaryingInput,
                                SAI::FT,
                                RAI::FT,
                                rooting_depth,
                                height) where {FT}

An outer constructor to help set up the PrescribedBiomassModel from 
LAI, SAI, and RAI directly, instead of requiring the user to make the
area index object first; rooting_depth and height are also required.

Height can be either:
- A scalar FT value (uniform height across domain)
- A ClimaCore.Fields.Field (spatially-varying height)
"""
function PrescribedBiomassModel{FT}(;
    LAI::AbstractTimeVaryingInput,
    SAI,
    RAI,
    rooting_depth,
    height,
) where {FT}
    plant_area_index = PrescribedAreaIndices(LAI, SAI, RAI)
    args = (plant_area_index, rooting_depth, height)
    PrescribedBiomassModel{FT, typeof.(args)...}(args...)
end

ClimaLand.auxiliary_vars(model::PrescribedBiomassModel) = (:area_index,)
ClimaLand.auxiliary_types(model::PrescribedBiomassModel{FT}) where {FT} =
    (NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},)
ClimaLand.auxiliary_domain_names(::PrescribedBiomassModel) = (:surface,)

function clip(x::FT, threshold::FT) where {FT}
    x > threshold ? x : FT(0)
end

"""
    update_biomass!(
        p,
        Y,
        t,
        component::PrescribedBiomassModel{FT},
        canopy,
    ) where {FT}

Sets the area indices pertaining to their values at time t.

Note that we clip all values of LAI below 0.05 to zero, and SAI to zero
where the plant area index is below 0.05 ([`clip_absent_plant_area!`](@ref)).
This is because we currently run into issues when LAI is
of order eps(FT) in the SW radiation code.
Please see Issue #644
or PR #645 for details.
For now, this clipping is similar to what CLM and NOAH MP do.
"""
function update_biomass!(
    p,
    Y,
    t,
    component::PrescribedBiomassModel{FT},
    canopy,
) where {FT}
    (; LAI, SAI, RAI) = component.plant_area_index
    evaluate!(p.canopy.biomass.area_index.leaf, LAI, t)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    @. p.canopy.biomass.area_index.stem = SAI
    @. p.canopy.biomass.area_index.root = RAI
    clip_absent_plant_area!(p.canopy.biomass.area_index)
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
end

"""
    clip_absent_plant_area!(area_index)

Sets the stem area index in `area_index` to zero where the plants are absent
([`plants_present`](@ref)), after the leaf area index has been clipped below
0.05. An absent canopy then neither absorbs radiation nor stores heat, so the
relaxation of its temperature to the air temperature
(`make_compute_imp_tendency(::BigLeafEnergyModel, canopy)`) changes no
energy.
"""
function clip_absent_plant_area!(area_index)
    @. area_index.stem = ifelse(
        plants_present(area_index.leaf, area_index.stem),
        area_index.stem,
        zero(area_index.stem),
    )
    return nothing
end

"""
    mask_biomass!(p, prognostic_land_components)

Default method of setting LAI/RAI/SAI to zero where there
cannot be canopy; does nothing.

Currently, this is only does something when a lake model
is included in integrated models, as we cannot have 
vegetation over a lake, and the lake masks may not be consistent with
the biomass model.
"""
mask_biomass!(p, prognostic_land_components) = nothing

"""
    root_distribution(z::FT, rooting_depth::FT)

Computes value of rooting probability density function at `z`.

The rooting probability density function is derived from the
cumulative distribution function F(z) = 1 - β^(100z), which is described
by Equation 2.23 of
Bonan, "Climate Change and Terrestrial Ecosystem Modeling", 2019 Cambridge University Press.
This probability distribution function is equivalent to the derivative of the
cumulative distribution function with respect to z,
where `rooting_depth` replaces (-1)/(100ln(β)) and z is expected to be negative.
"""
function root_distribution(z::FT, rooting_depth::FT) where {FT <: AbstractFloat}
    return (1 / rooting_depth) * exp(z / rooting_depth) # 1/m
end

"""
    area_index_from_max_lai(max_LAI, tree_share, coeff_tree, coeff_grass)

Stem or root area index (m^2 m^-2) proportional to the maximum LAI `max_LAI`, with
the area index per unit of maximum LAI weighted between that of trees (`coeff_tree`)
and grasses (`coeff_grass`) by the tree share `tree_share`.
"""
area_index_from_max_lai(max_LAI, tree_share, coeff_tree, coeff_grass) =
    (tree_share * coeff_tree + (1 - tree_share) * coeff_grass) * max_LAI

"""
    AreaIndexFromMaxLAI{FT}(tree, grass)
    AreaIndexFromMaxLAI{FT}(toml_dict, name)

A stem or root area index of the `ZhouOptimalLAIModel` that follows the model's own
seasonal maximum LAI (`LAI_max`) and the tree share that sets the unit cost of leaves
(`area_index_from_max_lai`): `tree` and `grass` are the area index per unit of
`LAI_max` of trees and grasses. From `toml_dict`, they are the parameters
`"\$(name)_coeff_tree"` and `"\$(name)_coeff_grass"`, with `name` `"SAI"` or `"RAI"`.
"""
struct AreaIndexFromMaxLAI{FT <: AbstractFloat}
    "Area index per unit of maximum LAI of trees"
    tree::FT
    "Area index per unit of maximum LAI of grasses"
    grass::FT
end

AreaIndexFromMaxLAI{FT}(
    toml_dict::CP.ParamDict,
    name,
) where {FT <: AbstractFloat} = AreaIndexFromMaxLAI{FT}(
    toml_dict["$(name)_coeff_tree"],
    toml_dict["$(name)_coeff_grass"],
)

#####################################################################
# ZhouOptimalLAIModel - Optimal LAI model based on Zhou et al. (2025)
#####################################################################

"""
    ZhouOptimalLAIModel{FT, OLPT <: OptimalLAIParameters{FT}, GD, RDTH, HTH} <: AbstractBiomassModel{FT}

An implementation of the optimal LAI model from Zhou et al. (2025) as a biomass model.

This model computes LAI dynamically based on optimality principles, balancing energy and
water constraints. LAI is prognostic, in `Y.canopy.biomass.LAI`, and is mirrored into
`p.canopy.biomass.area_index.leaf` for the rest of the canopy, consistent with
`PrescribedBiomassModel`.

# Fields
- `parameters`: Required parameters for the optimal LAI model
- `SAI`: Stem area index (m^2 m^-2): prescribed (scalar or Field), or
  `AreaIndexFromMaxLAI` to follow the model's seasonal maximum LAI and tree share
- `RAI`: Root area index (m^2 m^-2), prescribed or `AreaIndexFromMaxLAI` like `SAI`
- `rooting_depth`: Rooting depth parameter (m) - a characteristic depth below which 1/e of the root mass lies
- `height`: Canopy height (m) - can be scalar (uniform) or spatially-varying Field
- `tree_share`: Tree share of the vegetation (dimensionless), which sets the unit cost of
  leaves (`leaf_cost`): prescribed (scalar or Field), or `ClimateTreeShare()` to
  compute it from the simulated climate (`climate_tree_share`)

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125
"""
struct ZhouOptimalLAIModel{
    FT,
    OLPT <: OptimalLAIParameters{FT},
    SA <: Union{FT, ClimaCore.Fields.Field, AreaIndexFromMaxLAI{FT}},
    RA <: Union{FT, ClimaCore.Fields.Field, AreaIndexFromMaxLAI{FT}},
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
    TSH <: Union{FT, ClimaCore.Fields.Field, ClimateTreeShare},
    T,
} <: AbstractBiomassModel{FT}
    "Required parameters for the optimal LAI model"
    parameters::OLPT
    "Stem area index (m^2 m^-2): a scalar or Field, or `AreaIndexFromMaxLAI`"
    SAI::SA
    "Root area index (m^2 m^-2): a scalar or Field, or `AreaIndexFromMaxLAI`"
    RAI::RA
    "Rooting depth parameter (m)"
    rooting_depth::RDTH
    "Canopy height (m) - can be scalar (uniform) or spatially-varying Field"
    height::HTH
    "Tree share of the vegetation (dimensionless): a scalar or Field, or `ClimateTreeShare()`"
    tree_share::TSH
    "Time integrated prognostic vars"
    time_integrated_vars::T
end

Base.eltype(::ZhouOptimalLAIModel{FT}) where {FT} = FT

"""
    ZhouOptimalLAIModel{FT}(
        parameters::OptimalLAIParameters{FT};
        SAI,
        RAI,
        rooting_depth,
        height,
        tree_share,
    ) where {FT <: AbstractFloat}

Outer constructor for the ZhouOptimalLAIModel struct.

# Arguments
- `parameters`: OptimalLAIParameters for the model
- `SAI`: Stem area index (m^2 m^-2); scalar or spatially-varying Field, or
  `AreaIndexFromMaxLAI` to follow the model's seasonal maximum LAI and tree share
- `RAI`: Root area index (m^2 m^-2); scalar or spatially-varying Field, or
  `AreaIndexFromMaxLAI`
- `rooting_depth`: Rooting depth parameter (m)
- `height`: Canopy height (m) - can be scalar or spatially-varying Field
- `tree_share`: Tree share of the vegetation (dimensionless), which weights the unit
  cost of leaves between `z_tree` and `z_grass`; scalar or spatially-varying Field, or
  `ClimateTreeShare()` to compute it from the simulated climate

Declares the prognostic time integrated variables: the 1-day potential-GPP total
`A0_daily`, the 1-year totals `A0_annual` and `precip_annual` as `RunningSum`s of
the instantaneous rate, and `LAI` as a `RunningMean` relaxing toward the
instantaneous steady-state target `L_opt`:

    dA0_daily/dt      = (day·A0 - A0_daily) / τ_day,             τ_day  = 3 days,
    dA0_annual/dt     = (year·A0 - A0_annual) / τ_long,          τ_long = tau_long_term,
    dprecip_annual/dt = (year·P_inst  - precip_annual) / τ_long,
    dLAI/dt           = (L_opt - LAI) / τ_LAI,                   τ_LAI  = 1 day / α.

Six further 1-year `RunningSum`s carry the climate the LAI formulas respond to:
`PET_annual` (with `precip_annual`, the aridity index behind `f0`), `VPDgs_annual`
(the VPD summed while the air is above freezing, which with `growing_days` gives the
VPD of the growing season, used where it has no moist season), `growing_days` (the
growing-season length `GSL`),
`A0c3_annual`/`A0c4_annual`, the per-pathway
potential GPP the C3/C4 competition compares, and `GPPc3_annual`, the C3 potential
GPP scaled by the realized fAPAR and soil-moisture stress `βm`, from which the
competition estimates tree cover.

Each `RunningSum` holds a total over its own window (1 day, 1 year) whatever the
smoothing timescale τ_long: only the smoothing changes with τ_long, not the magnitude,
as the LAI_max/steady-state formulas require.

Nine more variables carry the climate of the moist growing season, of the leaf cost
and of the tree share: a degree-day snow store (`snow_store`), which holds the
precipitation that falls below freezing until it melts; the 30-day totals of the water
input, rain and snowmelt (`water_30d`), and of the PET (`PET_30d`); the yearly sums of
the VPD (`VPD_moist_annual`) and of the days (`moist_days`) of the moist growing
season, when the air is above freezing and the 30-day water input is at least half the
PET; the yearly sums of the degree-days (`degree_days`), of the days (`warm_days`) and
of the maintenance load (`maintenance_days`, see `leaf_cost_scale`) above freezing; and
their `age`, the time since the start. The yearly ones average all of their history
until it reaches τ_long, so they hardly depend on their initial values.
"""
function ZhouOptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT};
    SAI,
    RAI,
    rooting_depth,
    height,
    tree_share,
) where {FT <: AbstractFloat}
    seconds_per_day = IP.day(IP.InsolationParameters(FT))
    tau_long_term = parameters.tau_long_term
    base_tivs = (
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :A0_daily,
            reduction = ClimaLand.RunningSum(
                seconds_per_day,
                3 * seconds_per_day,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :A0_annual,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :precip_annual,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :PET_annual,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :VPDgs_annual,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :growing_days,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :A0c3_annual,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :A0c4_annual,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :GPPc3_annual,
            reduction = ClimaLand.RunningSum(
                365 * seconds_per_day,
                tau_long_term,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :LAI,
            reduction = ClimaLand.RunningMean(
                seconds_per_day / parameters.alpha,
            ),
        ),
    )
    month = 30 * seconds_per_day
    year = 365 * seconds_per_day
    climate_tivs = (
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :snow_store,
            reduction = ClimaLand.TimeIntegral(),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :water_30d,
            reduction = ClimaLand.RunningSum(month, month),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :PET_30d,
            reduction = ClimaLand.RunningSum(month, month),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :VPD_moist_annual,
            reduction = ClimaLand.RunningSum(year, tau_long_term),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :moist_days,
            reduction = ClimaLand.RunningSum(year, tau_long_term),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :degree_days,
            reduction = ClimaLand.RunningSum(year, tau_long_term),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :warm_days,
            reduction = ClimaLand.RunningSum(year, tau_long_term),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :maintenance_days,
            reduction = ClimaLand.RunningSum(year, tau_long_term),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :age,
            reduction = ClimaLand.TimeIntegral(),
        ),
    )
    tiv = ClimaLand.time_integrated_variables(base_tivs..., climate_tivs...)
    return ZhouOptimalLAIModel{
        FT,
        typeof(parameters),
        typeof(SAI),
        typeof(RAI),
        typeof(rooting_depth),
        typeof(height),
        typeof(tree_share),
        typeof(tiv),
    }(
        parameters,
        SAI,
        RAI,
        rooting_depth,
        height,
        tree_share,
        tiv,
    )
end

"""
    ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel)
    ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel)
    ClimaLand.auxiliary_domain_names(model::ZhouOptimalLAIModel)

Defines the auxiliary variables for the ZhouOptimalLAIModel:
- `area_index`: NamedTuple{(:root, :stem, :leaf)} containing area indices (m^2 m^-2)
- `OptVars.A0, OptVars.χ`: instantaneous potential GPP (mol CO2 m^-2 s^-1) and ci/ca ratio at the growing-season VPD (`compute_A0_and_χ`), from the PModel with the model's own unit cost ratios and cost of Jmax (`β_c3`, `β_c4`, `cstar`); `LAI_max` uses the ci/ca ratio of the growing-season temperature instead
- `OptVars.A0_c3, OptVars.A0_c4`: the same potential GPP for a pure-C3 and a pure-C4 canopy, which the C3/C4 competition compares
- `L_opt`: steady-state LAI of Zhou et al. (2025), capped by `LAI_max`, toward which the prognostic `LAI` relaxes
- `GSL`: growing season length (days), the trailing-year count of days above freezing
- `vpd_gs`: mean VPD of the moist growing season (Pa), for the water-limitation term of LAI_max (`moist_season_vpd`)
- `f0`: fraction of precipitation available for transpiration (dimensionless), from the aridity index
- `composition.tree, composition.c3_grass, composition.c4_grass`: shares of the canopy
  from the C3/C4 competition (dimensionless, summing to one); see
  `canopy_composition_from_competition`
- `LAI_max`: seasonal maximum LAI (m^2 m^-2) of Zhou et al. (2025), which sets and caps
  the steady-state LAI, the leaf area trees retain and the stem and root area indices
  of `AreaIndexFromMaxLAI`

`GSL`, `vpd_gs`, `f0`, `composition` and `LAI_max` are derived in `update_biomass!`
from the trailing totals in `Y`.
"""
ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel) =
    (:area_index, :OptVars, :L_opt, :GSL, :vpd_gs, :f0, :composition, :LAI_max)
ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel{FT}) where {FT} = (
    NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    NamedTuple{(:A0, :A0_c3, :A0_c4, :χ), NTuple{4, FT}},
    FT,
    FT,
    FT,
    FT,
    NamedTuple{(:tree, :c3_grass, :c4_grass), NTuple{3, FT}},
    FT,
)
ClimaLand.auxiliary_domain_names(::ZhouOptimalLAIModel) = (
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
)

ClimaLand.prognostic_vars(m::ZhouOptimalLAIModel) =
    ClimaLand.time_integrated_prognostic_vars(m.time_integrated_vars)
ClimaLand.prognostic_types(m::ZhouOptimalLAIModel) =
    ClimaLand.time_integrated_prognostic_types(m.time_integrated_vars)
ClimaLand.prognostic_domain_names(m::ZhouOptimalLAIModel) =
    ClimaLand.time_integrated_prognostic_domain_names(m.time_integrated_vars)

"""
    update_biomass!(
        p,
        Y,
        t,
        component::ZhouOptimalLAIModel{FT},
        canopy,
    ) where {FT}

Updates the optimal-LAI cache from the prognostic state in `Y`: derives `f0`, `vpd_gs`
and `GSL` from the trailing climate totals; sets the canopy `composition` (tree, C3
grass, C4 grass shares) from the C3/C4 competition on the trailing per-pathway
potential GPP (`A0c3_annual`, `A0c4_annual`) and realized C3 GPP (`GPPc3_annual`), or
the climate tree share; sets the seasonal maximum LAI `LAI_max`
(`set_max_leaf_area_index!`); sets the leaf area index used by the rest of the canopy
from the prognostic `LAI` and the canopy trees retain (`set_leaf_area_index!`),
clipped below 0.05; and sets SAI and RAI (`set_area_index!`), zeroing SAI where plants
are absent (`clip_absent_plant_area!`) and all area indices where a lake is present
(`mask_biomass!`). This runs first in `update_aux`, before radiative transfer and
photosynthesis read the area index and C3 fraction.
"""
function update_biomass!(
    p,
    Y,
    t,
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    (; SAI, RAI, parameters) = component
    @. p.canopy.biomass.f0 = f0_from_aridity(
        Y.canopy.biomass.PET_annual,
        Y.canopy.biomass.precip_annual,
        parameters.f0_max,
    )
    seconds_per_day = IP.day(IP.InsolationParameters(FT))
    b = Y.canopy.biomass
    @. p.canopy.biomass.vpd_gs = moist_season_vpd(
        b.VPD_moist_annual,
        b.moist_days,
        b.VPDgs_annual,
        b.growing_days,
        seconds_per_day,
    )
    @. p.canopy.biomass.GSL = b.growing_days
    update_composition!(p, Y, component.tree_share, component, canopy)
    set_max_leaf_area_index!(p, Y, component, canopy)
    set_leaf_area_index!(p, Y, component)
    # Apply clipping to LAI (same as PrescribedBiomassModel)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    tree = leaf_cost_tree_share(component.tree_share, p)
    (; area_index, LAI_max) = p.canopy.biomass
    set_area_index!(area_index.stem, SAI, LAI_max, tree)
    set_area_index!(area_index.root, RAI, LAI_max, tree)
    clip_absent_plant_area!(area_index)
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
end

# Memory (s) of a growing running sum at its start: its initial value weighs as much as
# this much of its history.
const RUNNING_SUM_START_MEMORY = 86400

# Melt of the snow store of the optimal-LAI climate (m of water per K per day above
# freezing), and the shortest time (s) over which it can empty
const SNOWMELT_DEGREE_DAY_FACTOR = 3e-3
const SNOWMELT_MIN_TIMESCALE = 3600

"""
    growing_running_sum_tendency(f, X, age, reduction::ClimaLand.RunningSum)

Tendency of a running sum `X` of `f` whose memory grows with its `age` (s) until it
reaches the timescale of `reduction`: `dX/dt = (f τ − X)/min(age + τ₀, τ_long)`, with
τ₀ = 1 day. While the memory grows, `X` is the mean of all of its history (and of its
initial value, with weight τ₀) scaled to the window τ, so it hardly depends on its
initial value.
"""
growing_running_sum_tendency(f, X, age, reduction::ClimaLand.RunningSum) =
    (f * reduction.τ - X) /
    min(age + oftype(age, RUNNING_SUM_START_MEMORY), reduction.τ_long)

"""
    moist_season_vpd(VPD_moist_annual, moist_days, VPDgs_annual, growing_days, seconds_per_day)

Mean VPD of the moist growing season (Pa), when the air is above freezing and the
30-day water input (rain and snowmelt) is at least half the PET: the VPD the water limit of LAI_max
uses. Where the moist season is shorter than a month (deserts), it blends into the
VPD of the whole growing season (above freezing).
"""
function moist_season_vpd(
    VPD_moist_annual::FT,
    moist_days::FT,
    VPDgs_annual::FT,
    growing_days::FT,
    seconds_per_day,
) where {FT}
    vpd_moist = VPD_moist_annual / max(moist_days * seconds_per_day, eps(FT))
    vpd_warm = VPDgs_annual / max(growing_days * seconds_per_day, eps(FT))
    w = min(moist_days / 30, one(FT))
    return w * vpd_moist + (1 - w) * vpd_warm
end

"""
    leaf_cost_tree_share(tree_share, p)

Tree share that sets the leaf cost `z`: the prescribed one, or the climate tree share
of the canopy composition with `ClimateTreeShare()`.
"""
leaf_cost_tree_share(tree_share, p) = tree_share
leaf_cost_tree_share(::ClimateTreeShare, p) = p.canopy.biomass.composition.tree

"""
    set_max_leaf_area_index!(p, Y, component::ZhouOptimalLAIModel, canopy)

Sets the seasonal maximum LAI `p.canopy.biomass.LAI_max` from the trailing totals in
`Y`, with the χ of the growing-season temperature and VPD.
"""
function set_max_leaf_area_index!(
    p,
    Y,
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    parameters = component.parameters
    lai_pmodel_parameters = optimal_lai_pmodel_parameters(
        canopy.photosynthesis.parameters,
        parameters,
    )
    constants = canopy.photosynthesis.constants
    T_freeze = LP.T_freeze(canopy.earth_param_set)
    b = Y.canopy.biomass
    tree = leaf_cost_tree_share(component.tree_share, p)
    T_growing = @. lazy(b.degree_days / max(b.warm_days, eps(FT)))
    fractional_c3 = @. lazy(1 - p.canopy.biomass.composition.c4_grass)
    χ_growing = @. lazy(
        fractional_c3 * optimal_chi(
            T_freeze + T_growing,
            p.drivers.P,
            p.drivers.c_co2,
            p.canopy.biomass.vpd_gs,
            lai_pmodel_parameters.β_c3,
            constants,
        ) +
        (1 - fractional_c3) * optimal_chi(
            T_freeze + T_growing,
            p.drivers.P,
            p.drivers.c_co2,
            p.canopy.biomass.vpd_gs,
            lai_pmodel_parameters.β_c4,
            constants,
        ),
    )
    @. p.canopy.biomass.LAI_max = compute_L_max(
        b.A0_annual,
        parameters.k,
        leaf_cost(tree, parameters.z_tree, parameters.z_grass) *
        leaf_cost_scale(b.maintenance_days, parameters.maintenance_share),
        b.precip_annual,
        p.canopy.biomass.f0,
        p.drivers.c_co2 * p.drivers.P,
        χ_growing,
        p.canopy.biomass.vpd_gs,
    )
end

"""
    set_leaf_area_index!(p, Y, component::ZhouOptimalLAIModel)

Sets the leaf area index from the prognostic `LAI`, which follows the seasonal optimum,
and the canopy trees keep through the unfavourable season: the tree share keeps at
least `tree_retention` of `LAI_max`.
"""
function set_leaf_area_index!(p, Y, component::ZhouOptimalLAIModel)
    tree = leaf_cost_tree_share(component.tree_share, p)
    LAI = Y.canopy.biomass.LAI
    retained_LAI =
        @. lazy(component.parameters.tree_retention * p.canopy.biomass.LAI_max)
    @. p.canopy.biomass.area_index.leaf =
        LAI + tree * max(retained_LAI - LAI, 0)
end

"""
    set_area_index!(area_index, value, LAI_max, tree_share)

Sets a stem or root `area_index` to its prescribed `value` (scalar or Field), or, for
an `AreaIndexFromMaxLAI`, to its multiple of the seasonal maximum LAI `LAI_max` for
the tree share `tree_share`.
"""
set_area_index!(area_index, value, LAI_max, tree_share) = (area_index .= value)
set_area_index!(area_index, value::AreaIndexFromMaxLAI, LAI_max, tree_share) =
    @. area_index =
        area_index_from_max_lai(LAI_max, tree_share, value.tree, value.grass)

"""
    update_composition!(p, Y, tree_share, component::ZhouOptimalLAIModel, canopy)

Sets the canopy `composition` from the C3/C4 competition. With a prescribed tree share,
the competition's tree share is estimated from the realized C3 GPP
(`canopy_composition_from_competition`); with `ClimateTreeShare()`, it is the climate
tree share (`climate_tree_share`) that also sets the leaf cost, from the LAI_max of a
tree canopy, the number of dry months of the growing season (its warm days that are not
moist) and the growing-season temperature.
"""
function update_composition!(p, Y, tree_share, component, canopy)
    @. p.canopy.biomass.composition = canopy_composition_from_competition(
        Y.canopy.biomass.A0c3_annual,
        Y.canopy.biomass.A0c4_annual,
        Y.canopy.biomass.GPPc3_annual,
        Y.canopy.biomass.growing_days,
        canopy.photosynthesis.constants.Mc,
        component.parameters,
    )
end

function update_composition!(
    p,
    Y,
    ::ClimateTreeShare,
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    parameters = component.parameters
    lai_pmodel_parameters = optimal_lai_pmodel_parameters(
        canopy.photosynthesis.parameters,
        parameters,
    )
    constants = canopy.photosynthesis.constants
    T_freeze = LP.T_freeze(canopy.earth_param_set)
    b = Y.canopy.biomass
    # mean air temperature above freezing (°C)
    T_growing = @. lazy(b.degree_days / max(b.warm_days, eps(FT)))
    χ_growing = @. lazy(
        optimal_chi(
            T_freeze + T_growing,
            p.drivers.P,
            p.drivers.c_co2,
            p.canopy.biomass.vpd_gs,
            lai_pmodel_parameters.β_c3,
            constants,
        ),
    )
    L_tree = @. lazy(
        compute_L_max(
            b.A0c3_annual,
            parameters.k,
            parameters.z_tree * leaf_cost_scale(
                b.maintenance_days,
                parameters.maintenance_share,
            ),
            b.precip_annual,
            p.canopy.biomass.f0,
            p.drivers.c_co2 * p.drivers.P,
            χ_growing,
            p.canopy.biomass.vpd_gs,
        ),
    )
    @. p.canopy.biomass.composition = canopy_composition(
        climate_tree_share(
            L_tree,
            max(b.warm_days - b.moist_days, 0) * 12 / 365,
            T_growing,
            parameters,
        ),
        open_canopy_c4_share(b.A0c3_annual, b.A0c4_annual, parameters),
    )
end

"""
    get_fractional_c3(p, canopy)
    get_fractional_c3(p, biomass::AbstractBiomassModel, photosynthesis)

C3 fraction of the canopy (1 = all C3), the weight photosynthesis blends its C3 and
C4 pathways with. A biomass model that predicts the canopy composition
(`ZhouOptimalLAIModel`) sets it; otherwise it is the photosynthesis model's static
value.
"""
get_fractional_c3(p, canopy) =
    get_fractional_c3(p, canopy.biomass, canopy.photosynthesis)
get_fractional_c3(p, ::AbstractBiomassModel, photosynthesis) =
    static_fractional_c3(photosynthesis)
get_fractional_c3(p, ::ZhouOptimalLAIModel, photosynthesis) =
    @. lazy(1 - p.canopy.biomass.composition.c4_grass)

"""
    ClimaLand.make_compute_exp_tendency(component::ZhouOptimalLAIModel, canopy)

Advances the optimal-LAI model's time-integrated variables.
"""
function ClimaLand.make_compute_exp_tendency(
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    ρ_m_liq = LP.ρ_m_liq(canopy.earth_param_set)  # mol H2O m^-3 (precip volume flux → molar flux)
    tivs = component.time_integrated_vars
    seconds_per_day = IP.day(IP.InsolationParameters(FT))
    # mol H2O m^-2 s^-1 K^-1
    snowmelt_factor = FT(SNOWMELT_DEGREE_DAY_FACTOR) * ρ_m_liq / seconds_per_day
    earth_param_set = canopy.earth_param_set
    σ = LP.Stefan(earth_param_set)
    M_w = LP.molar_mass_water(earth_param_set)  # kg mol^-1
    T_freeze = LP.T_freeze(earth_param_set)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    ϵ_sfc = canopy.radiative_transfer.parameters.ϵ_canopy
    parameters = component.parameters
    lai_pmodel_parameters = optimal_lai_pmodel_parameters(
        canopy.photosynthesis.parameters,
        parameters,
    )
    pmodel_constants = canopy.photosynthesis.constants
    function compute_exp_tendency!(dY, Y, p, t)
        fractional_c3 = get_fractional_c3(p, canopy)
        # Supersaturated forcing gives a negative VPD.
        VPD = @. lazy(
            max(
                Thermodynamics.vapor_pressure_deficit(
                    thermo_params,
                    p.drivers.T,
                    p.drivers.P,
                    p.drivers.q,
                ),
                zero(FT),
            ),
        )

        # A0 is a potential GPP, so βm = 1: water limitation enters once, through
        # the f0·P/A0 term of LAI_max.
        @. p.canopy.biomass.OptVars = compute_A0_and_χ(
            fractional_c3,
            lai_pmodel_parameters,
            pmodel_constants,
            earth_param_set,
            p.drivers.T,
            p.drivers.P,
            p.drivers.q,
            p.drivers.c_co2,
            compute_PPFD(
                p.canopy.radiative_transfer.par_d,
                canopy.radiative_transfer.parameters.λ_γ_PAR,
                pmodel_constants.lightspeed,
                pmodel_constants.planck_h,
                pmodel_constants.N_a,
            ),
            one(FT),
            p.canopy.biomass.vpd_gs,
        )

        (; GSL, LAI_max) = p.canopy.biomass
        @. p.canopy.biomass.L_opt = compute_steady_state_LAI(
            Y.canopy.biomass.A0_daily,
            compute_m(
                GSL,
                LAI_max,
                Y.canopy.biomass.A0_annual,
                parameters.sigma,
                parameters.k,
            ),
            parameters.k,
            LAI_max,
        )

        @. dY.canopy.biomass.A0_daily = apply_time_reduction(
            p.canopy.biomass.OptVars.A0,
            Y.canopy.biomass.A0_daily,
            tivs.A0_daily.reduction,
        )
        @. dY.canopy.biomass.A0_annual = apply_time_reduction(
            p.canopy.biomass.OptVars.A0,
            Y.canopy.biomass.A0_annual,
            tivs.A0_annual.reduction,
        )
        # P_liq/P_snow are negative-downward volume fluxes (m/s); negate for a positive total.
        precip = @. lazy(-(p.drivers.P_liq + p.drivers.P_snow) * ρ_m_liq)
        @. dY.canopy.biomass.precip_annual = apply_time_reduction(
            precip,
            Y.canopy.biomass.precip_annual,
            tivs.precip_annual.reduction,
        )
        PET = @. lazy(
            potential_evaporation(
                p.drivers.SW_d,
                p.drivers.LW_d,
                p.drivers.T,
                p.drivers.P,
                ϵ_sfc,
                σ,
                M_w,
                thermo_params,
            ),
        )
        # PET_annual / precip_annual is the aridity index behind f0.
        @. dY.canopy.biomass.PET_annual = apply_time_reduction(
            PET,
            Y.canopy.biomass.PET_annual,
            tivs.PET_annual.reduction,
        )
        # VPD while the air is above freezing; with growing_days, the VPD of the
        # growing season where it has no moist season.
        @. dY.canopy.biomass.VPDgs_annual = apply_time_reduction(
            ifelse(p.drivers.T > T_freeze, VPD, zero(FT)),
            Y.canopy.biomass.VPDgs_annual,
            tivs.VPDgs_annual.reduction,
        )
        # 1/day while air T is above freezing, so the yearly total is GSL in days.
        @. dY.canopy.biomass.growing_days = apply_time_reduction(
            ifelse(p.drivers.T > T_freeze, 1 / seconds_per_day, zero(FT)),
            Y.canopy.biomass.growing_days,
            tivs.growing_days.reduction,
        )
        @. dY.canopy.biomass.A0c3_annual = apply_time_reduction(
            p.canopy.biomass.OptVars.A0_c3,
            Y.canopy.biomass.A0c3_annual,
            tivs.A0c3_annual.reduction,
        )
        @. dY.canopy.biomass.A0c4_annual = apply_time_reduction(
            p.canopy.biomass.OptVars.A0_c4,
            Y.canopy.biomass.A0c4_annual,
            tivs.A0c4_annual.reduction,
        )
        # The tree-cover relation is fitted to annual realized GPP, so the C3
        # potential is scaled by the realized fAPAR and soil-moisture stress.
        @. dY.canopy.biomass.GPPc3_annual = apply_time_reduction(
            p.canopy.biomass.OptVars.A0_c3 *
            (1 - exp(-parameters.k * Y.canopy.biomass.LAI)) *
            p.canopy.soil_moisture_stress.βm,
            Y.canopy.biomass.GPPc3_annual,
            tivs.GPPc3_annual.reduction,
        )
        @. dY.canopy.biomass.LAI = apply_time_reduction(
            p.canopy.biomass.L_opt,
            Y.canopy.biomass.LAI,
            tivs.LAI.reduction,
        )
        b = Y.canopy.biomass
        # Precipitation below freezing accumulates as snow until it melts.
        frozen = @. lazy(p.drivers.T <= T_freeze)
        melt = @. lazy(
            ifelse(
                frozen,
                zero(FT),
                min(
                    snowmelt_factor * (p.drivers.T - T_freeze),
                    b.snow_store / SNOWMELT_MIN_TIMESCALE,
                ),
            ),
        )
        @. dY.canopy.biomass.snow_store = apply_time_reduction(
            ifelse(frozen, precip, zero(FT)) - melt,
            b.snow_store,
            tivs.snow_store.reduction,
        )
        @. dY.canopy.biomass.water_30d = apply_time_reduction(
            ifelse(frozen, zero(FT), precip) + melt,
            b.water_30d,
            tivs.water_30d.reduction,
        )
        @. dY.canopy.biomass.PET_30d =
            apply_time_reduction(PET, b.PET_30d, tivs.PET_30d.reduction)
        # Moist growing season: above freezing, and the last 30 days not dry.
        moist =
            @. lazy((p.drivers.T > T_freeze) & (b.water_30d >= b.PET_30d / 2))
        @. dY.canopy.biomass.VPD_moist_annual = growing_running_sum_tendency(
            ifelse(moist, VPD, zero(FT)),
            b.VPD_moist_annual,
            b.age,
            tivs.VPD_moist_annual.reduction,
        )
        # 1/day while moist (warm), so the yearly totals are in days.
        @. dY.canopy.biomass.moist_days = growing_running_sum_tendency(
            ifelse(moist, 1 / seconds_per_day, zero(FT)),
            b.moist_days,
            b.age,
            tivs.moist_days.reduction,
        )
        @. dY.canopy.biomass.degree_days = growing_running_sum_tendency(
            max(p.drivers.T - T_freeze, zero(FT)) / seconds_per_day,
            b.degree_days,
            b.age,
            tivs.degree_days.reduction,
        )
        @. dY.canopy.biomass.warm_days = growing_running_sum_tendency(
            ifelse(p.drivers.T > T_freeze, 1 / seconds_per_day, zero(FT)),
            b.warm_days,
            b.age,
            tivs.warm_days.reduction,
        )
        @. dY.canopy.biomass.maintenance_days = growing_running_sum_tendency(
            maintenance_rate(
                p.drivers.T,
                T_freeze,
                parameters.maintenance_q10,
            ) / seconds_per_day,
            b.maintenance_days,
            b.age,
            tivs.maintenance_days.reduction,
        )
        @. dY.canopy.biomass.age =
            apply_time_reduction(one(FT), b.age, tivs.age.reduction)
    end
    return compute_exp_tendency!
end

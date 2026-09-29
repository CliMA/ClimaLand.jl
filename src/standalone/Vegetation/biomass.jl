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
    PrognosticCarbonParameters,
    PrognosticCarbonModel,
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
    prescribed_lai_input(model::AbstractBiomassModel)

The prescribed LAI `TimeVaryingInput` of a biomass model, including one wrapped by
`PrognosticCarbonModel`.
"""
prescribed_lai_input(model::PrescribedBiomassModel) = model.plant_area_index.LAI

"""
    update_biomass!(
        p,
        Y,
        t,
        component::PrescribedBiomassModel{FT},
        canopy,
    ) where {FT}

Sets the area indices pertaining to their values at time t.

Note that we clip all values of LAI below 0.05 to zero.
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
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
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
- `SAI`: Prescribed stem area index (m^2 m^-2)
- `RAI`: Prescribed root area index (m^2 m^-2)
- `rooting_depth`: Rooting depth parameter (m) - a characteristic depth below which 1/e of the root mass lies
- `height`: Canopy height (m) - can be scalar (uniform) or spatially-varying Field

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125
"""
struct ZhouOptimalLAIModel{
    FT,
    OLPT <: OptimalLAIParameters{FT},
    FS <: Union{FT, ClimaCore.Fields.Field},
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
    T,
} <: AbstractBiomassModel{FT}
    "Required parameters for the optimal LAI model"
    parameters::OLPT
    "Prescribed stem area index (m^2 m^-2), scalar or Field"
    SAI::FS
    "Prescribed root area index (m^2 m^-2), scalar or Field"
    RAI::FS
    "Rooting depth parameter (m)"
    rooting_depth::RDTH
    "Canopy height (m) - can be scalar (uniform) or spatially-varying Field"
    height::HTH
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
    ) where {FT <: AbstractFloat}

Outer constructor for the ZhouOptimalLAIModel struct.

# Arguments
- `parameters`: OptimalLAIParameters for the model
- `SAI`: Prescribed stem area index (m^2 m^-2); scalar or spatially-varying Field
- `RAI`: Prescribed root area index (m^2 m^-2); scalar or spatially-varying Field
- `rooting_depth`: Rooting depth parameter (m)
- `height`: Canopy height (m) - can be scalar or spatially-varying Field

Declares the prognostic time integrated variables: the 1-day potential-GPP total
`A0_daily`, the 1-year totals `A0_annual` and `precip_annual` as `RunningSum`s of
the instantaneous rate, and `LAI` as a `RunningMean` relaxing toward the
instantaneous steady-state target `L_opt`:

    dA0_daily/dt      = (day·A0 - A0_daily) / τ_day,             τ_day  = 3 days,
    dA0_annual/dt     = (year·A0 - A0_annual) / τ_long,          τ_long = tau_long_term,
    dprecip_annual/dt = (year·P_inst  - precip_annual) / τ_long,
    dLAI/dt           = (L_opt - LAI) / τ_LAI,                   τ_LAI  = 1 day / α.

Six further 1-year `RunningSum`s carry the climate the LAI formulas respond to:
`PET_annual` (with `precip_annual`, the aridity index behind `f0`), `VPDA0_annual`
(with `A0_annual`, the A0-weighted growing-season VPD `vpd_gs`), `growing_days`
(the growing-season length `GSL`), `A0c3_annual`/`A0c4_annual`, the per-pathway
potential GPP the C3/C4 competition compares, and `GPPc3_annual`, the C3 potential
GPP scaled by the realized fAPAR, from which the competition estimates tree cover.

Each `RunningSum` holds a total over its own window (1 day, 1 year) whatever the
smoothing timescale τ_long: only the smoothing changes with τ_long, not the magnitude,
as the LAI_max/steady-state formulas require.
"""
function ZhouOptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT};
    SAI,
    RAI,
    rooting_depth,
    height,
) where {FT <: AbstractFloat}
    seconds_per_day = IP.day(IP.InsolationParameters(FT))
    tau_long_term = parameters.tau_long_term
    tiv = ClimaLand.time_integrated_variables(
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
            name = :VPDA0_annual,
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
    return ZhouOptimalLAIModel{
        FT,
        typeof(parameters),
        typeof(SAI),
        typeof(rooting_depth),
        typeof(height),
        typeof(tiv),
    }(
        parameters,
        SAI,
        RAI,
        rooting_depth,
        height,
        tiv,
    )
end

"""
    ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel)
    ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel)
    ClimaLand.auxiliary_domain_names(model::ZhouOptimalLAIModel)

Defines the auxiliary variables for the ZhouOptimalLAIModel:
- `area_index`: NamedTuple{(:root, :stem, :leaf)} containing area indices (m^2 m^-2)
- `OptVars.A0, OptVars.χ`: instantaneous potential GPP (mol CO2 m^-2 s^-1) and ci/ca ratio computed using the optimal values from the PModel
- `OptVars.A0_c3, OptVars.A0_c4`: the same potential GPP for a pure-C3 and a pure-C4 canopy, which the C3/C4 competition compares
- `L_opt`: Optimal LAI predicted by Zhou et al.
- `GSL`: growing season length (days), the trailing-year count of days above freezing
- `vpd_gs`: A0-weighted mean VPD (Pa), for the water-limitation term of LAI_max
- `f0`: fraction of precipitation available for transpiration (dimensionless), from the aridity index
- `composition.tree, composition.c3_grass, composition.c4_grass`: shares of the canopy
  from the C3/C4 competition (dimensionless, summing to one); see
  `canopy_composition_from_competition`

`GSL`, `vpd_gs`, `f0` and `composition` are derived in `update_biomass!` from the
trailing totals in `Y`.
"""
ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel) =
    (:area_index, :OptVars, :L_opt, :GSL, :vpd_gs, :f0, :composition)
ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel{FT}) where {FT} = (
    NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    NamedTuple{(:A0, :A0_c3, :A0_c4, :χ), NTuple{4, FT}},
    FT,
    FT,
    FT,
    FT,
    NamedTuple{(:tree, :c3_grass, :c4_grass), NTuple{3, FT}},
)
ClimaLand.auxiliary_domain_names(::ZhouOptimalLAIModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface, :surface)

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

Updates the optimal-LAI cache from the prognostic state in `Y`: sets SAI and RAI to
their prescribed values; derives `f0`, `vpd_gs` and `GSL` from the trailing climate
totals; sets the canopy `composition` (tree, C3 grass, C4 grass shares) from the
C3/C4 competition on the trailing per-pathway potential GPP (`A0c3_annual`,
`A0c4_annual`) and realized C3 GPP (`GPPc3_annual`); and sets the leaf area index
used by the rest of the canopy from the prognostic `LAI`, clipped below 0.05 and
zeroed where a lake is present (`mask_biomass!`). This runs first in `update_aux`,
before radiative transfer and photosynthesis read the area index and C3 fraction.
"""
function update_biomass!(
    p,
    Y,
    t,
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    (; SAI, RAI, parameters) = component
    @. p.canopy.biomass.area_index.stem = SAI
    @. p.canopy.biomass.area_index.root = RAI
    @. p.canopy.biomass.f0 = f0_from_aridity(
        Y.canopy.biomass.PET_annual,
        Y.canopy.biomass.precip_annual,
        parameters.f0_max,
    )
    @. p.canopy.biomass.vpd_gs =
        Y.canopy.biomass.VPDA0_annual / max(Y.canopy.biomass.A0_annual, eps(FT))
    @. p.canopy.biomass.GSL = Y.canopy.biomass.growing_days
    @. p.canopy.biomass.composition = canopy_composition_from_competition(
        Y.canopy.biomass.A0c3_annual,
        Y.canopy.biomass.A0c4_annual,
        Y.canopy.biomass.GPPc3_annual,
        canopy.photosynthesis.constants.Mc,
        parameters,
    )
    @. p.canopy.biomass.area_index.leaf = Y.canopy.biomass.LAI
    # Apply clipping to LAI (same as PrescribedBiomassModel)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
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

Advances the optimal-LAI model's ten time-integrated variables.
"""
function ClimaLand.make_compute_exp_tendency(
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    ρ_m_liq = LP.ρ_m_liq(canopy.earth_param_set)  # mol H2O m^-3 (precip volume flux → molar flux)
    tivs = component.time_integrated_vars
    seconds_per_day = IP.day(IP.InsolationParameters(FT))
    earth_param_set = canopy.earth_param_set
    σ = LP.Stefan(earth_param_set)
    M_w = LP.molar_mass_water(earth_param_set)  # kg mol^-1
    T_freeze = LP.T_freeze(earth_param_set)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    h_atmos = canopy.boundary_conditions.atmos.h
    ϵ_sfc = canopy.radiative_transfer.parameters.ϵ_canopy
    parameters = component.parameters
    pmodel_parameters = canopy.photosynthesis.parameters
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
            pmodel_parameters,
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

        @. p.canopy.biomass.L_opt = compute_L_steady_target(
            Y.canopy.biomass.A0_daily,
            parameters.k,
            Y.canopy.biomass.A0_annual,
            parameters.z,
            p.canopy.biomass.GSL,
            parameters.sigma,
            Y.canopy.biomass.precip_annual,
            p.canopy.biomass.f0,
            p.drivers.c_co2 * p.drivers.P,  # ca_pa: CO2 partial pressure (Pa)
            p.canopy.biomass.OptVars.χ,
            p.canopy.biomass.vpd_gs,
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
        @. dY.canopy.biomass.precip_annual = apply_time_reduction(
            -(p.drivers.P_liq + p.drivers.P_snow) * ρ_m_liq,
            Y.canopy.biomass.precip_annual,
            tivs.precip_annual.reduction,
        )
        # PET_annual / precip_annual is the aridity index behind f0.
        @. dY.canopy.biomass.PET_annual = apply_time_reduction(
            potential_evaporation(
                p.drivers.SW_d,
                p.drivers.LW_d,
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                p.drivers.u,
                h_atmos,
                ϵ_sfc,
                σ,
                M_w,
                thermo_params,
            ),
            Y.canopy.biomass.PET_annual,
            tivs.PET_annual.reduction,
        )
        # VPDA0_annual / A0_annual is vpd_gs, the A0-weighted mean VPD.
        @. dY.canopy.biomass.VPDA0_annual = apply_time_reduction(
            VPD * p.canopy.biomass.OptVars.A0,
            Y.canopy.biomass.VPDA0_annual,
            tivs.VPDA0_annual.reduction,
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
        # potential is scaled by the realized fAPAR before the yearly total.
        @. dY.canopy.biomass.GPPc3_annual = apply_time_reduction(
            p.canopy.biomass.OptVars.A0_c3 *
            (1 - exp(-parameters.k * Y.canopy.biomass.LAI)),
            Y.canopy.biomass.GPPc3_annual,
            tivs.GPPc3_annual.reduction,
        )
        @. dY.canopy.biomass.LAI = apply_time_reduction(
            p.canopy.biomass.L_opt,
            Y.canopy.biomass.LAI,
            tivs.LAI.reduction,
        )
    end
    return compute_exp_tendency!
end

#####################################################################
# PrognosticCarbonModel - live carbon pools wrapping an LAI model
#####################################################################

"""
    PrognosticCarbonParameters{FT <: AbstractFloat}

Parameters of the live carbon pools of [`PrognosticCarbonModel`](@ref). The leaf and
stem allocation fractions and the stem turnover time are blended between C3 and C4
values by the C3 fraction of the canopy; roots take the rest of the allocation.

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrognosticCarbonParameters{FT <: AbstractFloat}
    "Construction efficiency (-): the fraction of the sugar allocated to growth that becomes structure; the rest is growth respiration"
    a::FT
    "Leaf allocation fraction of C3 vegetation (-)"
    f_leaf_c3::FT
    "Stem allocation fraction of C3 vegetation (-)"
    f_stem_c3::FT
    "Leaf allocation fraction of C4 vegetation (-)"
    f_leaf_c4::FT
    "Stem allocation fraction of C4 vegetation (-)"
    f_stem_c4::FT
    "Leaf turnover time (s)"
    τ_leaf::FT
    "Stem turnover time of C3 vegetation at or above `T_ref_τ_stem` (s)"
    τ_stem_c3::FT
    "Stem turnover time of C4 vegetation at or above `T_ref_τ_stem` (s)"
    τ_stem_c4::FT
    "Fine-root turnover time (s)"
    τ_root::FT
    "Sapwood maintenance respiration rate at `T_ref` (s^-1)"
    r_stem::FT
    "Fine-root maintenance respiration rate at `T_ref` (s^-1)"
    r_root::FT
    "Sapwood carbon as the stem grows large (kg C m^-2); sapwood is half the stem when `C_stem` equals it"
    C_sap_half::FT
    "Target sugar pool as a fraction of the living biomass (-)"
    c_nsc::FT
    "Timescale of allocation from the sugar pool (s)"
    τ_alloc::FT
    "Exponent of the allocation ramp (-)"
    n_alloc::FT
    "Sugar pool below which maintenance respiration shuts down (kg C m^-2)"
    C_sugar_ref::FT
    "Q10 of sapwood and fine-root maintenance respiration (-)"
    Q10::FT
    "Reference temperature of the maintenance respiration rates (K)"
    T_ref::FT
    "Factor by which stem turnover time lengthens per 10 K of mean annual temperature below `T_ref_τ_stem` (-); 1 disables"
    q_τ_stem::FT
    "Mean annual temperature below which stem turnover time lengthens (K)"
    T_ref_τ_stem::FT
    "Mean annual precipitation at which stem allocation is halved (m yr^-1); 0 disables"
    map_half_woody::FT
    "Exponent of the precipitation limit on stem allocation (-)"
    n_map_woody::FT
    "Memory timescale of the mean annual temperature and precipitation (s)"
    τ_climate::FT
    "E-folding depth of the leaf and stem litter input to soil carbon (m)"
    soil_litter_depth::FT
    "Molar mass of carbon (kg mol^-1)"
    M_C::FT
end

Base.eltype(::PrognosticCarbonParameters{FT}) where {FT} = FT

"""
    PrognosticCarbonParameters(toml_dict::CP.ParamDict; kwargs...)

Constructs `PrognosticCarbonParameters` from a TOML dictionary; any parameter can be
overridden by keyword argument.
"""
function PrognosticCarbonParameters(
    toml_dict::CP.ParamDict;
    a = toml_dict["carbon_construction_efficiency"],
    f_leaf_c3 = toml_dict["carbon_f_leaf_c3"],
    f_stem_c3 = toml_dict["carbon_f_stem_c3"],
    f_leaf_c4 = toml_dict["carbon_f_leaf_c4"],
    f_stem_c4 = toml_dict["carbon_f_stem_c4"],
    τ_leaf = toml_dict["carbon_tau_leaf"],
    τ_stem_c3 = toml_dict["carbon_tau_stem_c3"],
    τ_stem_c4 = toml_dict["carbon_tau_stem_c4"],
    τ_root = toml_dict["carbon_tau_root"],
    r_stem = toml_dict["carbon_r_stem"],
    r_root = toml_dict["carbon_r_root"],
    C_sap_half = toml_dict["carbon_C_sap_half"],
    c_nsc = toml_dict["carbon_c_nsc"],
    τ_alloc = toml_dict["carbon_tau_alloc"],
    n_alloc = toml_dict["carbon_alloc_ramp_n"],
    C_sugar_ref = toml_dict["carbon_C_sugar_ref"],
    Q10 = toml_dict["carbon_Q10"],
    T_ref = toml_dict["carbon_T_ref"],
    q_τ_stem = toml_dict["carbon_tau_stem_q"],
    T_ref_τ_stem = toml_dict["carbon_tau_stem_T_ref"],
    map_half_woody = toml_dict["carbon_map_half_woody"],
    n_map_woody = toml_dict["carbon_n_map_woody"],
    τ_climate = toml_dict["carbon_tau_climate"],
    soil_litter_depth = toml_dict["carbon_soil_litter_depth"],
    M_C = toml_dict["molar_mass_carbon"],
)
    FT = CP.float_type(toml_dict)
    return PrognosticCarbonParameters{FT}(;
        a,
        f_leaf_c3,
        f_stem_c3,
        f_leaf_c4,
        f_stem_c4,
        τ_leaf,
        τ_stem_c3,
        τ_stem_c4,
        τ_root,
        r_stem,
        r_root,
        C_sap_half,
        c_nsc,
        τ_alloc,
        n_alloc,
        C_sugar_ref,
        Q10,
        T_ref,
        q_τ_stem,
        T_ref_τ_stem,
        map_half_woody,
        n_map_woody,
        τ_climate,
        soil_litter_depth,
        M_C,
    )
end

"""
    PrognosticCarbonModel{FT, LM, PCP, RDTH, HTH, TIV} <: AbstractBiomassModel{FT}

Live vegetation carbon in four prognostic pools (kg C m^-2): `C_sugar` (non-structural
carbon), `C_leaf`, `C_stem` and `C_root`, driven by the GPP of the canopy:

    dC_sugar/dt = GPP - Rm - S
    dC_leaf/dt  = a f_leaf S - C_leaf/τ_leaf
    dC_stem/dt  = a f_stem S - C_stem/τ_stem
    dC_root/dt  = a f_root S - C_root/τ_root

`S` is the sugar allocated to growth, `(1 - a) S` the growth respiration, `Rm` the
maintenance respiration, and the turnover terms are the litter passed to soil carbon;
see [`update_carbon_fluxes!`](@ref). Stem allocation decreases in dry climates
([`woody_fraction`](@ref)) and stem turnover time increases in cold ones
([`tau_stem_scale`](@ref)), from the mean annual precipitation and temperature carried
as time-integrated variables (`P_annual`, `T_annual`).

The pools wrap an LAI model, `lai_model` (`PrescribedBiomassModel` or
`ZhouOptimalLAIModel`), which still sets the area indices: GPP and LAI are the same
as without the pools. The canopy respiration comes from the pools, so this model
requires [`PoolBasedAutotrophicRespirationModel`](@ref).

$(DocStringExtensions.FIELDS)
"""
struct PrognosticCarbonModel{
    FT,
    LM <: AbstractBiomassModel{FT},
    PCP <: PrognosticCarbonParameters{FT},
    RDTH,
    HTH,
    TIV,
} <: AbstractBiomassModel{FT}
    "The LAI model setting the area indices"
    lai_model::LM
    "Parameters of the carbon pools"
    parameters::PCP
    "Rooting depth (m), that of `lai_model`"
    rooting_depth::RDTH
    "Canopy height (m), that of `lai_model`"
    height::HTH
    "Time-integrated variables: mean annual temperature and precipitation"
    time_integrated_vars::TIV
end

Base.eltype(::PrognosticCarbonModel{FT}) where {FT} = FT

"""
    PrognosticCarbonModel{FT}(lai_model, parameters::PrognosticCarbonParameters{FT})
    PrognosticCarbonModel{FT}(lai_model, toml_dict::CP.ParamDict; kwargs...)

Wraps `lai_model` in the carbon pools, with parameters given directly or read from
`toml_dict` (keyword arguments override them).
"""
function PrognosticCarbonModel{FT}(
    lai_model::AbstractBiomassModel{FT},
    parameters::PrognosticCarbonParameters{FT},
) where {FT}
    year = FT(365 * 86400)
    tivs = ClimaLand.time_integrated_variables(
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :T_annual,
            reduction = ClimaLand.RunningMean(parameters.τ_climate),
        ),
        # A yearly total, in m yr^-1.
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :P_annual,
            reduction = ClimaLand.RunningSum(year, parameters.τ_climate),
        ),
    )
    args =
        (lai_model, parameters, lai_model.rooting_depth, lai_model.height, tivs)
    return PrognosticCarbonModel{FT, typeof.(args)...}(args...)
end

function PrognosticCarbonModel{FT}(
    lai_model::AbstractBiomassModel{FT},
    toml_dict::CP.ParamDict;
    kwargs...,
) where {FT}
    parameters = PrognosticCarbonParameters(toml_dict; kwargs...)
    return PrognosticCarbonModel{FT}(lai_model, parameters)
end

ClimaLand.prognostic_vars(m::PrognosticCarbonModel) = (
    ClimaLand.prognostic_vars(m.lai_model)...,
    :C_sugar,
    :C_leaf,
    :C_stem,
    :C_root,
    ClimaLand.time_integrated_prognostic_vars(m.time_integrated_vars)...,
)
ClimaLand.prognostic_types(m::PrognosticCarbonModel{FT}) where {FT} = (
    ClimaLand.prognostic_types(m.lai_model)...,
    FT,
    FT,
    FT,
    FT,
    ClimaLand.time_integrated_prognostic_types(m.time_integrated_vars)...,
)
ClimaLand.prognostic_domain_names(m::PrognosticCarbonModel) = (
    ClimaLand.prognostic_domain_names(m.lai_model)...,
    :surface,
    :surface,
    :surface,
    :surface,
    ClimaLand.time_integrated_prognostic_domain_names(
        m.time_integrated_vars,
    )...,
)

"""
    ClimaLand.auxiliary_vars(model::PrognosticCarbonModel)

Adds to the auxiliary variables of the LAI model:
- `carbon`: the carbon fluxes (kg C m^-2 s^-1) of [`update_carbon_fluxes!`](@ref):
  maintenance, growth and total respiration `Rm`, `Rg`, `Ra`, the allocation `S`,
  and the litter `L_leaf`, `L_stem`, `L_root`
- `cVeg`: the live carbon (kg C m^-2), the sum of the pools
- `σl_implied`: the carbon per unit leaf area implied by the leaf pool and the LAI
  model, `C_leaf/LAI` (kg C m^-2 leaf), zero where there is no leaf area
"""
ClimaLand.auxiliary_vars(model::PrognosticCarbonModel) =
    (ClimaLand.auxiliary_vars(model.lai_model)..., :carbon, :cVeg, :σl_implied)
ClimaLand.auxiliary_types(model::PrognosticCarbonModel{FT}) where {FT} = (
    ClimaLand.auxiliary_types(model.lai_model)...,
    NamedTuple{(:Rm, :Rg, :Ra, :S, :L_leaf, :L_stem, :L_root), NTuple{7, FT}},
    FT,
    FT,
)
ClimaLand.auxiliary_domain_names(model::PrognosticCarbonModel) = (
    ClimaLand.auxiliary_domain_names(model.lai_model)...,
    :surface,
    :surface,
    :surface,
)

prescribed_lai_input(model::PrognosticCarbonModel) =
    prescribed_lai_input(model.lai_model)

get_fractional_c3(p, biomass::PrognosticCarbonModel, photosynthesis) =
    get_fractional_c3(p, biomass.lai_model, photosynthesis)

"""
    sapwood_carbon(C_stem, C_sap_half)

Living (sapwood) carbon of a stem pool `C_stem`, `C_stem/(1 + C_stem/C_sap_half)`:
linear in a small stem and saturating at `C_sap_half` in a large one, whose wood is
mostly dead heartwood.
"""
function sapwood_carbon(C_stem::FT, C_sap_half::FT) where {FT}
    return C_stem / (1 + C_stem / C_sap_half)
end

"""
    allocation_ramp(x, n)

The smooth ramp `x^n/(1 + x^n)` (zero for `x ≤ 0`), which switches a flux on as the
sugar pool rises past a reference, `x` being the ratio of the two.
"""
function allocation_ramp(x::FT, n::FT) where {FT}
    xn = max(x, zero(FT))^n
    return xn / (1 + xn)
end

# Cap on `tau_stem_scale`: 300 years with the default 30-year stem turnover.
const MAX_TAU_STEM_SCALE = 10

"""
    tau_stem_scale(MAT, T_ref, q)

Factor on the stem turnover time for a mean annual temperature `MAT`:
`q^((T_ref - MAT)/10)` below `T_ref` and 1 above it, capped at `MAX_TAU_STEM_SCALE`.
Cold-climate trees live longer; `q ≤ 1` disables the scaling.
"""
function tau_stem_scale(MAT::FT, T_ref::FT, q::FT) where {FT}
    q <= 1 && return one(FT)
    return min(q^(max(T_ref - MAT, zero(FT)) / 10), FT(MAX_TAU_STEM_SCALE))
end

"""
    woody_fraction(MAP, half, n)

Factor on the stem allocation fraction for a mean annual precipitation `MAP`,
`x^n/(1 + x^n)` with `x = MAP/half`: dry climates build little wood (Sankaran et al.,
2005). The allocation it withholds goes to roots. `half ≤ 0` disables it.
"""
function woody_fraction(MAP::FT, half::FT, n::FT) where {FT}
    half <= 0 && return one(FT)
    return allocation_ramp(MAP / half, n)
end

"""
    update_biomass!(p, Y, t, component::PrognosticCarbonModel, canopy)

Updates the area indices with the LAI model, then the live carbon `cVeg` and
`σl_implied` from the pools.
"""
function update_biomass!(
    p,
    Y,
    t,
    component::PrognosticCarbonModel{FT},
    canopy,
) where {FT}
    update_biomass!(p, Y, t, component.lai_model, canopy)
    (; C_sugar, C_leaf, C_stem, C_root) = Y.canopy.biomass
    @. p.canopy.biomass.cVeg = C_sugar + C_leaf + C_stem + C_root
    LAI = p.canopy.biomass.area_index.leaf
    @. p.canopy.biomass.σl_implied =
        ifelse(LAI > 0, C_leaf / max(LAI, eps(FT)), zero(FT))
    return nothing
end

"""
    update_carbon_fluxes!(p, Y, biomass, canopy)

Updates the carbon fluxes of the pools, `p.canopy.biomass.carbon` (kg C m^-2 s^-1);
does nothing for a biomass model without pools. Called in `update_aux` after
photosynthesis, which supplies GPP and the leaf respiration `Rd`.

Maintenance respiration is that of the leaves, `Rd`, plus that of the sapwood and fine
roots, at rates `r_stem` and `r_root` scaled by a Q10 of the canopy temperature:

    Rm = g(C_sugar/C_sugar_ref) [Rd + Q10^((T - T_ref)/10) (r_stem C_sap + r_root C_root)]

with `C_sap` the sapwood carbon (`sapwood_carbon`) and `g` the ramp `allocation_ramp`,
which stops respiration as the sugar pool empties. Allocation draws the sugar pool toward a
target `c_nsc (C_leaf + C_sap + C_root)`:

    S = C_sugar/τ_alloc g(C_sugar/(c_nsc (C_leaf + C_sap + C_root)))

The litter is the turnover of each structural pool, with the stem turnover time scaled
by [`tau_stem_scale`](@ref).
"""
update_carbon_fluxes!(p, Y, biomass::AbstractBiomassModel, canopy) = nothing

function update_carbon_fluxes!(
    p,
    Y,
    biomass::PrognosticCarbonModel{FT},
    canopy,
) where {FT}
    (;
        a,
        τ_leaf,
        τ_stem_c3,
        τ_stem_c4,
        τ_root,
        r_stem,
        r_root,
        C_sap_half,
        c_nsc,
        τ_alloc,
        n_alloc,
        C_sugar_ref,
        Q10,
        T_ref,
        q_τ_stem,
        T_ref_τ_stem,
        M_C,
    ) = biomass.parameters
    (; C_sugar, C_leaf, C_stem, C_root, T_annual) = Y.canopy.biomass
    carbon = p.canopy.biomass.carbon
    # Accessors resolved outside `@.`, which would otherwise broadcast over `p`.
    fractional_c3 = get_fractional_c3(p, canopy)
    Rd = get_Rd_canopy(p, canopy.photosynthesis)
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    C_sap = @. lazy(sapwood_carbon(C_stem, C_sap_half))

    # Rd has its own temperature response. The roots use the canopy temperature,
    # as soil temperature is not available to a standalone canopy.
    @. carbon.Rm =
        allocation_ramp(C_sugar / C_sugar_ref, n_alloc) * (
            M_C * Rd +
            Q10^((T_canopy - T_ref) / 10) * (r_stem * C_sap + r_root * C_root)
        )
    @. carbon.S =
        max(C_sugar, 0) / τ_alloc * allocation_ramp(
            C_sugar / max(c_nsc * (C_leaf + C_sap + C_root), eps(FT)),
            n_alloc,
        )
    @. carbon.Rg = (1 - a) * carbon.S
    @. carbon.Ra = carbon.Rm + carbon.Rg
    @. carbon.L_leaf = C_leaf / τ_leaf
    @. carbon.L_stem =
        C_stem / (
            blend(τ_stem_c3, τ_stem_c4, fractional_c3) *
            tau_stem_scale(T_annual, T_ref_τ_stem, q_τ_stem)
        )
    @. carbon.L_root = C_root / τ_root
    return nothing
end

"""
    equilibrium_carbon_pools(parameters, GPP, Rd, f_T, MAT, MAP, fractional_c3)

Steady state of the structural carbon pools of [`PrognosticCarbonModel`](@ref), as
`(; C_leaf, C_stem, C_root)` (kg C m^-2), under constant forcing: GPP and leaf
respiration `Rd` (mol CO2 m^-2 s^-1), the Q10 factor `f_T = Q10^((T - T_ref)/10)` of
sapwood and fine-root respiration, the mean annual temperature `MAT` (K) and
precipitation `MAP` (m yr^-1), and the C3 fraction. For periodic forcing, pass the means
over the period.

In steady state, each structural pool is its allocation times its turnover time,
`C_i = a f_i τ_i S`, and the allocation `S` is the GPP left after maintenance
respiration, `S = M_C (GPP - Rd) - f_T (r_stem C_sap + r_root C_root)`. As the sapwood
`C_sap` saturates with `C_stem = a f_stem τ_stem S`, this is a quadratic in `S`. The
small sugar pool, and the limitation of respiration by an empty sugar pool, are
neglected.
"""
function equilibrium_carbon_pools(
    parameters::PrognosticCarbonParameters{FT},
    GPP,
    Rd,
    f_T,
    MAT,
    MAP,
    fractional_c3,
) where {FT}
    (; a, f_leaf_c3, f_leaf_c4, f_stem_c3, f_stem_c4) = parameters
    (; τ_leaf, τ_stem_c3, τ_stem_c4, τ_root, r_stem, r_root) = parameters
    (; C_sap_half, q_τ_stem, T_ref_τ_stem, map_half_woody, n_map_woody, M_C) =
        parameters
    f_leaf = blend(f_leaf_c3, f_leaf_c4, fractional_c3)
    f_stem =
        blend(f_stem_c3, f_stem_c4, fractional_c3) *
        woody_fraction(MAP, map_half_woody, n_map_woody)
    f_root = 1 - f_leaf - f_stem
    τ_stem =
        blend(τ_stem_c3, τ_stem_c4, fractional_c3) *
        tau_stem_scale(MAT, T_ref_τ_stem, q_τ_stem)
    A = max(M_C * (GPP - Rd), zero(FT))
    k = a * f_stem * τ_stem # stem carbon per unit allocation
    ρ = f_T * r_root * a * f_root * τ_root # root respiration per unit allocation
    # (1 + ρ) S + f_T r_stem k S/(1 + k S/C_sap_half) = A, i.e. α S^2 + β S - A = 0
    α = (1 + ρ) * k / C_sap_half
    β = 1 + ρ + f_T * r_stem * k - A * k / C_sap_half
    S = 2 * A / (β + sqrt(β^2 + 4 * α * A))
    return (;
        C_leaf = a * f_leaf * τ_leaf * S,
        C_stem = k * S,
        C_root = a * f_root * τ_root * S,
    )
end

"""
    ClimaLand.make_compute_exp_tendency(component::PrognosticCarbonModel, canopy)

Advances the prognostic variables of the LAI model, then the carbon pools (see
[`PrognosticCarbonModel`](@ref)) and the mean annual temperature and precipitation.
The fluxes are those computed by [`update_carbon_fluxes!`](@ref).
"""
function ClimaLand.make_compute_exp_tendency(
    component::PrognosticCarbonModel{FT},
    canopy,
) where {FT}
    lai_tendency! =
        ClimaLand.make_compute_exp_tendency(component.lai_model, canopy)
    (; a, f_leaf_c3, f_stem_c3, f_leaf_c4, f_stem_c4) = component.parameters
    (; map_half_woody, n_map_woody, M_C) = component.parameters
    tivs = component.time_integrated_vars
    function compute_exp_tendency!(dY, Y, p, t)
        lai_tendency!(dY, Y, p, t)
        (; S, Rm, L_leaf, L_stem, L_root) = p.canopy.biomass.carbon
        fractional_c3 = get_fractional_c3(p, canopy)
        GPP = get_GPP(p, canopy.photosynthesis)
        f_leaf = @. lazy(blend(f_leaf_c3, f_leaf_c4, fractional_c3))
        f_stem = @. lazy(
            blend(f_stem_c3, f_stem_c4, fractional_c3) * woody_fraction(
                Y.canopy.biomass.P_annual,
                map_half_woody,
                n_map_woody,
            ),
        )
        @. dY.canopy.biomass.C_sugar = M_C * GPP - Rm - S
        @. dY.canopy.biomass.C_leaf = a * f_leaf * S - L_leaf
        @. dY.canopy.biomass.C_stem = a * f_stem * S - L_stem
        @. dY.canopy.biomass.C_root = a * (1 - f_leaf - f_stem) * S - L_root
        @. dY.canopy.biomass.T_annual = apply_time_reduction(
            p.drivers.T,
            Y.canopy.biomass.T_annual,
            tivs.T_annual.reduction,
        )
        # P_liq and P_snow are volume fluxes (m s^-1), negative downward.
        @. dY.canopy.biomass.P_annual = apply_time_reduction(
            -(p.drivers.P_liq + p.drivers.P_snow),
            Y.canopy.biomass.P_annual,
            tivs.P_annual.reduction,
        )
    end
    return compute_exp_tendency!
end

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
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
                                         time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
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
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
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

The prescribed LAI `TimeVaryingInput` a biomass model reads, for the consistency
check in the prescribed-LAI `CanopyModel` constructor. Defined for models that
have one, and forwarded by any model that wraps such a model, so the check does
not need to know how deeply the LAI model is nested.
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

Five further 1-year `RunningSum`s carry the climate the LAI formulas respond to,
replacing what were static maps: `PET_annual` (over `precip_annual`, the aridity
index behind `f0`), `VPDA0_annual` (over `A0_annual`, the A0-weighted growing-season
VPD `vpd_gs`), `growing_days` (the growing-season length `GSL`), and the per-pathway
pair `A0c3_annual`/`A0c4_annual` driving the C3/C4 competition.

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
- `GSL`: growing season length (days), the trailing-year count of growing days
- `vpd_gs`: mean VPD during growing season (Pa), for water limitation WUE factor in LAI_max
- `f0`: fraction of precipitation available for transpiration (dimensionless), from the aridity index

`GSL`, `vpd_gs` and `f0` are recomputed each step in `update_biomass!` from the
trailing climate totals in `Y`.
"""
ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel) =
    (:area_index, :OptVars, :L_opt, :GSL, :vpd_gs, :f0)
ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel{FT}) where {FT} = (
    NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    NamedTuple{(:A0, :A0_c3, :A0_c4, :χ), NTuple{4, FT}},
    FT,
    FT,
    FT,
    FT,
)
ClimaLand.auxiliary_domain_names(::ZhouOptimalLAIModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)

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

Updates the optimal-LAI cache: sets SAI and RAI to their prescribed constant
values, derives the growing-season inputs of the LAI formulas (`f0`, `vpd_gs`,
`GSL`) from the trailing climate totals in `Y`, and mirrors the prognostic `LAI`
into the cache area index read by the rest of the canopy. LAI must be updated here
first in `update_aux`, before radiative transfer reads the area index.
"""
function update_biomass!(
    p,
    Y,
    t,
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    (; SAI, RAI) = component
    @. p.canopy.biomass.area_index.stem = SAI
    @. p.canopy.biomass.area_index.root = RAI
    # The water-limitation inputs track the simulated climate through the trailing
    # totals in `Y`, so they are recomputed here in update_aux, before the tendency
    # reads them for `L_opt`.
    @. p.canopy.biomass.f0 = f0_from_aridity(
        Y.canopy.biomass.PET_annual,
        Y.canopy.biomass.precip_annual,
    )
    @. p.canopy.biomass.vpd_gs =
        Y.canopy.biomass.VPDA0_annual / max(Y.canopy.biomass.A0_annual, eps(FT))
    @. p.canopy.biomass.GSL = Y.canopy.biomass.growing_days
    @. p.canopy.biomass.area_index.leaf = Y.canopy.biomass.LAI
    # Apply clipping to LAI (same as PrescribedBiomassModel)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
end

"""
    update_fractional_c3!(p, Y, biomass::ZhouOptimalLAIModel, canopy)

Resolves the C3 fraction from the C3/C4 competition on the trailing per-pathway
potential GPP (`A0c3_annual`, `A0c4_annual`), rather than taking the P-model's
static value. Setting `optimal_lai_online_c3c4 = 0` recovers that static value.
"""
function update_fractional_c3!(
    p,
    Y,
    biomass::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    online_c3c4 = biomass.parameters.online_c3c4
    k = biomass.parameters.k
    Mc = canopy.photosynthesis.constants.Mc
    static_c3 = canopy.photosynthesis.fractional_c3
    @. p.canopy.photosynthesis.fractional_c3 =
        online_c3c4 * c3_fraction_from_competition(
            Y.canopy.biomass.A0c3_annual,
            Y.canopy.biomass.A0c4_annual,
            Mc,
            1 - exp(-k * Y.canopy.biomass.LAI),  # realized fAPAR
        ) + (1 - online_c3c4) * static_c3
    return nothing
end

"""
    ClimaLand.make_compute_exp_tendency(component::ZhouOptimalLAIModel, canopy)

Advances the optimal-LAI model's nine time-integrated variables.
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
    λv = LP.LH_v0(earth_param_set)              # J kg^-1
    M_w = LP.molar_mass_water(earth_param_set)  # kg mol^-1
    parameters = component.parameters
    pmodel_parameters = canopy.photosynthesis.parameters
    pmodel_constants = canopy.photosynthesis.constants
    function compute_exp_tendency!(dY, Y, p, t)
        # The C3 fraction resolved by the competition in `update_fractional_c3!`, so
        # the potential GPP, chi and the C3/C4 parameter blends all track it.
        fractional_c3 = p.canopy.photosynthesis.fractional_c3
        VPD = @. lazy(
            max(
                Thermodynamics.vapor_pressure_deficit(
                    LP.thermodynamic_parameters(earth_param_set),
                    p.drivers.T,
                    p.drivers.P,
                    p.drivers.q,
                ),
                sqrt(eps(FT)),
            ),
        )

        # A0 is a *potential* GPP, so it carries no soil-moisture stress (βm = 1);
        # water limitation enters once, through the f0·P/A0 term of LAI_max.
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
            a0_mapped_z(
                pft_blend(fractional_c3, parameters.z, parameters.z_c4),
                Y.canopy.biomass.A0_annual,
                parameters.z_a0,
            ),
            p.canopy.biomass.GSL,
            pft_blend(fractional_c3, parameters.sigma, parameters.sigma_c4),
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
        # Over `precip_annual` this is the aridity index behind `f0`.
        @. dY.canopy.biomass.PET_annual = apply_time_reduction(
            potential_evaporation(
                p.drivers.SW_d,
                p.drivers.LW_d,
                p.drivers.T,
                σ,
                λv,
                M_w,
            ),
            Y.canopy.biomass.PET_annual,
            tivs.PET_annual.reduction,
        )
        # Over `A0_annual` this is `vpd_gs`, a VPD mean weighted toward the
        # productive period rather than diluted by the dormant season.
        @. dY.canopy.biomass.VPDA0_annual = apply_time_reduction(
            VPD * p.canopy.biomass.OptVars.A0,
            Y.canopy.biomass.VPDA0_annual,
            tivs.VPDA0_annual.reduction,
        )
        # One growing-day per day with air T > 0 C, so the trailing-year total is the
        # growing-season length in days.
        @. dY.canopy.biomass.growing_days = apply_time_reduction(
            ifelse(p.drivers.T > FT(273.15), 1 / seconds_per_day, zero(FT)),
            Y.canopy.biomass.growing_days,
            tivs.growing_days.reduction,
        )
        # The pure-C3/C4 pair drives the competition through (A0c4 - A0c3)/A0c3.
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

Parameters for the prognostic live-carbon pools.

Allocation fractions and turnover times are global constants, optionally split
C3/C4 and blended by the modelled C3 fraction. They never vary with a plant
functional type. `f_root` is not stored: it is `1 - f_leaf - f_stem`, so the
fractions sum to one by construction.

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrognosticCarbonParameters{FT <: AbstractFloat}
    "Construction efficiency a (-); (1-a) of the sugar drawn for structure is growth respiration"
    a::FT
    "Allocation fraction to leaves, C3 (-)"
    f_leaf_c3::FT
    "Allocation fraction to stem, C3 (-)"
    f_stem_c3::FT
    "Allocation fraction to leaves, C4 (-)"
    f_leaf_c4::FT
    "Allocation fraction to stem, C4 (-)"
    f_stem_c4::FT
    "Leaf lifespan (s)"
    τ_leaf::FT
    "Stem turnover time, C3 (s)"
    τ_stem_c3::FT
    "Stem turnover time, C4 (s)"
    τ_stem_c4::FT
    "Fine root turnover time (s)"
    τ_root::FT
    "Sapwood specific maintenance respiration rate (s^-1)"
    r_stem::FT
    "Root specific maintenance respiration rate (s^-1)"
    r_root::FT
    "Sapwood saturation scale (kg C m^-2)"
    C_sap_half::FT
    "Target sugar store as a fraction of live biomass (-)"
    c_nsc::FT
    "Allocation timescale (s)"
    τ_alloc::FT
    "Sharpness of the allocation ramp (-)"
    n_alloc::FT
    "Q10 temperature sensitivity of maintenance respiration (-)"
    Q10::FT
    "Reference temperature for the Q10 factor (K)"
    T_ref::FT
end

Base.eltype(::PrognosticCarbonParameters{FT}) where {FT} = FT

"""
    PrognosticCarbonParameters(toml_dict::CP.ParamDict; kwargs...)

Constructor for `PrognosticCarbonParameters` from a TOML dictionary. Any
parameter can be overridden by keyword argument.
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
    Q10 = toml_dict["carbon_Q10"],
    T_ref = toml_dict["carbon_T_ref"],
)
    FT = CP.float_type(toml_dict)
    PrognosticCarbonParameters{FT}(;
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
        Q10,
        T_ref,
    )
end

"""
    PrognosticCarbonModel{FT, LM <: AbstractBiomassModel{FT}, PCP, RDTH, HTH} <: AbstractBiomassModel{FT}

Four prognostic live-carbon pools - `C_sugar`, `C_leaf`, `C_stem`, `C_root`, all
in kg C m^-2 of ground - wrapped around an existing LAI model.

The carbon model does not compute LAI. It *composes* with a `lai_model`, either
`PrescribedBiomassModel` or `ZhouOptimalLAIModel`, and delegates every area-index
and LAI decision to it. It consumes GPP, tracks the pools, and returns
respiration and litter. This is what makes the phase-1 coupling one-way and
verifiable: a run with the carbon model reproduces the same GPP and LAI as one
without it.

In this stage the pools are driven by the carbon model's own internal
maintenance and growth respiration, but that respiration is *not* yet wired into
`p.canopy.autotrophic_respiration.Ra` - the existing scheme still drives the
canopy fluxes, so behaviour outside the pools is unchanged. Wiring it in is a
separate, explicit step.

$(DocStringExtensions.FIELDS)
"""
struct PrognosticCarbonModel{
    FT,
    LM <: AbstractBiomassModel{FT},
    PCP <: PrognosticCarbonParameters{FT},
    RDTH,
    HTH,
} <: AbstractBiomassModel{FT}
    "The underlying LAI/area-index model the carbon pools compose with"
    lai_model::LM
    "Parameters of the carbon pools"
    parameters::PCP
    "Rooting depth (m); mirrored from `lai_model` so the rest of the canopy can read it directly"
    rooting_depth::RDTH
    "Canopy height (m); mirrored from `lai_model`"
    height::HTH
end

Base.eltype(::PrognosticCarbonModel{FT}) where {FT} = FT

"""
    PrognosticCarbonModel{FT}(lai_model, parameters)

Outer constructor wrapping `lai_model` in the prognostic carbon pools.

`rooting_depth` and `height` are copied from `lai_model` because much of the
canopy reads `canopy.biomass.rooting_depth` and `canopy.biomass.height`
directly. In phase 1 both stay prescribed and neither is ever mutated, so the
copies cannot drift from their source.
"""
function PrognosticCarbonModel{FT}(
    lai_model::AbstractBiomassModel{FT},
    parameters::PrognosticCarbonParameters{FT},
) where {FT}
    rooting_depth = lai_model.rooting_depth
    height = lai_model.height
    args = (lai_model, parameters, rooting_depth, height)
    PrognosticCarbonModel{FT, typeof.(args)...}(args...)
end

"""
    PrognosticCarbonModel{FT}(lai_model, toml_dict::CP.ParamDict; kwargs...)

Convenience constructor building the parameters from a TOML dictionary.
"""
function PrognosticCarbonModel{FT}(
    lai_model::AbstractBiomassModel{FT},
    toml_dict::CP.ParamDict;
    kwargs...,
) where {FT}
    parameters = PrognosticCarbonParameters(toml_dict; kwargs...)
    PrognosticCarbonModel{FT}(lai_model, parameters)
end

# The pools are appended to whatever the wrapped LAI model already carries, so
# a prescribed-LAI wrap adds four prognostic variables and a Zhou wrap adds the
# same four alongside its nine time-integrated ones.
const CARBON_POOLS = (:C_sugar, :C_leaf, :C_stem, :C_root)

ClimaLand.prognostic_vars(m::PrognosticCarbonModel) =
    (ClimaLand.prognostic_vars(m.lai_model)..., CARBON_POOLS...)
ClimaLand.prognostic_types(m::PrognosticCarbonModel{FT}) where {FT} =
    (ClimaLand.prognostic_types(m.lai_model)..., FT, FT, FT, FT)
ClimaLand.prognostic_domain_names(m::PrognosticCarbonModel) = (
    ClimaLand.prognostic_domain_names(m.lai_model)...,
    :surface,
    :surface,
    :surface,
    :surface,
)

"""
    ClimaLand.auxiliary_vars(model::PrognosticCarbonModel)

Adds to the wrapped LAI model's cache:
- `carbon`: the carbon fluxes (kg C m^-2 s^-1) - maintenance respiration `Rm`,
  growth respiration `Rg`, total autotrophic respiration `Ra` of the carbon
  model, the structural allocation rate `S`, and the leaf, stem and root litter
  fluxes `L_leaf`, `L_stem`, `L_root`
- `cVeg`: total live carbon (kg C m^-2), the sum of the four pools
- `σl_implied`: `C_leaf/LAI` (kg C m^-2 leaf), the diagnostic that reveals
  whether the constant allocation fractions and the LAI model agree
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

"""
    sapwood_carbon(C_stem::FT, C_sap_half::FT) where {FT}

The living (sapwood) fraction of the stem pool,
`C_sap = C_stem/(1 + C_stem/C_sap_half)`, which saturates at `C_sap_half` as
`C_stem` grows. Living wood does not scale linearly with total wood, and this is
what keeps a heavy tropical stem pool from respiring itself to death.
"""
function sapwood_carbon(C_stem::FT, C_sap_half::FT) where {FT}
    return C_stem / (1 + C_stem / C_sap_half)
end

"""
    allocation_ramp(x::FT, n::FT) where {FT}

The smooth ramp `g(x) = x^n/(1 + x^n)` regulating the allocation rate, where
`x = C_sugar/C_sugar_target`. Sugar well above target gives `g -> 1` and growth
runs at its maximum, drawing the pool down; sugar well below target gives
`g -> 0` and growth stops while maintenance respiration continues. That
proportional control is what makes the sugar pool oscillate seasonally rather
than pin at zero, and it needs no hard clamp.
"""
function allocation_ramp(x::FT, n::FT) where {FT}
    xn = max(x, zero(FT))^n
    return xn / (1 + xn)
end

"""
    pft_free_blend(fractional_c3::FT, v_c3::FT, v_c4::FT) where {FT}

Blends a C3 and a C4 parameter value by the modelled C3 fraction. This is a
photosynthetic-pathway blend computed from GPP, not a plant functional type.
"""
function pft_free_blend(fractional_c3::FT, v_c3::FT, v_c4::FT) where {FT}
    return fractional_c3 * v_c3 + (1 - fractional_c3) * v_c4
end

"""
    update_biomass!(p, Y, t, component::PrognosticCarbonModel, canopy)

Delegates the area indices and LAI entirely to the wrapped LAI model, then
updates the carbon diagnostics. Delegating first is what guarantees LAI is
unchanged by the presence of the carbon model.
"""
function update_biomass!(
    p,
    Y,
    t,
    component::PrognosticCarbonModel{FT},
    canopy,
) where {FT}
    update_biomass!(p, Y, t, component.lai_model, canopy)
    @. p.canopy.biomass.cVeg =
        Y.canopy.biomass.C_sugar +
        Y.canopy.biomass.C_leaf +
        Y.canopy.biomass.C_stem +
        Y.canopy.biomass.C_root
    # Reported against the prescribed specific leaf density; far from ~0.03-0.1
    # kg C m^-2 leaf means the allocation fractions and the LAI model disagree.
    @. p.canopy.biomass.σl_implied =
        Y.canopy.biomass.C_leaf / max(p.canopy.biomass.area_index.leaf, eps(FT))
    return nothing
end

"""
    update_fractional_c3!(p, Y, biomass::PrognosticCarbonModel, canopy)

Forwards to the wrapped LAI model. Without this the Zhou model's C3/C4
competition would silently stop being applied when the carbon model wraps it,
which would change GPP.
"""
function update_fractional_c3!(p, Y, biomass::PrognosticCarbonModel, canopy)
    update_fractional_c3!(p, Y, biomass.lai_model, canopy)
end

"""
    ClimaLand.make_compute_exp_tendency(component::PrognosticCarbonModel, canopy)

Advances the four carbon pools, after advancing whatever prognostic variables
the wrapped LAI model carries.

Following the model specification, with all fluxes in kg C m^-2 s^-1:

    dC_sugar/dt = GPP - Rm - S
    dC_leaf/dt  = a*f_leaf*S - C_leaf/τ_leaf
    dC_stem/dt  = a*f_stem*S - C_stem/τ_stem
    dC_root/dt  = a*f_root*S - C_root/τ_root

where `S = (C_sugar/τ_alloc)*g(C_sugar/C_sugar_target)` is the sugar drawn for
structure, `Rg = (1-a)*S` is growth respiration and `Ra = Rm + Rg`. Summing the
four gives `d(ΣC)/dt = GPP - Ra - litter`, the balance the conservation test
checks.

Maintenance respiration follows the specification,

    Rm = f_T(T_canopy)*(M_C*Rd_canopy + r_stem*C_sap) + f_T(T_soil)*r_root*C_root

reusing the leaf dark respiration the photosynthesis scheme already computes,
which avoids double-counting `Rd` (it is already inside `An = GPP - Rd`). The
root term uses ground rather than canopy temperature.

Phenological leaf shedding (the `C_leaf*max(-dLAI/dt, 0)/LAI` term) is not yet
included; only the background `C_leaf/τ_leaf` turnover is. Deciduous leaf carbon
therefore lags LAI through autumn.
"""
function ClimaLand.make_compute_exp_tendency(
    component::PrognosticCarbonModel{FT},
    canopy,
) where {FT}
    lai_tendency! =
        ClimaLand.make_compute_exp_tendency(component.lai_model, canopy)
    (;
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
        Q10,
        T_ref,
    ) = component.parameters
    M_C = FT(0.012011)  # kg C per mol
    function compute_exp_tendency!(dY, Y, p, t)
        lai_tendency!(dY, Y, p, t)

        fractional_c3 = p.canopy.photosynthesis.fractional_c3
        # Resolved outside the broadcasts below: these select a field from the
        # cache and must not themselves be broadcast over it.
        GPP_mol = get_GPP(p, canopy.photosynthesis)
        Rd_mol = get_Rd_canopy(p, canopy.photosynthesis)
        GPP = @. lazy(M_C * GPP_mol)
        Rd_canopy = @. lazy(M_C * Rd_mol)
        T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
        T_soil = p.drivers.T_ground

        C_sap = @. lazy(sapwood_carbon(Y.canopy.biomass.C_stem, C_sap_half))
        # A floor at zero in the respiration term is the only clamp: it keeps the
        # explicit step from drawing a pool negative. Genuine carbon starvation
        # (sugar at zero and staying there) remains expressible.
        @. p.canopy.biomass.carbon.Rm =
            Q10^((T_canopy - T_ref) / 10) * (Rd_canopy + r_stem * C_sap) +
            Q10^((T_soil - T_ref) / 10) *
            r_root *
            max(Y.canopy.biomass.C_root, zero(FT))

        # The sugar buffer's target is set by the standing live biomass, so
        # allocation throttles smoothly as reserves are drawn down.
        @. p.canopy.biomass.carbon.S =
            max(Y.canopy.biomass.C_sugar, zero(FT)) / τ_alloc * allocation_ramp(
                Y.canopy.biomass.C_sugar / max(
                    c_nsc *
                    (Y.canopy.biomass.C_leaf + C_sap + Y.canopy.biomass.C_root),
                    eps(FT),
                ),
                n_alloc,
            )
        @. p.canopy.biomass.carbon.Rg = (1 - a) * p.canopy.biomass.carbon.S
        @. p.canopy.biomass.carbon.Ra =
            p.canopy.biomass.carbon.Rm + p.canopy.biomass.carbon.Rg

        @. p.canopy.biomass.carbon.L_leaf = Y.canopy.biomass.C_leaf / τ_leaf
        @. p.canopy.biomass.carbon.L_stem =
            Y.canopy.biomass.C_stem /
            pft_free_blend(fractional_c3, τ_stem_c3, τ_stem_c4)
        @. p.canopy.biomass.carbon.L_root = Y.canopy.biomass.C_root / τ_root

        @. dY.canopy.biomass.C_sugar =
            GPP - p.canopy.biomass.carbon.Rm - p.canopy.biomass.carbon.S
        @. dY.canopy.biomass.C_leaf =
            a *
            pft_free_blend(fractional_c3, f_leaf_c3, f_leaf_c4) *
            p.canopy.biomass.carbon.S - p.canopy.biomass.carbon.L_leaf
        @. dY.canopy.biomass.C_stem =
            a *
            pft_free_blend(fractional_c3, f_stem_c3, f_stem_c4) *
            p.canopy.biomass.carbon.S - p.canopy.biomass.carbon.L_stem
        # f_root is the remainder, so the three fractions sum to one exactly.
        @. dY.canopy.biomass.C_root =
            a *
            (
                1 - pft_free_blend(fractional_c3, f_leaf_c3, f_leaf_c4) -
                pft_free_blend(fractional_c3, f_stem_c3, f_stem_c4)
            ) *
            p.canopy.biomass.carbon.S - p.canopy.biomass.carbon.L_root
    end
    return compute_exp_tendency!
end

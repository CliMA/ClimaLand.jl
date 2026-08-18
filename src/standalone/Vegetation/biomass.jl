import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
using Dates
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import Interpolations: Constant
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
using DocStringExtensions

export prescribed_lai_era5,
    TwoComponentResNetLAIModel,
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
- `optimal_lai_inputs`: NamedTuple with the spatially varying inputs of the LAI_max and
  steady-state-LAI formulas: GSL (growing season length in days), vpd_gs (growing-season
  mean VPD in Pa), and f0 (fraction of precip for transpiration). Typically created using
  `optimal_lai_inputs`, which reads them alongside the initial conditions.
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
    GD,
    FS <: Union{FT, ClimaCore.Fields.Field},
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
    T,
} <: AbstractBiomassModel{FT}
    "Required parameters for the optimal LAI model"
    parameters::OLPT
    "Spatially varying inputs of the LAI formulas (GSL, vpd_gs, f0)"
    optimal_lai_inputs::GD
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
        parameters::OptimalLAIParameters{FT},
        optimal_lai_inputs;
        SAI,
        RAI,
        rooting_depth,
        height,
    ) where {FT <: AbstractFloat}

Outer constructor for the ZhouOptimalLAIModel struct.

# Arguments
- `parameters`: OptimalLAIParameters for the model
- `optimal_lai_inputs`: NamedTuple with the spatially varying GSL, vpd_gs and f0 fields,
  typically created using `optimal_lai_inputs`.
- `SAI`: Prescribed stem area index (m^2 m^-2); scalar or spatially-varying Field
- `RAI`: Prescribed root area index (m^2 m^-2); scalar or spatially-varying Field
- `rooting_depth`: Rooting depth parameter (m)
- `height`: Canopy height (m) - can be scalar or spatially-varying Field

Declares the prognostic time integrated variables: the 1-day potential-GPP total
`A0_daily` and the 1-year totals `A0_annual` and `precip_annual` as `RunningSum`s of
the instantaneous rate, and `LAI` as a `RunningMean` relaxing toward the
instantaneous steady-state target `L_opt`:

    dA0_daily/dt      = (day·A0_inst - A0_daily) / τ_day,        τ_day  = 3 days,
    dA0_annual/dt     = (year·A0_inst - A0_annual) / τ_long,     τ_long = tau_long_term,
    dprecip_annual/dt = (year·P_inst  - precip_annual) / τ_long,
    dLAI/dt           = (L_opt - LAI) / τ_LAI,                   τ_LAI  = 1 day / α.

Each `RunningSum` holds a total over its own window (1 day, 1 year) whatever the
smoothing timescale τ_long: only the smoothing changes with τ_long, not the magnitude,
as the LAI_max/steady-state formulas require.
"""
function ZhouOptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT},
    optimal_lai_inputs;
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
            name = :LAI,
            reduction = ClimaLand.RunningMean(
                seconds_per_day / parameters.alpha,
            ),
        ),
    )
    return ZhouOptimalLAIModel{
        FT,
        typeof(parameters),
        typeof(optimal_lai_inputs),
        typeof(SAI),
        typeof(rooting_depth),
        typeof(height),
        typeof(tiv),
    }(
        parameters,
        optimal_lai_inputs,
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
- `L_opt`: Optimal LAI predicted by Zhou et al.
"""
ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel) =
    (:area_index, :OptVars, :L_opt)
ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel{FT}) where {FT} = (
    NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    NamedTuple{(:A0, :χ), Tuple{FT, FT}},
    FT,
)
ClimaLand.auxiliary_domain_names(::ZhouOptimalLAIModel) =
    (:surface, :surface, :surface)

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
values and mirrors the prognostic `LAI` into the cache area index read by the rest
of the canopy. LAI must be updated here first in `update_aux`, before radiative
transfer reads the area index.
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
    @. p.canopy.biomass.area_index.leaf = Y.canopy.biomass.LAI
    # Apply clipping to LAI (same as PrescribedBiomassModel)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
end

"""
    ClimaLand.make_compute_exp_tendency(component::ZhouOptimalLAIModel, canopy)

Advances the optimal-LAI model's four time-integrated variables.
"""
function ClimaLand.make_compute_exp_tendency(
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    ρ_m_liq = LP.ρ_m_liq(canopy.earth_param_set)  # mol H2O m^-3 (precip volume flux → molar flux)
    tivs = component.time_integrated_vars
    function compute_exp_tendency!(dY, Y, p, t)
        zhou_model = component
        parameters = zhou_model.parameters
        static_inputs = zhou_model.optimal_lai_inputs
        pmodel_parameters = canopy.photosynthesis.parameters
        pmodel_constants = canopy.photosynthesis.constants
        fractional_c3 = canopy.photosynthesis.fractional_c3
        earth_param_set = canopy.earth_param_set

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
            p.canopy.soil_moisture_stress.βm,
            static_inputs.vpd_gs,
        )

        @. p.canopy.biomass.L_opt = compute_L_steady_target(
            Y.canopy.biomass.A0_daily,
            parameters.k,
            Y.canopy.biomass.A0_annual,
            parameters.z,
            static_inputs.GSL,
            parameters.sigma,
            Y.canopy.biomass.precip_annual,
            static_inputs.f0,
            p.drivers.c_co2 * p.drivers.P,  # ca_pa: CO2 partial pressure (Pa)
            p.canopy.biomass.OptVars.χ,
            static_inputs.vpd_gs,
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
        @. dY.canopy.biomass.LAI = apply_time_reduction(
            p.canopy.biomass.L_opt,
            Y.canopy.biomass.LAI,
            tivs.LAI.reduction,
        )
    end
    return compute_exp_tendency!
end

struct TwoComponentResNetLAIModel{
    FT,
    PLPT <: TwoComponentResNetLAIParameters{FT},
    GD,
    FS <: Union{FT, ClimaCore.Fields.Field},
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
    T,
} <: AbstractBiomassModel{FT}
    "Required parameters for the model"
    parameters::PLPT
    "Spatially varying environmental inputs"
    optimal_lai_inputs::GD
    "Prescribed stem area index (m^2 m^-2)"
    SAI::FS
    "Prescribed root area index (m^2 m^-2)"
    RAI::FS
    "Rooting depth parameter (m)"
    rooting_depth::RDTH
    "Canopy height (m)"
    height::HTH
    "Prognostic time-integrated variables (A0_daily, A0_annual, precip_annual, LAI_wood, LAI_herb)"
    time_integrated_vars::T
end

Base.eltype(::TwoComponentResNetLAIModel{FT}) where {FT} = FT

function TwoComponentResNetLAIModel{FT}(
    parameters::TwoComponentResNetLAIParameters{FT},
    optimal_lai_inputs;
    SAI,
    RAI,
    rooting_depth,
    height,
) where {FT <: AbstractFloat}
    seconds_per_day = IP.day(IP.InsolationParameters(FT))
    tau_long_term = FT(365 * 86400.0)
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
            name = :LAI_wood,
            reduction = ClimaLand.RunningMean(
                seconds_per_day / parameters.alpha_wood,
            ),
        ),
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :LAI_herb,
            reduction = ClimaLand.RunningMean(
                seconds_per_day / parameters.alpha_herb,
            ),
        ),
    )
    return TwoComponentResNetLAIModel{
        FT,
        typeof(parameters),
        typeof(optimal_lai_inputs),
        typeof(SAI),
        typeof(rooting_depth),
        typeof(height),
        typeof(tiv),
    }(
        parameters,
        optimal_lai_inputs,
        SAI,
        RAI,
        rooting_depth,
        height,
        tiv,
    )
end

ClimaLand.auxiliary_vars(model::TwoComponentResNetLAIModel) =
    (:area_index, :OptVars, :L_opt_wood, :L_opt_herb)
ClimaLand.auxiliary_types(model::TwoComponentResNetLAIModel{FT}) where {FT} = (
    NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    NamedTuple{(:A0, :χ), Tuple{FT, FT}},
    FT,
    FT,
)
ClimaLand.auxiliary_domain_names(::TwoComponentResNetLAIModel) =
    (:surface, :surface, :surface, :surface)

ClimaLand.prognostic_vars(m::TwoComponentResNetLAIModel) =
    ClimaLand.time_integrated_prognostic_vars(m.time_integrated_vars)
ClimaLand.prognostic_types(m::TwoComponentResNetLAIModel) =
    ClimaLand.time_integrated_prognostic_types(m.time_integrated_vars)
ClimaLand.prognostic_domain_names(m::TwoComponentResNetLAIModel) =
    ClimaLand.time_integrated_prognostic_domain_names(m.time_integrated_vars)

function update_biomass!(
    p,
    Y,
    t,
    component::TwoComponentResNetLAIModel{FT},
    canopy,
) where {FT}
    (; SAI, RAI, parameters, optimal_lai_inputs) = component
    @. p.canopy.biomass.area_index.stem = SAI
    @. p.canopy.biomass.area_index.root = RAI

    # Eco-climatic indices for active optical greenness gating
    p_mm_ann = @. Y.canopy.biomass.precip_annual / FT(55.55)
    si = compute_seasonality_index.(optimal_lai_inputs.GSL)
    ai = compute_aridity_index.(
        Y.canopy.biomass.precip_annual,
        Y.canopy.biomass.A0_annual,
        optimal_lai_inputs.vpd_gs,
    )
    rainforest_idx = @. clamp((FT(1) - si) * tanh(p_mm_ann / FT(1500.0)), FT(0), FT(1))
    evergreen_idx = @. clamp((FT(1) - si) * tanh(ai), FT(0), FT(1))
    deciduous_idx = @. clamp(si * (FT(1) - rainforest_idx), FT(0), FT(1))
    dryland_pheno_idx = compute_dryland_phenology_index.(rainforest_idx, evergreen_idx, p_mm_ann)

    # Thermal and cryogenic phenology gates
    t_air_k = p.drivers.T
    f_freeze = compute_thermal_freeze_gate.(t_air_k, parameters.theta_freeze, parameters.kappa_cold)
    pheno_gate = @. FT(1.0) - deciduous_idx * (FT(1.0) - f_freeze)
    cryo_gate = @. FT(1.0) - evergreen_idx * (FT(1.0) - f_freeze)

    # Root-zone moisture stress
    soil_stress = p.canopy.soil_moisture_stress.βm

    # Active optical green leaf fraction (accounts for drought curing & dormancy)
    f_green = compute_hydro_optical_greenness.(
        deciduous_idx,
        evergreen_idx,
        dryland_pheno_idx,
        pheno_gate,
        cryo_gate,
        f_freeze,
        soil_stress,
        parameters.green_min,
        parameters.hydro_weight,
    )

    # Dual-trigger dormancy gate for structural baseline floor
    f_drought = compute_drought_dormancy_gate.(soil_stress, dryland_pheno_idx, parameters.beta_dry)
    f_dormant = compute_dual_trigger_dormancy.(f_freeze, f_drought)
    struct_floor = compute_structural_canopy_floor.(
        rainforest_idx,
        evergreen_idx,
        deciduous_idx,
        ai,
        f_dormant,
        parameters.theta_ai_scale,
        parameters.enf_retention_scale,
        parameters.enf_cold_gain,
        t_air_k,
    )

    @. p.canopy.biomass.area_index.leaf =
        max(struct_floor, (Y.canopy.biomass.LAI_wood + Y.canopy.biomass.LAI_herb) * f_green)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
end

function ClimaLand.make_compute_exp_tendency(
    component::TwoComponentResNetLAIModel{FT},
    canopy,
) where {FT}
    ρ_m_liq = LP.ρ_m_liq(canopy.earth_param_set)
    tivs = component.time_integrated_vars
    function compute_exp_tendency!(dY, Y, p, t)
        model = component
        parameters = model.parameters
        static_inputs = model.optimal_lai_inputs
        pmodel_parameters = canopy.photosynthesis.parameters
        pmodel_constants = canopy.photosynthesis.constants
        fractional_c3 = canopy.photosynthesis.fractional_c3
        earth_param_set = canopy.earth_param_set

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
            p.canopy.soil_moisture_stress.βm,
            static_inputs.vpd_gs,
        )

        # Compute total target LAI
        L_opt_total = compute_L_steady_target.(
            Y.canopy.biomass.A0_daily,
            parameters.k,
            Y.canopy.biomass.A0_annual,
            parameters.z,
            static_inputs.GSL,
            parameters.sigma,
            Y.canopy.biomass.precip_annual,
            static_inputs.f0,
            p.drivers.c_co2 .* p.drivers.P,
            p.canopy.biomass.OptVars.χ,
            static_inputs.vpd_gs,
        )

        # Dynamic ResNet-PINN carbon allocation multiplier
        A0_norm = clamp.(p.canopy.biomass.OptVars.A0 ./ max.(Y.canopy.biomass.A0_daily, eps(FT)), FT(0), FT(2))
        vpd_norm = clamp.(static_inputs.vpd_gs ./ FT(2000.0), FT(0), FT(2))
        sw_norm = clamp.(p.canopy.soil_moisture_stress.βm, FT(0), FT(1))
        alloc_mult = compute_resnet_carbon_allocation.(Ref(model.parameters), A0_norm, vpd_norm, sw_norm)

        # Eco-climatic partitioning
        p_mm_ann = @. Y.canopy.biomass.precip_annual / FT(55.55)
        ai = compute_aridity_index.(Y.canopy.biomass.precip_annual, Y.canopy.biomass.A0_annual, static_inputs.vpd_gs)
        si = compute_seasonality_index.(static_inputs.GSL)
        f_wood = compute_woody_fraction.(ai, si)

        @. p.canopy.biomass.L_opt_wood = L_opt_total * alloc_mult * f_wood
        @. p.canopy.biomass.L_opt_herb = L_opt_total * alloc_mult * (FT(1) - f_wood)

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
        @. dY.canopy.biomass.precip_annual = apply_time_reduction(
            -(p.drivers.P_liq + p.drivers.P_snow) * ρ_m_liq,
            Y.canopy.biomass.precip_annual,
            tivs.precip_annual.reduction,
        )
        @. dY.canopy.biomass.LAI_wood = apply_time_reduction(
            p.canopy.biomass.L_opt_wood,
            Y.canopy.biomass.LAI_wood,
            tivs.LAI_wood.reduction,
        )
        @. dY.canopy.biomass.LAI_herb = apply_time_reduction(
            p.canopy.biomass.L_opt_herb,
            Y.canopy.biomass.LAI_herb,
            tivs.LAI_herb.reduction,
        )
    end
    return compute_exp_tendency!
end

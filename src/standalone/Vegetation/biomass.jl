import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import Interpolations: Constant
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
using DocStringExtensions

export prescribed_lai_era5, prescribed_lai_modis, PrescribedAreaIndices



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
by passing interpolation_method = Interpolations.Linear().
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
                          time_interpolation_method =
                                        LinearInterpolation(PeriodicCalendar()))
                          regridder_type = :InterpolationsRegridder,
                          interpolation_method = Interpolations.Constant(),
                          context = ClimaComms.context(surface_space))

A helper function which constructs the TimeVaryingInput object for Leaf Area
Index using MODIS LAI data; requires the
surface_space, the start and stop dates, and the earth_param_set.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported
by passing interpolation_method = Interpolations.Linear().
"""
function prescribed_lai_modis(
    surface_space,
    start_date,
    stop_date;
    time_interpolation_method = LinearInterpolation(),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
    context = ClimaComms.context(surface_space),
)
    modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
        context,
        start_date,
        stop_date,
    )
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
    AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for plant hydraulics models.
"""
abstract type AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(::AbstractBiomassModel) = :biomass

abstract type AbstractAreaIndexModel{FT <: AbstractFloat} end

"""
   PrescribedAreaIndices{FT <:AbstractFloat, F <: AbstractTimeVaryingInput}

A struct containing the area indices of the plants at a specific site;
LAI varies in time, while SAI and RAI are fixed in space and time.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedAreaIndices{
    FT <: AbstractFloat,
    F <: AbstractTimeVaryingInput,
} <: AbstractAreaIndexModel{FT}
    "A function of simulation time `t` giving the leaf area index (LAI; m2/m2)"
    LAIfunction::F
    "The constant stem area index (SAI; m2/m2)"
    SAI::FT
    "The constant root area index (RAI; m2/m2)"
    RAI::FT
end

function PrescribedAreaIndices{FT}(
    LAIfunction::AbstractTimeVaryingInput,
    SAI::FT,
    RAI::FT,
) where {FT <: AbstractFloat}
    PrescribedAreaIndices{FT, typeof(LAIfunction)}(LAIfunction, SAI, RAI)
end

function PrescribedAreaIndices(
    LAIfunction::AbstractTimeVaryingInput,
    toml_dict::CP.AbstractTOMLDict;
    SAI = toml_dict["SAI"],
    RAI = toml_dict["RAI"],
)
    FT = CP.float_type(toml_dict)
    return PrescribedAreaIndices{FT, typeof(LAIfunction)}(LAIfunction, SAI, RAI)
end

"""
    lai_consistency_check(
        n_stem::Int64,
        n_leaf::Int64,
        area_index::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    ) where {FT}

Carries out consistency checks using the area indices supplied and the number of
stem elements being modeled.

Note that it is possible to have a plant with no stem compartments
but with leaf compartments, and that a plant must have leaf compartments
(even if LAI = 0).

Specifically, this checks that:
1. n_leaf > 0
2. if LAI is nonzero or SAI is nonzero, RAI must be nonzero.
3. if SAI > 0, n_stem must be > 0 for consistency. If SAI == 0, n_stem must
be zero.
"""
function lai_consistency_check(
    n_stem::Int64,
    n_leaf::Int64,
    area_index::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
) where {FT}
    @assert n_leaf > 0
    if area_index.leaf > eps(FT) || area_index.stem > eps(FT)
        @assert area_index.root > eps(FT)
    end
    # If there SAI is zero, there should be no stem compartment
    if area_index.stem < eps(FT)
        @assert n_stem == FT(0)
    else
        # if SAI is > 0, n_stem should be > 0 for consistency
        @assert n_stem > 0
    end

end

"""
    struct PrescribedBiomassModel{FT, PSAI <: PrescribedAreaIndices, RDTH <: Union{FT, ClimaCore.Fields.Field}} <: AbstractBiomassModel{FT}

A prescribed biomass model where LAI, SAI, RAI, and rooting depth are prescribed.

In  global run with patches
of bare soil, you can "turn off" the canopy model (to get zero root extraction, zero absorption and
emission, zero transpiration and sensible heat flux from the canopy), by setting:
- n_leaf = 1
- n_stem = 0
- LAI = SAI = RAI = 0.

If run with PlantHydraulics, this must be consistent with the plant hydraulics model:a
 plant model can have leaves but no stem, but not vice versa. If n_stem = 0, SAI must be zero.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedBiomassModel{
    FT,
    PSAI <: PrescribedAreaIndices,
    RDTH <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractBiomassModel{FT}
    "The area index model for LAI, SAI, RAI"
    ai_parameterization::PSAI
    "Rooting depth parameter (m) - a characteristic depth below which 1/e of the root mass lies"
    rooting_depth::RDTH
end

auxiliary_vars(model::PrescribedBiomassModel) = (:area_index,)
"""
    ClimaLand.auxiliary_types(model::PlantHydraulicsModel{FT}) where {FT}

Defines the auxiliary types for the PlantHydraulicsModel.
"""
ClimaLand.auxiliary_types(model::PrescribedBiomassModel{FT}) where {FT} =
    (NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},)
ClimaLand.auxiliary_domain_names(::PrescribedBiomassModel) = (:surface,)

function clip(x::FT, threshold::FT) where {FT}
    x > threshold ? x : FT(0)
end

"""
    set_canopy_prescribed_field!(component::PlantHydraulics{FT},
                                 p,
                                 t,
                                 ) where {FT}


Sets the canopy prescribed fields pertaining to the PlantHydraulics
component (the area indices) with their values at time t.

Note that we clip all values of LAI below 0.05 to zero.
This is because we currently run into issues when LAI is
of order eps(FT) in the SW radiation code.
Please see Issue #644
or PR #645 for details.
For now, this clipping is similar to what CLM and NOAH MP do.
"""
function ClimaLand.Canopy.set_canopy_prescribed_field!(
    component::PrescribedBiomassModel{FT},
    p,
    t,
) where {FT}
    (; LAIfunction, SAI, RAI) = component.ai_parameterization
    evaluate!(p.canopy.biomass.area_index.leaf, LAIfunction, t)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    @. p.canopy.biomass.area_index.stem = SAI
    @. p.canopy.biomass.area_index.root = RAI
end

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
    function make_biomass_callback(
        ::PrescribedBiomassModel,
        canopy,
    )

"""
function make_biomass_callback(
    model::PrescribedBiomassModel,
    canopy; kwargs...
)
    return FrequencyBasedCallback(
        dt,         # period of this callback
        start_date; # initial datetime, UTC
        func = (integrator) -> evaluate!(
            integrator.p.canopy.biomass.area_index.leaf,
            canopy.biomass.ai_parameterization.LAI,
            integrator.t,
        ),
    )
end

"""
    get_model_callbacks(component::AbstractCanopyComponent, canopy; kwargs...)

Creates the biomass callback and returns it as a single element tuple of model callbacks.
"""
function get_model_callbacks(
    component::PrescribedBiomassModel,
    canopy;
    kwargs...
)
    lai_cb = make_biomass_callback(component, canopy)
    return (lai_cb,)
end

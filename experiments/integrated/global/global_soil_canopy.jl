import SciMLBase
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations
using Insolation

import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using CairoMakie
using Statistics
using Dates
import NCDatasets

regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
context = ClimaComms.context()
outdir = joinpath(pkgdir(ClimaLand), "experiments/integrated/global/plots")
!ispath(outdir) && mkpath(outdir)

device_suffix =
    typeof(ClimaComms.context().device) <: ClimaComms.CPUSingleThreaded ?
    "cpu" : "gpu"

FT = Float64
earth_param_set = LP.LandParameters(FT)
radius = FT(6378.1e3);
depth = FT(50)
nelements = (50, 10)
domain = ClimaLand.Domains.SphericalShell(;
    radius = radius,
    depth = depth,
    nelements = nelements,
    npolynomial = 1,
    dz_tuple = FT.((10.0, 0.1)),
);
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface

ref_time = DateTime(2021);
t_start = 0.0

# Forcing data
era5_artifact_path =
    ClimaLand.Artifacts.era5_land_forcing_data2021_folder_path(; context)
precip = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
    "rf",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
    file_reader_kwargs = (; preprocess_func = (data) -> -data / 3600,),
)

snow_precip = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
    "sf",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
    file_reader_kwargs = (; preprocess_func = (data) -> -data / 3600,),
)

u_atmos = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
    "ws",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
)
q_atmos = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
    "q",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
)
P_atmos = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
    "sp",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
)

T_atmos = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
    "t2m",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
)
h_atmos = FT(10);

atmos = PrescribedAtmosphere(
    precip,
    snow_precip,
    T_atmos,
    u_atmos,
    q_atmos,
    P_atmos,
    ref_time,
    h_atmos,
    earth_param_set,
);

# Prescribed radiation
SW_d = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
    "ssrd",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
    file_reader_kwargs = (; preprocess_func = (data) -> data / 3600,),
)
LW_d = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
    "strd",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
    file_reader_kwargs = (; preprocess_func = (data) -> data / 3600,),
)

function zenith_angle(
    t,
    ref_time;
    latitude = ClimaCore.Fields.coordinate_field(surface_space).lat,
    longitude = ClimaCore.Fields.coordinate_field(surface_space).long,
    insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
) where {FT}
    # This should be time in UTC
    current_datetime = ref_time + Dates.Second(round(t))

    # Orbital Data uses Float64, so we need to convert to our sim FT
    d, δ, η_UTC =
        FT.(
            Insolation.helper_instantaneous_zenith_angle(
                current_datetime,
                ref_time,
                insol_params,
            )
        )

    Insolation.instantaneous_zenith_angle.(d, δ, η_UTC, longitude, latitude).:1
end
radiation =
    PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time; θs = zenith_angle);

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments/integrated/global/global_parameters.jl",
    ),
)
soil_args = (domain = domain, parameters = soil_params)
soil_model_type = Soil.EnergyHydrology{FT}

# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

# soil microbes args
Csom = (z, t) -> eltype(z)(5.0)

# Set the soil CO2 BC to being atmospheric CO2
soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
soilco2_sources = (Soil.Biogeochemistry.MicrobeProduction{FT}(),)

soilco2_boundary_conditions = (; top = soilco2_top_bc, bottom = soilco2_bot_bc)

soilco2_drivers = Soil.Biogeochemistry.SoilDrivers(
    Soil.Biogeochemistry.PrognosticMet{FT}(),
    Soil.Biogeochemistry.PrescribedSOC{FT}(Csom),
    atmos,
)

soilco2_args = (;
    boundary_conditions = soilco2_boundary_conditions,
    sources = soilco2_sources,
    domain = domain,
    parameters = soilco2_ps,
    drivers = soilco2_drivers,
)

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
    radiative_transfer = Canopy.TwoStreamModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
    energy = Canopy.BigLeafEnergyModel{FT},
)
# Individual Component arguments
# Set up autotrophic respiration
autotrophic_respiration_args =
    (; parameters = Canopy.AutotrophicRespirationParameters(FT))
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = Canopy.TwoStreamParameters(
        FT;
        Ω,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
    )
)
# Set up conductance
conductance_args = (; parameters = Canopy.MedlynConductanceParameters(FT; g1))
# Set up photosynthesis
photosynthesis_args = (;
    parameters = Canopy.FarquharParameters(FT, Canopy.C3(); Vcmax25 = Vcmax25)
)
# Set up plant hydraulics
# Not ideal
LAIfunction = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
    "lai",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type,
    file_reader_kwargs = (;
        preprocess_func = (data) -> data > 0.05 ? data : 0.0,
    ),
)
ai_parameterization = Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    root_distribution = root_distribution,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
)
plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

# Canopy component args
canopy_component_args = (;
    autotrophic_respiration = autotrophic_respiration_args,
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
    energy = energy_args,
)

# Other info needed
shared_params = Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (;
    parameters = shared_params,
    domain = ClimaLand.obtain_surface_domain(domain),
)

# Integrated plant hydraulics and soil model
land_input = (atmos = atmos, radiation = radiation)
land = SoilCanopyModel{FT}(;
    soilco2_type = soilco2_type,
    soilco2_args = soilco2_args,
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
)

Y, p, cds = initialize(land)

t0 = 0.0
dt = 180.0
tf = 3600

init_soil(ν, θ_r) = θ_r + (ν - θ_r) / 2
Y.soil.ϑ_l .= init_soil.(ν, θ_r)
Y.soil.θ_i .= FT(0.0)
T = FT(276.85)
ρc_s =
    Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        soil_params.ρc_ds,
        soil_params.earth_param_set,
    )
Y.soil.ρe_int .=
    Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T,
        soil_params.earth_param_set,
    )
Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
# On the Caltech cluster, this is not permitted, so we allowscalar.
# This operation is fine on clima cluster
Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
evaluate!(Y.canopy.energy.T, atmos.T, t0)

set_initial_cache! = make_set_initial_cache(land)
exp_tendency! = make_exp_tendency(land);
imp_tendency! = ClimaLand.make_imp_tendency(land);
tendency_jacobian! = ClimaLand.make_tendency_jacobian(land);
set_initial_cache!(p, Y, t0)
stepper = CTS.ARS343()
norm_condition = CTS.MaximumError(FT(1e-8))
conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

# set up jacobian info
jac_kwargs =
    (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = tendency_jacobian!)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
)

saveat = [t0, tf]
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
updateat = Array(t0:(3600 * 3):tf)
updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
@time sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = cb,
    adaptive = false,
    saveat = saveat,
)

if device_suffix == "cpu"
    longpts = range(-180.0, 180.0, 101)
    latpts = range(-90.0, 90.0, 101)
    hcoords = [
        ClimaCore.Geometry.LatLongPoint(lat, long) for long in longpts,
        lat in reverse(latpts)
    ]
    remapper = ClimaCore.Remapping.Remapper(surface_space, hcoords)
    S = ClimaLand.Domains.top_center_to_surface(
        (sol.u[end].soil.ϑ_l .- θ_r) ./ (ν .- θ_r),
    )
    S_ice = ClimaLand.Domains.top_center_to_surface(sol.u[end].soil.θ_i ./ ν)
    T_soil = ClimaLand.Domains.top_center_to_surface(sv.saveval[end].soil.T)
    SW = sv.saveval[end].drivers.SW_d

    GPP = sv.saveval[end].canopy.photosynthesis.GPP .* 1e6
    T_canopy = sol.u[end].canopy.energy.T
    fields = [S, S_ice, T_soil, GPP, T_canopy, SW]
    titles = [
        "Effective saturation",
        "Effective ice saturation",
        "Temperature (K) - Soil",
        "GPP",
        "Temperature (K) - Canopy",
        "Incident SW",
    ]
    plotnames = ["S", "Sice", "temp", "gpp", "temp_canopy", "sw"]
    mask_top = ClimaLand.Domains.top_center_to_surface(soil_params_mask)
    mask_remap = ClimaCore.Remapping.interpolate(remapper, mask_top)
    for (id, x) in enumerate(fields)
        title = titles[id]
        plotname = plotnames[id]
        x_remap = ClimaCore.Remapping.interpolate(remapper, x)

        fig = Figure(size = (600, 400))
        ax = Axis(
            fig[1, 1],
            xlabel = "Longitude",
            ylabel = "Latitude",
            title = title,
        )
        clims = extrema(x_remap)
        CairoMakie.heatmap!(
            ax,
            longpts,
            latpts,
            oceans_to_value.(x_remap, mask_remap, 0.0),
            colorrange = clims,
        )
        Colorbar(fig[:, end + 1], colorrange = clims)
        outfile = joinpath(outdir, "$plotname.png")
        CairoMakie.save(outfile, fig)
    end
end

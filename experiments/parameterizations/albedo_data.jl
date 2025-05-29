import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities
import ClimaUtilities.Regridders: InterpolationsRegridder, regrid
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
using ClimaCore
using Plots
using Insolation
import ClimaParams as CP
using Dates
using Statistics
using NCDatasets
using DelimitedFiles


using ClimaLand
import ClimaLand.Parameters as LP
start_date = DateTime(2008)
time_interpolation_method = LinearInterpolation()
regridder_type = :InterpolationsRegridder
nelements = (101, 15)
FT = Float32
earth_param_set = LP.LandParameters(FT)
thermo_params = LP.thermodynamic_parameters(earth_param_set);
domain = ClimaLand.ModelSetup.global_domain(FT; nelements, mask_threshold = 0.99)
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface
outpath = "experiments/parameterizations"

# Spatially varying canopy parameters from CLM
clm_parameters = ClimaLand.clm_canopy_parameters(surface_space)
    (;
     is_c3,
     α_PAR_leaf,
     τ_PAR_leaf,
     α_NIR_leaf,
     τ_NIR_leaf,
     ) = clm_parameters



spatially_varying_soil_params =
    ClimaLand.default_spatially_varying_soil_parameters(
        subsurface_space,
        surface_space,
        FT,
    )

(; ν, ν_ss_om, ν_ss_quartz, ν_ss_gravel, hydrology_cm, K_sat, θ_r, PAR_albedo_wet,
 NIR_albedo_wet,
 PAR_albedo_dry,
 NIR_albedo_dry) = 
     spatially_varying_soil_params

function fillmissing(x)
    tmp = x isa Missing ? 0.0 : x
end

soc = SpaceVaryingInput(
    "../ClimaArtifacts/soilgrids/soilgrids_nc/soilgrids_nc/soc_0-5cm_mean_5000.nc",
    "Band1",
    surface_space;
    regridder_type,
    compose_function = (x) -> fillmissing.(x),
)
bdod = SpaceVaryingInput(
    "../ClimaArtifacts/soilgrids/soilgrids_nc/soilgrids_nc/bdod_0-5cm_mean_5000.nc",
    "Band1",
    surface_space;
    regridder_type,
    compose_function = (x) -> fillmissing.(x),
)
clay = SpaceVaryingInput(
    "../ClimaArtifacts/soilgrids/soilgrids_nc/soilgrids_nc/clay_0-5cm_mean_5000.nc",
    "Band1",
    surface_space;
    regridder_type,
    compose_function = (x) -> fillmissing.(x),
)
silt = SpaceVaryingInput(
    "../ClimaArtifacts/soilgrids/soilgrids_nc/soilgrids_nc/silt_0-5cm_mean_5000.nc",
    "Band1",
    surface_space;
    regridder_type,
    compose_function = (x) -> fillmissing.(x),
)
sand = SpaceVaryingInput(
    "../ClimaArtifacts/soilgrids/soilgrids_nc/soilgrids_nc/sand_0-5cm_mean_5000.nc",
    "Band1",
    surface_space;
    regridder_type,
    compose_function = (x) -> fillmissing.(x),
)
cfvo = SpaceVaryingInput(
    "../ClimaArtifacts/soilgrids/soilgrids_nc/soilgrids_nc/cfvo_0-5cm_mean_5000.nc",
    "Band1",
    surface_space;
    regridder_type,
    compose_function = (x) -> fillmissing.(x),
)

function replace_missing_with_mean(x)
    absent = parent(x) .≈ 0.0
    parent(x)[absent] .= mean(parent(x)[.!absent])
end
replace_missing_with_mean(soc)
replace_missing_with_mean(bdod)
replace_missing_with_mean(clay)
replace_missing_with_mean(silt)
replace_missing_with_mean(sand)
replace_missing_with_mean(cfvo)

modis_lai_artifact_path = ClimaLand.Artifacts.modis_lai_forcing_data_path()
modis_lai_ncdata_path =
    joinpath(modis_lai_artifact_path, "Yuan_et_al_2008_1x1.nc")


LAIfunction = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date;
    time_interpolation_method = time_interpolation_method,
)
LAI_data = NCDataset(modis_lai_ncdata_path)
maxLAI = maximum(LAI_data["lai"][:,:,:], dims = 3)[:,:,1]
rg = InterpolationsRegridder(
           surface_space,
       )
maxLAI_field = regrid(rg, maxLAI,(LAI_data["lon"][:], LAI_data["lat"][:]) )
close(LAI_data)
era5_path = "../ClimaArtifacts/era5_soil_albedo_data/soil_a_data_2008_1.0x1.0_inst.nc";
snow_present = TimeVaryingInput(
    era5_path,
    "sd",
    surface_space;
    start_date,
    regridder_type,
    method = time_interpolation_method,
    compose_function = (sd) -> sd .> 0.0001,
)

θ = TimeVaryingInput(
    era5_path,
    "swvl1",
    surface_space;
    start_date,
    regridder_type,
    method = time_interpolation_method,
)
#=
vpd(Td, T, P; params = earth_param_set) =
    ClimaLand.vapor_pressure_deficit.(T, P, ClimaLand.specific_humidity_from_dewpoint.(Td, T, P, params), thermo_params)
vpd_atmos = TimeVaryingInput(
        [era5_path, era5_path, era5_path],
        ["d2m", "t2m", "sp"],
        surface_space;
        start_date,
        regridder_type,
        compose_function = vpd,
        method = time_interpolation_method,
)
T_atmos = TimeVaryingInput(
        era5_path,
        "t2m",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
)

P_atmos = TimeVaryingInput(
        era5_path,
        "sp",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
)
=#
era5_average_path = "../ClimaArtifacts/era5_soil_albedo_data/averaged_atmospheric_data.nc"
vpd_atmos = TimeVaryingInput(
    era5_average_path,
    "mean_vpd",
    surface_space;
    start_date,
    regridder_type,
    method = time_interpolation_method,
)
T_atmos = TimeVaryingInput(
        era5_average_path,
        "mean_t2m",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
)

P_atmos = TimeVaryingInput(
    era5_average_path,
    "mean_sp",
    surface_space;
    start_date,
    regridder_type,
    method = time_interpolation_method,
)

land_mask(lsm, cl) = cl < 0.1 && lsm >= 0.5 # 1 if land
era5_mask_tvi = TimeVaryingInput(
    [era5_path, era5_path],
    ["lsm", "cl"],
    surface_space;
    start_date,
    regridder_type,
    compose_function = (lsm, cl) -> land_mask.(lsm, cl),
    method = time_interpolation_method,
)

era5_ncdata_path = "../ClimaArtifacts/era5_soil_albedo_data/soil_a_data_2008_1.0x1.0_accum.nc";

compute_albedo(SSR, SSRD) = 1 - SSR/max(SSRD, eps(Float32))
land_albedo = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["ssr", "ssrd"],
        surface_space;
        start_date,
        regridder_type,
        compose_function = (ssr, ssrd) -> compute_albedo.(ssr, ssrd),
        method = time_interpolation_method,
)
diffuse_fraction = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["ssrd", "fdir"],
        surface_space;
        start_date,
        regridder_type,
        compose_function = (ssrd, fdir) -> (ssrd .- fdir)./max.(ssrd, eps(Float32)),
        method = time_interpolation_method,
)

ssrd = TimeVaryingInput(
        era5_ncdata_path,
        "ssrd",
        surface_space;
        start_date,
        regridder_type,
        method = time_interpolation_method,
)



albedo_field = ClimaCore.Fields.zeros(surface_space)
θ_field = ClimaCore.Fields.zeros(surface_space)
snow_field = ClimaCore.Fields.zeros(surface_space)
ssrd_field = ClimaCore.Fields.zeros(surface_space)
mask = ClimaCore.Fields.zeros(surface_space)
era5_mask = ClimaCore.Fields.zeros(surface_space)
lat = ClimaCore.Fields.coordinate_field(surface_space).lat
lon=  ClimaCore.Fields.coordinate_field(surface_space).long
N = 0
field = ClimaCore.Fields.zeros(surface_space)


function zenith_angle(
    current_datetime, latitude::FT, longitude::FT, start_date, insol_params) where {FT}

        # Orbital Data uses Float64, so we need to convert to our sim FT
        d, δ, η_UTC =
            FT.(
                Insolation.helper_instantaneous_zenith_angle(
                    current_datetime,
                    start_date,
                    insol_params,
                )
            )

        Insolation.instantaneous_zenith_angle.(
            d,
            δ,
            η_UTC,
            longitude,
            latitude,
        ).:1
    end

for (i, t) in enumerate(θ.data_handler.available_dates)
    evaluate!(albedo_field, land_albedo, t)
    evaluate!(θ_field, θ, t)
    evaluate!(ssrd_field, ssrd, t)
    evaluate!(snow_field, snow_present, t)
    evaluate!(era5_mask, era5_mask_tvi, t)
    μ = max.(cos.(zenith_angle.(t, lat, lon, start_date, Ref(earth_param_set.insol_params))), 0)
    @. mask = !(
        (snow_field ≈ 1.0) ||
        (albedo_field < 0.05) ||
        (albedo_field > 0.95) ||
        (μ < 0.09) ||
        (era5_mask <0.99)
    )
    array_mask = parent(mask)[:] .≈ 1
    N_i = sum(array_mask)
    @show i
    @show N_i
    N += sum(array_mask)
end


net_features = zeros(N, 23)
nonnormed_input = zeros(N, 4)
clm_albedos = zeros(N, 8)
metadata = zeros(N, 2)

id = 1
for (i, t) in enumerate(θ.data_handler.available_dates)
    evaluate!(albedo_field, land_albedo, t)
    evaluate!(θ_field, θ, t)
    evaluate!(snow_field, snow_present, t)
    evaluate!(era5_mask, era5_mask_tvi, t)
    μ = max.(cos.(zenith_angle.(t, lat, lon, start_date, Ref(earth_param_set.insol_params))), 0)
    @. mask = !(
        (snow_field ≈ 1.0) ||
        (albedo_field < 0.05) ||
        (albedo_field > 0.95) ||
        (μ < 0.09) ||
        (era5_mask <0.99)
    )
    array_mask = parent(mask)[:] .≈ 1
    N_i = sum(array_mask)
    
    albedo_i = parent(albedo_field)[:][array_mask]
    ν_i = parent(ν)[end, :, :, :, :][array_mask]
    ν_ss_om_i = parent(ν_ss_om)[end, :, :, :, :][array_mask]
    ν_ss_quartz_i = parent(ν_ss_quartz)[end, :, :, :, :][array_mask]
    ν_ss_gravel_i = parent(ν_ss_gravel)[end, :, :, :, :][array_mask]
    vg_α_i = parent(hydrology_cm.α)[end, :, :, :, :][array_mask]
    vg_n_i = parent(hydrology_cm.n)[end, :, :, :, :][array_mask]
    K_sat_i = parent(K_sat)[end, :, :, :, :][array_mask]
    θ_i = parent(θ_field)[:][array_mask]
    # net features
    net_features[id:(id + N_i - 1), 1] .= min.(θ_i ./ ν_i, 1.0)
    net_features[id:(id + N_i - 1), 2] .= 1.0 ./ vg_α_i
    net_features[id:(id + N_i - 1), 3] .= vg_n_i .- 1
    net_features[id:(id + N_i - 1), 4] .= log10.(K_sat_i)
    net_features[id:(id + N_i - 1), 5] .= ν_i
    net_features[id:(id + N_i - 1), 6] .= parent(cfvo)[:][array_mask]
    net_features[id:(id + N_i - 1), 7] .= parent(clay)[:][array_mask]
    net_features[id:(id + N_i - 1), 8] .= parent(sand)[:][array_mask]
    net_features[id:(id + N_i - 1), 9] .= parent(silt)[:][array_mask]
    net_features[id:(id + N_i - 1), 10] .= parent(bdod)[:][array_mask]
    net_features[id:(id + N_i - 1), 11] .= log10.(parent(soc)[:][array_mask])
    net_features[id:(id + N_i - 1), 12] .= parent(maxLAI_field)[:][array_mask]
    net_features[id:(id + N_i - 1), 13] .= ν_ss_om_i
    net_features[id:(id + N_i - 1), 14] .= ν_ss_gravel_i
    net_features[id:(id + N_i - 1), 15] .= ν_ss_quartz_i
    
    
    evaluate!(field, LAIfunction, t)
    net_features[id:(id + N_i - 1), 16] .= parent(field)[:][array_mask]
    evaluate!(field, P_atmos, t)
    net_features[id:(id + N_i - 1), 17] .= parent(field)[:][array_mask]
    evaluate!(field, T_atmos, t)
    net_features[id:(id + N_i - 1), 18] .= parent(field)[:][array_mask]
    evaluate!(field, vpd_atmos, t)
    net_features[id:(id + N_i - 1), 19] .= parent(field)[:][array_mask]
    evaluate!(field, diffuse_fraction, t)
    net_features[id:(id + N_i - 1), 20] .= parent(field)[:][array_mask]
    field = is_c3;
    net_features[id:(id + N_i - 1), 21] .= parent(field)[:][array_mask]
    net_features[id:(id + N_i - 1), 22] .= parent(μ)[:][array_mask]
    net_features[id:(id + N_i - 1), 23] .= albedo_i
    # non-normed input
    nonnormed_input[id:(id + N_i - 1), 1] .= θ_i
    
    evaluate!(field, LAIfunction, t)
    nonnormed_input[id:(id + N_i - 1), 2] .= parent(field)[:][array_mask]
    
    evaluate!(field, diffuse_fraction, t)
    nonnormed_input[id:(id + N_i - 1), 3] .= parent(field)[:][array_mask]
    
    nonnormed_input[id:(id + N_i - 1), 4] .= parent(μ)[:][array_mask]

    # clm_albedos
    field = α_PAR_leaf;
    clm_albedos[id:(id + N_i - 1), 1] .= parent(field)[:][array_mask]
    field = α_NIR_leaf;
    clm_albedos[id:(id + N_i - 1), 2] .= parent(field)[:][array_mask]
    field = τ_PAR_leaf;
    clm_albedos[id:(id + N_i - 1), 3] .= parent(field)[:][array_mask]
    field = τ_NIR_leaf;
    clm_albedos[id:(id + N_i - 1), 4] .= parent(field)[:][array_mask]
    field =  PAR_albedo_wet;
    clm_albedos[id:(id + N_i - 1), 5] .= parent(field)[:][array_mask]
    field =  PAR_albedo_dry;
    clm_albedos[id:(id + N_i - 1), 6] .= parent(field)[:][array_mask]
    field =  NIR_albedo_wet;
    clm_albedos[id:(id + N_i - 1), 7] .= parent(field)[:][array_mask]
    field =  NIR_albedo_dry;
    clm_albedos[id:(id + N_i - 1), 8] .= parent(field)[:][array_mask]
    
    # met data
    metadata[id:(id + N_i - 1), 1] .= parent(lat)[:][array_mask]
    metadata[id:(id + N_i - 1), 2] .= parent(lon)[:][array_mask]
    id = id + N_i
end

# normalize the appropriate fields
net_features[:,1:end-1] .= (net_features[:,1:end-1] .- mean(net_features[:,1:end-1], dims = 1)) ./ std(net_features[:,1:end-1], dims = 1);

open(joinpath(outpath,"net_features.csv"); write=true) do f
    write(f, "θ, 1/α, n-1, log10(K_sat), ν, cfvo, clay, sand, silt, bdod, soc, maxLAI, ν_ss_om, ν_ss_gravel, ν_ss_quartz, LAI, P, T, vpd, diffuse_frac, c3,μ, albedo\n")
         writedlm(f, net_features,',')
end
open(joinpath(outpath, "nonnormed_input.csv"); write=true) do f
    write(f, "θ,  LAI,diffuse_frac,μ\n")
         writedlm(f, nonnormed_input,',')
end

open(joinpath(outpath, "clm_albedos.csv"); write=true) do f
    write(f, "α_PAR_leaf, α_NIR_leaf, τ_PAR_leaf, τ_NIR_leaf,PAR_α_wet, PAR_α_dry, NIR_α_wet, NIR_α_dry \n")
    writedlm(f, clm_albedos, ',')
end
open(joinpath(outpath, "metadata.csv"); write=true) do f
    write(f, "lat, lon \n")
    writedlm(f, metadata, ',')
end

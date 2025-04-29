# Try making heatmaps of variables - are they reasonable? do they correlate with soil type?
# Look at CLM prediction - is moisture variable really wrong?

import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
using ClimaCore
using Plots
using Insolation
import ClimaParams as CP
using Dates
using Statistics

using ClimaLand
import ClimaLand.Parameters as LP
start_date = DateTime(2008)
time_interpolation_method = LinearInterpolation()
regridder_type = :InterpolationsRegridder
nelements = (101, 15)
FT = Float32
earth_param_set = LP.LandParameters(FT)
thermo_params = LP.thermodynamic_parameters(earth_param_set);
domain = ClimaLand.global_domain(FT; nelements)
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface


spatially_varying_soil_params =
    ClimaLand.default_spatially_varying_soil_parameters(
        subsurface_space,
        surface_space,
        FT,
    )

(; ν, ν_ss_om, ν_ss_quartz, ν_ss_gravel, hydrology_cm, K_sat, θ_r) =
    spatially_varying_soil_params

fillmissing(x) = x isa Missing ? -1.0 : x
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
    absent = (parent(x) .≈ -1.0) .|| (parent(x) .≈ 0.0)
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
#=
longpts = range(-180.0, 180.0, 101)
latpts = range(-90.0, 90.0, 101)
hcoords = [
    ClimaCore.Geometry.LatLongPoint(lat, long) for long in longpts,
    lat in latpts
]
remapper = ClimaCore.Remapping.Remapper(surface_space, hcoords)
ClimaCore.Remapping.interpolate(remapper, field)
names = ["albedo", "lai", "θ", "snow", "mask"]
for (i,field) in enumerate([albedo_field, lai_field, θ_field, snow_field, mask])
    name = names[i]
    clims = extrema(field)
    fig = Figure()
    ax1 = Axis(fig[1,1])
    CairoMakie.heatmap!(ax1, longpts, latpts, ClimaCore.Remapping.interpolate(remapper, field), colorrange = clims)
    Colorbar(fig[:, 2], colorrange = clims)
    CairoMakie.save("tmp_$name.png", fig)
end
=#
lat = ClimaCore.Fields.coordinate_field(surface_space).lat
lon=  ClimaCore.Fields.coordinate_field(surface_space).long
N = 0
field = ClimaCore.Fields.zeros(surface_space)

for (i, t) in enumerate(θ.data_handler.available_dates)
    evaluate!(albedo_field, land_albedo, t)
    evaluate!(θ_field, θ, t)
    evaluate!(ssrd_field, ssrd, t)
    evaluate!(snow_field, snow_present, t)
    μ = max.(cos.(zenith_angle.(t, lat, lon, start_date, Ref(earth_param_set.insol_params))), 0)
    @. mask = !(
        (snow_field ≈ 1.0) ||
        (θ_field < eps(Float32)) ||
        (albedo_field < 0.05) ||
        (albedo_field > 0.95) ||
        (μ < eps(FT))
    )
    array_mask = parent(mask)[:] .≈ 1
    N_i = sum(array_mask)
    @show i
    @show N_i
    N += sum(array_mask)
end


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



data = zeros(N, 23)
id = 1
for (i, t) in enumerate(θ.data_handler.available_dates)
    evaluate!(albedo_field, land_albedo, t)
    evaluate!(θ_field, θ, t)
    evaluate!(snow_field, snow_present, t)
    μ = max.(cos.(zenith_angle.(t, lat, lon, start_date, Ref(earth_param_set.insol_params))), 0)
    @. mask = !(
        (snow_field ≈ 1.0) ||
        (θ_field < eps(Float32)) ||
        (albedo_field < 0.05) ||
        (albedo_field > 0.95) ||
        (μ < eps(FT))
    )
    array_mask = parent(mask)[:] .≈ 1
    albedo_i = parent(albedo_field)[:][array_mask]
    ν_i = parent(ν)[end, :, :, :, :][array_mask]
    ν_ss_om_i = parent(ν_ss_om)[end, :, :, :, :][array_mask]
    ν_ss_quartz_i = parent(ν_ss_quartz)[end, :, :, :, :][array_mask]
    ν_ss_gravel_i = parent(ν_ss_gravel)[end, :, :, :, :][array_mask]
    vg_α_i = parent(hydrology_cm.α)[end, :, :, :, :][array_mask]
    vg_n_i = parent(hydrology_cm.n)[end, :, :, :, :][array_mask]
    K_sat_i = parent(K_sat)[end, :, :, :, :][array_mask]
    θ_i = parent(θ_field)[:][array_mask]
    albedo_i = parent(albedo_field)[:][array_mask]
    N_i = sum(array_mask)
    if N_i > 0
        data[id:(id + N_i - 1), 1] .= θ_i
        data[id:(id + N_i - 1), 2] .= 1.0 ./ vg_α_i
        data[id:(id + N_i - 1), 3] .= vg_n_i .- 1
        data[id:(id + N_i - 1), 4] .= log10.(K_sat_i)
        data[id:(id + N_i - 1), 5] .= ν_i
        data[id:(id + N_i - 1), 6] .= parent(cfvo)[:][array_mask]
        data[id:(id + N_i - 1), 7] .= parent(clay)[:][array_mask]
        data[id:(id + N_i - 1), 8] .= parent(sand)[:][array_mask]
        data[id:(id + N_i - 1), 9] .= parent(silt)[:][array_mask]
        data[id:(id + N_i - 1), 10] .= parent(bdod)[:][array_mask]
        data[id:(id + N_i - 1), 11] .= parent(soc)[:][array_mask]
        data[id:(id + N_i - 1), 12] .= ν_ss_om_i
        data[id:(id + N_i - 1), 13] .= ν_ss_quartz_i
        data[id:(id + N_i - 1), 14] .= ν_ss_gravel_i

        evaluate!(field, LAIfunction, t)
        data[id:(id + N_i - 1), 15] .= parent(field)[:][array_mask]
        evaluate!(field, P_atmos, t)
        data[id:(id + N_i - 1), 16] .= parent(field)[:][array_mask]
        evaluate!(field, T_atmos, t)
        data[id:(id + N_i - 1), 17] .= parent(field)[:][array_mask]
        evaluate!(field, vpd_atmos, t)
        data[id:(id + N_i - 1), 18] .= parent(field)[:][array_mask]
        evaluate!(field, diffuse_fraction, t)
        data[id:(id + N_i - 1), 19] .= parent(field)[:][array_mask]
        data[id:(id + N_i - 1), 20] .= parent(μ)[:][array_mask]
        data[id:(id + N_i - 1), 21] .= parent(lat)[:][array_mask]
        data[id:(id + N_i - 1), 22] .= parent(lon)[:][array_mask]
        data[id:(id + N_i - 1), 23] .= albedo_i
    end

    id = id + N_i
end

using DelimitedFiles
open("raw_data.csv"; write=true) do f
    write(f, "# θ, 1/α, n-1, log10(K_sat), ν, cfvo, clay, sand, silt, bdod, soc, νom, ν_sand, ν_cf, LAI, P, T, vpd, diffuse_frac, ν, lat, lon, albedo\n")
         writedlm(f, data,',')
       end

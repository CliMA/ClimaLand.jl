"""
This script runs a standalone canopy model with the P-model for photosynthesis and stomatal
conductance with prescribed atmosphere and soil. 
"""

import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using Insolation
using NCDatasets
import Random

using ClimaLand
using ClimaLand.Domains: Column, Point
using ClimaLand.Soil
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

import ClimaUtilities.OutputPathGenerator: generate_output_path
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeManager: ITime, date
using DelimitedFiles
FluxnetSimulationsExt =
    Base.get_extension(ClimaLand, :FluxnetSimulationsExt).FluxnetSimulationsExt;

climaland_dir = pkgdir(ClimaLand)
const FT = Float64
earth_param_set = LP.LandParameters(FT)

# Read all site-specific domain parameters from the simulation file for the site
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/get_fluxnet_parameters.jl"))
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/get_fluxnet_domain.jl"))

######### Simulation setup #########
site_ID = length(ARGS) >= 1 ? ARGS[1] : "US-NR1"
photo_model = length(ARGS) >= 2 ? ARGS[2] : "pmodel" # "pmodel" or "farquhar"
soil_moisture_stress = length(ARGS) >= 3 ? ARGS[3] : "piecewise" # "piecewise" or "no_sms"

# If these are specified, they override the defaults which are inferred
# from data availability. Note that these are in UTC time 
start_date = nothing
end_date = nothing

saved_var_names = [
    "gpp", "msf", "vcmax25", "clhf" 
]
average_period = :halfhourly
Δt = Float64(600.0)
save_dir = "outputs/fluxnet2015/standalone/$(site_ID)"
####################################

if !isdir(save_dir)
    mkpath(save_dir)
end


function setup_standalone_canopy_model(
    FT, 
    site_ID,
    earth_param_set;
    use_default_param_maps,
    photo_model = "pmodel",
    soil_moisture_stress = "piecewise",
)
    # Get domain info
    (; lat, long, time_offset, h_canopy, atmospheric_sensor_height,
        swc_depths, ts_depths) = get_site_domain(site_ID)

    canopy_domain = Point(; z_sfc = FT(0.0), longlat = (FT(long), lat))

    # Get parameters
    try 
        site_params = get_site_parameters(site_ID)
    catch
        @warn "No site-specific params found in get_fluxnet_parameters.jl for site $(site_ID). Using default parameters specified by parameter maps"
        use_default_param_maps = true
    end

    if use_default_param_maps
        surface_space = canopy_domain.space.surface

        # RT 
        (; Ω, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf) =
            ClimaLand.Canopy.clm_canopy_radiation_parameters(surface_space)

        # medlyn conductance
        g1 = ClimaLand.Canopy.clm_medlyn_g1(surface_space)

        # farquhar photosynthesis
        (; is_c3, Vcmax25) =
            ClimaLand.Canopy.clm_photosynthesis_parameters(surface_space)

        # plant hydraulics
        plant_ν = FT(1.44e-4)
        plant_S_s = FT(1e-2 * 0.0098)
        rooting_depth = ClimaLand.Canopy.clm_rooting_depth(surface_space)
        K_sat_plant = FT(7e-8) # m/s
        ψ63 = FT(-4 / 0.0098) # / MPa to m
        Weibull_param = FT(4) # unitless
        a = FT(0.2 * 0.0098) # 1/m
        conductivity_model = Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
        retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
        SAI = FT(0.0)
        h_stem = FT(0.0)
        h_leaf = FT(h_canopy)
        n_stem = 0
        n_leaf = 1
        zmin = FT(0.0)
        zmax = FT(0.0)
        f_root_to_shoot = FT(3.5)

        # energy balance
        ac_canopy = FT(2.5e3)

        # soil parameters
        (; ν, hydrology_cm, K_sat, θ_r) =
            ClimaLand.Soil.soil_vangenuchten_parameters(surface_space, FT)
        S_s = ClimaCore.Fields.zeros(surface_space) .+ FT(1e-3)
    else
        site_params = get_site_parameters(site_ID)

        # RT 
        Ω = site_params.Ω
        α_PAR_leaf = site_params.α_PAR_leaf
        τ_PAR_leaf = site_params.τ_PAR_leaf
        α_NIR_leaf = site_params.α_NIR_leaf
        τ_NIR_leaf = site_params.τ_NIR_leaf
        G_Function = site_params.G_Function

        # medlyn conductance
        g1 = site_params.g1

        # farquhar photosynthesis
        is_c3 = site_params.is_c3
        Vcmax25 = site_params.Vcmax25

        # plant hydraulics
        plant_ν = site_params.plant_ν
        plant_S_s = site_params.plant_S_s
        rooting_depth = site_params.rooting_depth
        conductivity_model = site_params.conductivity_model
        retention_model = site_params.retention_model
        SAI = site_params.SAI
        h_stem = site_params.h_stem
        h_leaf = site_params.h_leaf
        n_stem = site_params.n_stem
        n_leaf = site_params.n_leaf
        zmin = site_params.zmin
        zmax = site_params.zmax
        f_root_to_shoot = site_params.f_root_to_shoot

        # Energy balance
        ac_canopy = site_params.ac_canopy

        # soil parameters
        ν = site_params.soil_ν
        K_sat = site_params.soil_K_sat
        S_s = site_params.soil_S_s
        soil_vg_n = site_params.soil_vg_n
        soil_vg_α = site_params.soil_vg_α
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n)
        θ_r = site_params.θ_r
    end

    # surface fluxes requires that h_atmos > d_sfc + z_0m ≈ 0.8 * h_canopy
    h_atmos = max(h_canopy, atmospheric_sensor_height[1])

    # simulation time 
    (start_date, end_date) =
        FluxnetSimulationsExt.get_data_dates(site_ID, time_offset)

    # Get the forcing
    (; atmos, radiation, soil) = FluxnetSimulationsExt.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        h_atmos,
        start_date,
        earth_param_set,
        FT;
        construct_soil_driver = true,
        soil_driver_args = (
            α_PAR = FT(0.2),
            α_NIR = FT(0.4),
            ϵ = FT(0.99),
        )
    )

    (; LAI, maxLAI) = FluxnetSimulationsExt.prescribed_LAI_fluxnet(site_ID, start_date)

    # Set up canopy model 
    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m
    shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z0_m,
        z0_b,
        earth_param_set,
    );

    # Set the radiative transfer model
    rt_params = TwoStreamParameters(
        FT;
        G_Function = G_Function,
        α_PAR_leaf = α_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        τ_PAR_leaf = τ_PAR_leaf,
        τ_NIR_leaf = τ_NIR_leaf,
        Ω = Ω,
    )
    radiative_transfer_model = TwoStreamModel{FT}(rt_params);

    # Set up the photosynthesis and stomatal conductance models
    if photo_model == "pmodel"
        # Set the conductance model 
        cond_params = PModelConductanceParameters(Drel = FT(1.6))
        conductance_model = PModelConductance{FT}(cond_params);

        # Set the photosynthesis model (P-model currently only supports C3) 
        photo_params = PModelParameters(
            cstar = FT(0.41),
            β = FT(146),
            ϕc = FT(0.087),
            ϕ0 = FT(NaN),
            ϕa0 = FT(0.352),
            ϕa1 = FT(0.022),
            ϕa2 = FT(-0.00034),
            α = FT(0.933),
        )
        photosynthesis_model = PModel{FT}(photo_params);
    else
        # Set the conductance model 
        cond_params = MedlynConductanceParameters(FT; g1)
        conductance_model = MedlynConductanceModel{FT}(cond_params);

        # Set the photosynthesis model 
        photo_params = FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25)
        photosynthesis_model = FarquharModel{FT}(photo_params);
    end

    # Set the autotrophic respiration model
    AR_params = AutotrophicRespirationParameters(FT)
    autotrophic_respiration_model = AutotrophicRespirationModel{FT}(AR_params);

    # Set the plant hydraulics model
    RAI = maxLAI * f_root_to_shoot
    ai_parameterization = PrescribedSiteAreaIndex{FT}(LAI, SAI, RAI)
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a);
    plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        rooting_depth = rooting_depth,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )
    compartment_midpoints =
        n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
    compartment_surfaces =
        n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]
    plant_hydraulics_model = PlantHydraulics.PlantHydraulicsModel{FT}(;
        parameters = plant_hydraulics_ps,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_midpoints = compartment_midpoints,
        compartment_surfaces = compartment_surfaces,
    );

    # Define the soil moisture stress model 
    if soil_moisture_stress == "piecewise"
        soil_moisture_stress_params = PiecewiseMoistureStressParametersFromHydrology(
            hydrology_cm,
            ν,
            θ_r,
            zmin;
            verbose = true
        )

        soil_moisture_stress_model = PiecewiseMoistureStressModel{FT}(soil_moisture_stress_params)
    elseif soil_moisture_stress == "no_sms"
        soil_moisture_stress_model = NoMoistureStressModel{FT}()
    else
        error("Unknown soil moisture stress model: $(soil_moisture_stress)")
    end

    # instantiate the canopy model with all the components defined above.
    canopy = ClimaLand.Canopy.CanopyModel{FT}(;
        parameters = shared_params,
        domain = canopy_domain,
        autotrophic_respiration = autotrophic_respiration_model,
        radiative_transfer = radiative_transfer_model,
        photosynthesis = photosynthesis_model,
        conductance = conductance_model,
        soil_moisture_stress = soil_moisture_stress_model,
        hydraulics = plant_hydraulics_model,
        boundary_conditions = Canopy.AtmosDrivenCanopyBC(
            atmos,
            radiation,
            soil,
        ),
    );

    return canopy
end

(; lat, long, time_offset, h_canopy, atmospheric_sensor_height,
    swc_depths, ts_depths) = get_site_domain(site_ID)

(start_date_from_data, end_date_from_data) =
    FluxnetSimulationsExt.get_data_dates(site_ID, time_offset)
start_date = isnothing(start_date) ? start_date_from_data : start_date
end_date = isnothing(end_date) ? end_date_from_data : end_date
start_date_local = start_date - Dates.Hour(time_offset)
data_dt = Float64(FluxnetSimulationsExt.get_data_dt(site_ID))

canopy = setup_standalone_canopy_model(
    FT,
    site_ID,
    earth_param_set;
    use_default_param_maps = true,
    photo_model = photo_model,
    soil_moisture_stress = soil_moisture_stress,
)

set_ic! = FluxnetSimulationsExt.make_set_fluxnet_initial_conditions_standalone_canopy(
    site_ID,
    start_date,
    time_offset,
    canopy,
);

dict_writer = ClimaDiagnostics.Writers.DictWriter()

simulation = LandSimulation(
    FT,
    start_date,
    end_date,
    Δt,
    canopy;
    set_ic! = set_ic!,
    user_callbacks = (),
    diagnostics = ClimaLand.default_diagnostics(
        canopy,
        start_date_local; # note that saved times will be indexed to local time
        output_writer = dict_writer,
        output_vars = saved_var_names,
        average_period = average_period,
    ),
    driver_update_period = data_dt,
)

@time ClimaLand.Simulations.solve!(simulation)

# Save the output to NetCDF files
time_label_dict = Dict(
    :halfhourly => "30m",
    :hourly => "1h",
    :daily => "1d",
)

for short_name in saved_var_names
    diag_name = short_name * "_$(time_label_dict[average_period])_average"
    nc_file = joinpath(save_dir, "$(diag_name)_$(photo_model)_$(soil_moisture_stress).nc")
    isfile(nc_file) && rm(nc_file)

    @info "Saving diagnostic $diag_name to $nc_file"
    time_vector, data_vector = ClimaLand.Diagnostics.diagnostic_as_vectors(dict_writer, diag_name)
    
    if time_vector[1] isa ITime
        epoch = start_date_local
        # TODO: figure out why time_vector is ahead by one timestep
        time_vector = Float64.(time_vector) .- Δt 
    end

    NCDatasets.Dataset(nc_file, "c") do ds
        NCDatasets.defDim(ds, "time", length(time_vector))
        time_var = NCDatasets.defVar(ds, "time", FT, ("time",))
        time_var[:] = time_vector

        nc_var = NCDatasets.defVar(ds, short_name, FT, ("time",))
        nc_var[:] = data_vector

        time_var.attrib["units"] = "seconds since $(epoch)"

        if soil_moisture_stress == "piecewise"
            ds.attrib["field_capacity"] = parent(canopy.soil_moisture_stress.parameters.θ_c)[1]
            ds.attrib["wilting_point"] = parent(canopy.soil_moisture_stress.parameters.θ_w)[1]
        end
    end
end
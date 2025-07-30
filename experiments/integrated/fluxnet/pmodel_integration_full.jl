"""
This script runs a standalone canopy model with the P-model for photosynthesis and stomatal
conductance. We use the US-MOz (Ozark) site from the FLUXNET dataset for atmospheric forcing. 
We integrate the canopy model for one year and save the outputs to a NetCDF file. 
"""

import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Statistics
using Dates
using Insolation
using StatsBase
using NCDatasets

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaUtilities.OutputPathGenerator: generate_output_path
using ClimaDiagnostics
using ClimaUtilities
const FT = Float64
earth_param_set = LP.LandParameters(FT)
climaland_dir = pkgdir(ClimaLand)

save_outputs = true 
verbose = true

include(joinpath(climaland_dir, "experiments/integrated/fluxnet/data_tools.jl"))
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/plot_utils.jl"))

# Read in the site to be run from the command line
if length(ARGS) < 3
    error("Provide args <site_ID> <photo_model> <exp_name> in the command line")
end

site_ID = ARGS[1]
photo_model = ARGS[2] 
exp_name = ARGS[3]

save_dir = joinpath(climaland_dir, "outputs/fluxnet/$(site_ID)_$(photo_model)_$(exp_name)/")
if save_outputs && !isdir(save_dir)
    mkpath(save_dir)
end

@assert photo_model in ("pmodel", "farquhar") "Photo model must be either 'pmodel' or 'farquhar'"

# Read all site-specific domain parameters from the simulation file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_simulation.jl",
    ),
)

include(
    joinpath(climaland_dir, "experiments/integrated/fluxnet/fluxnet_domain.jl"),
)

# Read all site-specific parameters from the parameter file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_parameters.jl",
    ),
)

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/fluxnet_simulation.jl",
    ),
)

include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
)

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)
soil_ps = Soil.EnergyHydrologyParameters(
    FT;
    ν = soil_ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    K_sat = soil_K_sat,
    S_s = soil_S_s,
    θ_r,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
    albedo = soil_albedo,
);

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}
soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
soilco2_args = (; domain = soil_domain, parameters = soilco2_ps)

# Now we set up the canopy model, which we set up by component

# photosynthesis model 
if photo_model == "pmodel"
    photosynthesis_model = Canopy.PModel{FT}
    conductance_model = Canopy.PModelConductance{FT}

    conductance_args = (; parameters = PModelConductanceParameters(Drel = FT(1.6)))
    photosynthesis_args = (; parameters = PModelParameters(
        cstar = FT(0.41),
        β = FT(146),
        ϕc = FT(0.087),
        ϕ0 = FT(NaN),
        ϕa0 = FT(0.352),
        ϕa1 = FT(0.022),
        ϕa2 = FT(-0.00034),
        α = FT(0.933)
    ))
else
    photosynthesis_model = Canopy.FarquharModel{FT}
    conductance_model = Canopy.MedlynConductanceModel{FT}

    conductance_args = (; parameters = MedlynConductanceParameters(FT; g1))
    is_c3 = FT(1) 
    photosynthesis_args =
        (; parameters = FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
end

# soil moisture stress model
SM_params = PiecewiseMoistureStressParametersFromHydrology(
    soil_ps.hydrology_cm,
    soil_ps.ν,
    soil_ps.θ_r;
    c = FT(1.0),
    β0 = FT(1.0),
    verbose = true
)

soil_moisture_stress_model = PiecewiseMoistureStressModel{FT}

soil_moisture_stress_args = (; parameters = SM_params)

# Component Types
canopy_component_types = (;
    autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
    radiative_transfer = Canopy.TwoStreamModel{FT},
    photosynthesis = photosynthesis_model,
    conductance = conductance_model,
    soil_moisture_stress = soil_moisture_stress_model,
    hydraulics = Canopy.PlantHydraulicsModel{FT},
    energy = Canopy.BigLeafEnergyModel{FT},
)
# Individual Component arguments
# Set up autotrophic respiration
autotrophic_respiration_args =
    (; parameters = AutotrophicRespirationParameters(FT))
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = TwoStreamParameters(
        FT;
        Ω,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        G_Function,
    )
)

# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    rooting_depth = rooting_depth,
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
    soil_moisture_stress = soil_moisture_stress_args,
    hydraulics = plant_hydraulics_args,
    energy = energy_args,
)

# Other info needed
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = canopy_domain)

# Snow model
snow_parameters = SnowParameters{FT}(dt; earth_param_set = earth_param_set);
snow_args = (; parameters = snow_parameters, domain = canopy_domain);
snow_model_type = Snow.SnowModel
# Integrated plant hydraulics and soil model
land_input = (
    atmos = atmos,
    radiation = radiation,
    soil_organic_carbon = Csom,
    runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
)
land = LandModel{FT}(;
    soilco2_type = soilco2_type,
    soilco2_args = soilco2_args,
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
    snow_args = snow_args,
    snow_model_type = snow_model_type,
)

Y, p, cds = initialize(land)
@info "Land model initialized" 

#Initial conditions
Y.soil.ϑ_l =
    drivers.SWC.status != absent ?
    drivers.SWC.values[1 + Int(round(t0 / DATA_DT))] : soil_ν / 2 # Get soil water content at t0
# Both data and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 =
    drivers.TS.status != absent ?
    drivers.TS.values[1 + Int(round(t0 / DATA_DT))] :
    drivers.TA.values[1 + Int(round(t0 / DATA_DT))] + 40# Get soil temperature at t0
ρc_s =
    volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        land.soil.parameters.ρc_ds,
        earth_param_set,
    )
Y.soil.ρe_int =
    volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T_0, earth_param_set)
Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
ψ_stem_0 = FT(-1e5 / 9800) # pressure in the leaf divided by rho_liquid*gravitational acceleration [m]
ψ_leaf_0 = FT(-2e5 / 9800)
ψ_comps = n_stem > 0 ? [ψ_stem_0, ψ_leaf_0] : ψ_leaf_0

S_l_ini =
    inverse_water_retention_curve.(retention_model, ψ_comps, plant_ν, plant_S_s)

for i in 1:(n_stem + n_leaf)
    Y.canopy.hydraulics.ϑ_l.:($i) .=
        augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

Y.canopy.energy.T = drivers.TA.values[1 + Int(round(t0 / DATA_DT))] # Get atmos temperature at t0

Y.snow.S .= 0.0
Y.snow.S_l .= 0.0
Y.snow.U .= 0.0

set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0);
@info "Initial cache set"

exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)
jacobian! = make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);


# Callbacks
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(save_dir)

start_date = UTC_DATETIME[1] - Dates.Hour(time_offset)
dict_writer = ClimaDiagnostics.Writers.DictWriter()

short_names_1D = [
    "gs",
    "gpp",
    "ct", # canopy temperature 
    "et",
    "msf",
    "lhf",
]

if photo_model == "pmodel"
    short_names_1D = vcat(short_names_1D, ["vcmax25", "vcmax", "ac", "aj"])
end

short_names_2D = ["swc"]
output_vars = [short_names_1D..., short_names_2D...]
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = dict_writer,
    output_vars,
    average_period = :halfhourly,
)

diagnostic_handler = ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0, dt = dt);

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler);

## How often we want to update the drivers. Note that this uses the defined `t0` and `tf`
## defined in the simulatons file
updateat = Array(t0:DATA_DT:tf)
model_drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

if photo_model == "pmodel" 
    pmodel_cb = ClimaLand.make_PModel_callback(FT, UTC_DATETIME[1], dt, land.canopy, long)
    cb = SciMLBase.CallbackSet(diag_cb, driver_cb, pmodel_cb)
else
    cb = SciMLBase.CallbackSet(diag_cb, driver_cb)
end

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
);

@time sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb);

ClimaLand.Diagnostics.close_output_writers(diags)
half_hourly_diag_name_2D = short_names_2D .* "_30m_average"

for short_name in vcat(short_names_1D, short_names_2D)
    diag_name = short_name * "_30m_average"
    nc_file = joinpath(save_dir, "$(diag_name).nc")

    isfile(nc_file) && rm(nc_file) # remove existing file to overwrite

    verbose && @info "Saving $diag_name to $nc_file"
    if short_name in short_names_2D
        time_vector, data_vector = ClimaLand.Diagnostics.diagnostic_as_vectors(dict_writer, diag_name, 
            layer=1:land_domain.nelements[1])
    else
        time_vector, data_vector = ClimaLand.Diagnostics.diagnostic_as_vectors(dict_writer, diag_name)
    end 

    NCDatasets.Dataset(nc_file, "c") do ds
        NCDatasets.defDim(ds, "time", length(time_vector))
        time_var = NCDatasets.defVar(ds, "time", FT, ("time",))
        time_var[:] = time_vector
        
        if short_name in short_names_2D
            NCDatasets.defDim(ds, "z", length(data_vector[1]))
            depth_var = NCDatasets.defVar(ds, "z", FT, ("z",))
            depth_var[:] = Array(parent(land_domain.fields.z))
            nc_var = NCDatasets.defVar(ds, short_name, FT, ("time", "z"))
            nc_var[:, :] = data_vector
        else 
            nc_var = NCDatasets.defVar(ds, short_name, FT, ("time",))
            nc_var[:] = data_vector
        end
    end
end

verbose && @info "All diagnostic variables saved to NetCDF files in $save_dir"


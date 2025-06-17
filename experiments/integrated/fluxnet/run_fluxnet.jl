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

include(joinpath(climaland_dir, "experiments/integrated/fluxnet/data_tools.jl"))
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/plot_utils.jl"))

# Read in the site to be run from the command line
#if length(ARGS) < 1
#    error("Must provide site ID on command line")
#end
#
site_ID = "US-MOz"
params = Dict("g1" => 1.0, "Vcmax25" => 1.0)


function zenith_angle(
        t,
        start_date;
        latitude = lat,
        longitude = long,
        insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
    ) where {FT}
    # This should be time in UTC
    current_datetime = start_date + Dates.Second(round(t))

    # Orbital Data uses Float64, so we need to convert to our sim FT
    d, δ, η_UTC =
    FT.(
        Insolation.helper_instantaneous_zenith_angle(
                                                     current_datetime,
                                                     start_date,
                                                     insol_params,
                                                    )
       )

    FT(
       Insolation.instantaneous_zenith_angle(
                                             d,
                                             δ,
                                             η_UTC,
                                             longitude,
                                             latitude,
                                            )[1],
      )
end

# Currently, site_ID can only be "US-MOz", "US-Var", "US-Ha1", "US-NR1"
# But we want ALL fluxnet sites
# to do: change the included files to read meta data, all same config otherwise
# (same parameters, same domains, same resolution, etc.)
# data files will be read directly from caltech hpc filepath

function run_single_site(params, site_ID) # e.g., run_single_site(params, "US-MOz")

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
                                             PAR_albedo = soil_α_PAR,
                                             NIR_albedo = soil_α_NIR,
                                            );

<<<<<<< HEAD
<<<<<<< HEAD
    soil_args = (domain = soil_domain, parameters = soil_ps)
    soil_model_type = Soil.EnergyHydrology{FT}

    # Soil microbes model
    soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

    soilco2_ps = SoilCO2ModelParameters(FT)

    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

    # Set the soil CO2 BC to being atmospheric CO2
    soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
    soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (MicrobeProduction{FT}(),)

    soilco2_boundary_conditions = (; top = soilco2_top_bc, bottom = soilco2_bot_bc)
=======
=======
>>>>>>> 1333c48a9 (rebased from origin (for real this time?))
# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}
soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
soilco2_args = (; domain = soil_domain, parameters = soilco2_ps)
<<<<<<< HEAD
>>>>>>> origin

=======
=======
    soil_args = (domain = soil_domain, parameters = soil_ps)
    soil_model_type = Soil.EnergyHydrology{FT}

    # Soil microbes model
    soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

    soilco2_ps = SoilCO2ModelParameters(FT)

    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

    # Set the soil CO2 BC to being atmospheric CO2
    soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
    soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (MicrobeProduction{FT}(),)

    soilco2_boundary_conditions = (; top = soilco2_top_bc, bottom = soilco2_bot_bc)
>>>>>>> ec7314fbe (fluxnet function)

>>>>>>> 1333c48a9 (rebased from origin (for real this time?))
    soilco2_args = (;
                    boundary_conditions = soilco2_boundary_conditions,
                    sources = soilco2_sources,
                    domain = soil_domain,
                    parameters = soilco2_ps,
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
    # Set up conductance
    conductance_args = (; parameters = MedlynConductanceParameters(FT; g1 = params["g1"]))
    # Set up photosynthesis
    is_c3 = FT(1) # set the photosynthesis mechanism to C3
    photosynthesis_args =
    (; parameters = FarquharParameters(FT, is_c3; Vcmax25 = params["Vcmax25"]))
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

    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = make_imp_tendency(land)
    jacobian! = make_jacobian(land);
    jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);


    # Callbacks
    outdir = joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/out")
    output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(outdir)

    d_writer = ClimaDiagnostics.Writers.DictWriter()

    ref_time = DateTime(2010)

    diags = ClimaLand.default_diagnostics(
                                          land,
                                          ref_time;
                                          output_writer = d_writer,
                                          output_vars = :long,
                                          average_period = :hourly,
                                         )

    diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0, dt = dt);

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler);

    ## How often we want to update the drivers. Note that this uses the defined `t0` and `tf`
    ## defined in the simulatons file
    updateat = Array(t0:DATA_DT:tf)
    model_drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(model_drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    cb = SciMLBase.CallbackSet(driver_cb, diag_cb)


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

    sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb);
    return sol
end

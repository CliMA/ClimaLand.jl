import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.OutputPathGenerator: generate_output_path
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

import CairoMakie
import GeoMakie
using Statistics
using Dates
import NCDatasets

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import ClimaUtilities
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
context = ClimaComms.context()
ClimaComms.init(context)
outdir = generate_output_path("experiments/integrated/global")

device_suffix =
    typeof(context.device) <: ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

FT = Float64
earth_param_set = LP.LandParameters(FT)
prognostic_land_components = (:canopy, :soil, :soilco2)

# Set up the domain
nelements = (50, 10)
dz_tuple = (10.0, 0.1)
domain = ClimaLand.Domains.global_domain(FT; nelements, dz_tuple)
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface

# Set up dates and times for the simulation
start_date = DateTime(2008);
t0 = 0.0
dt = 450.0
tf = 3600

# Forcing data
era5_ncdata_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT;
    time_interpolation_method = time_interpolation_method,
)

soil_forcing = (; atmos, radiation)
soil = Soil.EnergyHydrology{FT}(
    domain,
    soil_forcing,
    earth_param_set;
    prognostic_land_components,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
)

# Soil microbes model
soil_organic_carbon =
    ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers = Soil.Biogeochemistry.SoilDrivers(
    co2_prognostic_soil,
    soil_organic_carbon,
    atmos,
)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(; domain, drivers)

# Now we set up the canopy model, which mostly use defaults for:
ground = ClimaLand.PrognosticSoilConditions{FT}()
canopy_domain = ClimaLand.obtain_surface_domain(domain)
canopy_forcing = (; atmos, radiation, ground)

# Set up plant hydraulics
modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
    context = nothing,
    start_date = start_date + Second(t0),
    end_date = start_date + Second(t0) + Second(tf),
)
LAI = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date,
)

canopy = Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    earth_param_set;
    prognostic_land_components,
)

# Combine the soil and canopy models into a single prognostic land model
land = SoilCanopyModel{FT}(soilco2, soil, canopy)

Y, p, cds = initialize(land)

ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context)
set_initial_state! =
    ClimaLand.Simulations.make_set_initial_state_from_file(ic_path, land)
set_initial_cache! = make_set_initial_cache(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = ClimaLand.make_imp_tendency(land)
jacobian! = ClimaLand.make_jacobian(land)

set_initial_state!(Y, p, t0, land)
set_initial_cache!(p, Y, t0)

stepper = CTS.ARS343()
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewTimeStep),
    ),
)

# set up jacobian info
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!)

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

# ClimaDiagnostics
nc_writer =
    ClimaDiagnostics.Writers.NetCDFWriter(subsurface_space, outdir; start_date)

diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = nc_writer,
    average_period = :hourly,
)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = dt)

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

updateat = Array(t0:(3600 * 3):tf)
drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

sol = ClimaComms.@time ClimaComms.device() SciMLBase.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = SciMLBase.CallbackSet(driver_cb, diag_cb),
)

ClimaLand.Diagnostics.close_output_writers(diags)

# ClimaAnalysis
if ClimaComms.iamroot(context)
    simdir = ClimaAnalysis.SimDir(outdir)

    for short_name in ClimaAnalysis.available_vars(simdir)
        var = get(simdir; short_name)
        times = var.dims["time"]
        for t in times
            fig = CairoMakie.Figure(size = (800, 600))
            kwargs = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict()
            viz.heatmap2D_on_globe!(
                fig,
                ClimaAnalysis.slice(var, time = t; kwargs...),
            )
            CairoMakie.save(joinpath(outdir, "$short_name.png"), fig)
        end
    end
end

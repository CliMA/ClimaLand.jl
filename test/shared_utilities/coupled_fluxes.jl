using ClimaComms
ClimaComms.@import_required_backends
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name
using ClimaUtilities.ClimaArtifacts
using Dates
using Test
using StaticArrays
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Soil
using ClimaLand.Canopy
using ClimaLand.Snow
import ClimaLand.Parameters as LP
using ClimaCore
using ClimaUtilities.TimeManager: ITime

FT = Float64
context = ClimaComms.context()
nelements = (101, 15)
Δt = 450.0
FT = Float64
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)

domain = ClimaLand.Domains.global_domain(
    FT;
    nelements = nelements,
    mask_threshold = FT(0.99),
);
surface_space = domain.space.surface;
start_date = DateTime(2008);
stop_date = start_date + Second(Δt)
gustiness = FT(1)
Δh = ClimaCore.Fields.zeros(surface_space) .+ 10
atmos = ClimaLand.CoupledAtmosphere{FT, typeof(Δh)}(Δh, gustiness)
latitude = ClimaCore.Fields.coordinate_field(surface_space).lat
longitude = ClimaCore.Fields.coordinate_field(surface_space).long
radiation = ClimaLand.CoupledRadiativeFluxes{FT}(
    start_date;
    latitude,
    longitude,
    toml_dict,
)
forcing = (; atmos, radiation)
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt;
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
);

Y, p, cds = ClimaLand.initialize(land)
# Set drivers (We assume the coupler will have done this)
p.drivers.T .= FT(289)
p.drivers.q .= FT(0)
p.drivers.P .= FT(101325)
p.drivers.P_liq .= FT(0)
p.drivers.P_snow .= FT(0)
p.drivers.c_co2 .= FT(4e-4)
p.drivers.u.data.:1 .= FT(2)
p.drivers.u.data.:2 .= FT(3)
p.drivers.SW_d .= FT(500)
p.drivers.LW_d .= FT(100)
p.drivers.frac_diff .= FT(0.5)
for component in (:soil, :canopy, :snow)
    cache = getproperty(p, component).turbulent_fluxes
    @test :buoyancy_flux ∈ propertynames(cache)
    @test :ρτxz ∈ propertynames(cache)
    @test :ρτyz ∈ propertynames(cache)
end

if !(ClimaComms.device() isa ClimaComms.CUDADevice)
    ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context)
    set_ic! =
        ClimaLand.Simulations.make_set_initial_state_from_file(ic_path, land)
    t0 = ITime(0, Second(1), start_date)
    set_ic!(Y, p, t0, land)
    set_initial_cache! = make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0)

    #TODO: After the LW/SW split PR, update this test below to check LW and SW are not nan
    for component in (:soil, :canopy, :snow)
        cache = getproperty(p, component).turbulent_fluxes
        if component == :soil
            for field in (:lhf, :shf, :vapor_flux_ice, :vapor_flux_liq)
                @test sum(isnan.(parent(getproperty(cache, field)))) == 0
            end
        else
            for field in (:lhf, :shf, :vapor_flux)
                @test sum(isnan.(parent(getproperty(cache, field)))) == 0
            end
        end
    end
end


# Repeat with Bucket
albedo = ClimaLand.Bucket.PrescribedBaregroundAlbedo(toml_dict, surface_space);
τc = FT(3600);
parameters = ClimaLand.Bucket.BucketModelParameters(toml_dict; albedo, τc);
radiation = ClimaLand.Bucket.CoupledRadiativeFluxes{FT}();
bucket = ClimaLand.Bucket.BucketModel(;
    parameters,
    domain,
    atmosphere = atmos,
    radiation,
);

Y, p, cds = ClimaLand.initialize(bucket)
# Set drivers (We assume the coupler will have done this)
p.drivers.T .= FT(289)
p.drivers.q .= FT(0)
p.drivers.P .= FT(101325)
p.drivers.P_liq .= FT(0)
p.drivers.P_snow .= FT(0)
p.drivers.c_co2 .= FT(4e-4)
p.drivers.u.data.:1 .= FT(2)
p.drivers.u.data.:2 .= FT(3)
p.drivers.SW_d .= FT(500)
p.drivers.LW_d .= FT(100)

t0 = ITime(0, Second(1), start_date)
Y.bucket.W .= FT(0.1)
Y.bucket.T .= FT(290)
set_initial_cache! = make_set_initial_cache(bucket)
set_initial_cache!(p, Y, t0)

cache = p.bucket.turbulent_fluxes
@test :buoyancy_flux ∈ propertynames(cache)
@test :ρτxz ∈ propertynames(cache)
@test :ρτyz ∈ propertynames(cache)
for field in (:lhf, :shf, :vapor_flux)
    @test sum(isnan.(parent(getproperty(cache, field)))) == 0
end
@test sum(isnan.(parent(p.bucket.R_n))) == 0

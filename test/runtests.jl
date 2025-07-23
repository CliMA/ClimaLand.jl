using SafeTestsets

import ClimaComms

# Performance and code quality tests
@safetestset "Aqua tests" begin
    include("aqua.jl")
end

# Shared ClimaLand utilities tests
@safetestset "Richards model implicit timestepping tests" begin
    include("shared_utilities/implicit_timestepping/richards_model.jl")
end
@safetestset "Full soil model implicit timestepping tests" begin
    include("shared_utilities/implicit_timestepping/energy_hydrology_model.jl")
end
@safetestset "Domains module tests" begin
    include("shared_utilities/domains.jl")
end
@safetestset "General utilities tests" begin
    include("shared_utilities/utilities.jl")
end
@safetestset "Variable types tests" begin
    include("shared_utilities/variable_types.jl")
end
@safetestset "Driver tests" begin
    include("shared_utilities/drivers.jl")
end

# Standalone Bucket model tests
@safetestset "Bucket albedo types tests" begin
    include("standalone/Bucket/albedo_types.jl")
end
@safetestset "Bucket snow tests" begin
    include("standalone/Bucket/snow_bucket_tests.jl")
end
@safetestset "Bucket soil tests" begin
    include("standalone/Bucket/soil_bucket_tests.jl")
end
@safetestset "Restart tests" begin
    include("standalone/Bucket/restart.jl")
end

# Standalone Snow model tests
@safetestset "Snow parameterization tests" begin
    include("standalone/Snow/parameterizations.jl")
    include("standalone/Snow/snow.jl")
end
@safetestset "Neural Snow model tools tests" begin
    include("standalone/Snow/tool_tests.jl")
end
@safetestset "Snow integrated water and energy content" begin
    include("standalone/Snow/conservation.jl")
end

@safetestset "Snow parameter constructors" begin
    include("standalone/Snow/parameters.jl")
end


# Standalone Soil model tests
@safetestset "Soil Biogeochemistry module tests" begin
    include("standalone/Soil/Biogeochemistry/biogeochemistry_module.jl")
end
@safetestset "Soil CO2 parameterization tests" begin
    include("standalone/Soil/Biogeochemistry/co2_parameterizations.jl")
end

@safetestset "Soil climate drivers tests" begin
    include("standalone/Soil/climate_drivers.jl")
end
@safetestset "Soil runoff tests" begin
    include("standalone/Soil/runoff.jl")
end
@safetestset "Soil boundary condition tests" begin
    include("standalone/Soil/soil_bc.jl")
end
@safetestset "Soil parameterization tests" begin
    include("standalone/Soil/soil_parameterizations.jl")
end
@safetestset "Soil 3D domain tests" begin
    include("standalone/Soil/soil_test_3d.jl")
    include("standalone/Soil/mask_test.jl")
end
@safetestset "Soil integration tests" begin
    include("standalone/Soil/soiltest.jl")
end
@safetestset "Soil integrated water and energy content" begin
    include("standalone/Soil/conservation.jl")
end

@safetestset "Soil spatial parameters and parameter constructors" begin
    include("standalone/Soil/parameters.jl")
end

# Standalone Surface Water model tests
@safetestset "Pond module tests" begin
    include("standalone/SurfaceWater/pond_test.jl")
end

# Standalone Vegetation model tests
@safetestset "Canopy module tests" begin
    include("standalone/Vegetation/canopy_model.jl")
end
@safetestset "Canopy PlantHydraulics tests" begin
    include("standalone/Vegetation/plant_hydraulics_test.jl")
end
@safetestset "Big Leaf model tests" begin
    include("standalone/Vegetation/test_bigleaf_parameterizations.jl")
end
@safetestset "Two Stream model tests" begin
    include("standalone/Vegetation/test_two_stream.jl")
end
@safetestset "Canopy integrated water and energy content" begin
    include("standalone/Vegetation/conservation.jl")
end
@safetestset "Canopy spatial parameters" begin
    include("standalone/Vegetation/spatial_parameters.jl")
end
@safetestset "P model tests" begin
    include("standalone/Vegetation/test_pmodel.jl")
end

# Integrated LSM tests
@safetestset "Integrated LSM unit tests" begin
    include("integrated/lsms.jl")
end
@safetestset "Integrated soil/canopy unit tests" begin
    include("integrated/soil_canopy_lsm.jl")
end
@safetestset "Integrated pond/soil LSM tests" begin
    include("integrated/pond_soil_lsm.jl")
end
@safetestset "Integrated soil energy/hydrology/biogeochem LSM tests" begin
    include("integrated/soil_energy_hydrology_biogeochemistry.jl")
end
@safetestset "Integrated soil and snow" begin
    include("integrated/soil_snow.jl")
end
@safetestset "Full land" begin
    include("integrated/full_land.jl")
end

# Diagnostics
@safetestset "Diagnostics" begin
    include("diagnostics/diagnostics_tests.jl")
end

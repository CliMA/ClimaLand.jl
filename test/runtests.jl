using SafeTestsets

import ClimaComms

gpu_broken = ClimaComms.device() isa ClimaComms.CUDADevice

# Performance and code quality tests
@safetestset "Aqua tests" begin
    include("aqua.jl")
end

# Integrated LSM tests
@safetestset "Integrated LSM unit tests" begin
    include("integrated/lsms.jl")
end
@safetestset "Integrated LSM unit tests" begin
    include("integrated/soil_canopy_lsm.jl")
end
@safetestset "Integrated pond/soil LSM tests" begin
    include("integrated/pond_soil_lsm.jl")
end
gpu_broken ||
    @safetestset "Integrated soil energy/hydrology/biogeochem LSM tests" begin
        include("integrated/soil_energy_hydrology_biogeochemistry.jl")
    end

# Shared ClimaLSM utilities tests
@safetestset "Richards model implicit timestepping tests" begin
    include("shared_utilities/implicit_timestepping/richards_model.jl")
end
@safetestset "Domains module tests" begin
    include("shared_utilities/domains.jl")
end
gpu_broken || @safetestset "FileReader module tests" begin
    include("shared_utilities/file_reader.jl")
end
gpu_broken || @safetestset "Regridder module tests" begin
    include("shared_utilities/regridder.jl")
end
gpu_broken || @safetestset "General utilities tests" begin
    include("shared_utilities/utilities.jl")
end
@safetestset "Variable types tests" begin
    include("shared_utilities/variable_types.jl")
end

# Standalone Bucket model tests
gpu_broken || @safetestset "Bucket albedo types tests" begin
    include("standalone/Bucket/albedo_types.jl")
end
gpu_broken || @safetestset "Bucket snow tests" begin
    include("standalone/Bucket/snow_bucket_tests.jl")
end
gpu_broken || @safetestset "Bucket soil tests" begin
    include("standalone/Bucket/soil_bucket_tests.jl")
end

# Standalone Snow model tests
@safetestset "Snow parameterization tests" begin
    include("standalone/Snow/parameterizations.jl")
end

# Standalone Soil model tests
gpu_broken || @safetestset "Soil Biogeochemistry module tests" begin
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
end
gpu_broken || @safetestset "Soil integration tests" begin
    include("standalone/Soil/soiltest.jl")
end

# Standalone Surface Water model tests
@safetestset "Pond module tests" begin
    include("standalone/SurfaceWater/pond_test.jl")
end

# Standalone Vegetation model tests
gpu_broken || @safetestset "Canopy module tests" begin
    include("standalone/Vegetation/canopy_model.jl")
end
gpu_broken || @safetestset "Canopy PlantHydraulics tests" begin
    include("standalone/Vegetation/plant_hydraulics_test.jl")
end
@safetestset "Big Leaf model tests" begin
    include("standalone/Vegetation/test_bigleaf_parameterizations.jl")
end
@safetestset "Two Stream model tests" begin
    include("standalone/Vegetation/test_two_stream.jl")
end

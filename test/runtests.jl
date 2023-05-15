using SafeTestsets
#! format: off
@safetestset "parameterizations" begin @time include("./Snow/parameterizations.jl") end
@safetestset "richards_model" begin @time include("./implicit_timestepping/richards_model.jl") end
@safetestset "test_bigleaf_parameterizations" begin @time include("./Vegetation/test_bigleaf_parameterizations.jl") end
@safetestset "canopy_model" begin @time include("./Vegetation/canopy_model.jl") end
@safetestset "co2_parameterizations" begin @time include("./Soil/Biogeochemistry/co2_parameterizations.jl") end
@safetestset "biogeochemistry_module" begin @time include("./Soil/Biogeochemistry/biogeochemistry_module.jl") end
@safetestset "soil_bucket_tests" begin @time include("./Bucket/soil_bucket_tests.jl") end
@safetestset "snow_bucket_tests" begin @time include("./Bucket/snow_bucket_tests.jl") end
@safetestset "regridder" begin @time include("./Bucket/regridder.jl") end
@safetestset "albedo_map_bucket" begin @time include("./Bucket/albedo_map_bucket.jl") end
@safetestset "plant_hydraulics_test" begin @time include("./Vegetation/plant_hydraulics_test.jl") end
@safetestset "domains" begin @time include("./domains.jl") end
@safetestset "variable_types" begin @time include("./variable_types.jl") end
@safetestset "pond_test" begin @time include("./SurfaceWater/pond_test.jl") end
@safetestset "pond_soil_lsm" begin @time include("./LSM/pond_soil_lsm.jl") end
@safetestset "soil_energy_hydrology_biogeochemistry" begin @time include("./LSM/soil_energy_hydrology_biogeochemistry.jl") end
@safetestset "climate_drivers" begin @time include("./Soil/climate_drivers.jl") end
@safetestset "soil_test_3d" begin @time include("./Soil/soil_test_3d.jl") end
@safetestset "soiltest" begin @time include("./Soil/soiltest.jl") end
@safetestset "soil_parameterizations" begin @time include("./Soil/soil_parameterizations.jl") end
@safetestset "soil_bc" begin @time include("./Soil/soil_bc.jl") end
#! format: on

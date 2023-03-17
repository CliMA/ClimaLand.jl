#!/bin/bash

julia --project src/Soil/examples/imex_testing.jl exp mod_picard iters_1 dt_1
julia --project src/Soil/examples/imex_testing.jl exp mod_picard iters_1 dt_4
julia --project src/Soil/examples/imex_testing.jl exp mod_picard iters_1 dt_16
julia --project src/Soil/examples/imex_testing.jl exp mod_picard iters_1 dt_64
julia --project src/Soil/examples/imex_testing.jl exp mod_picard iters_1 dt_128
julia --project src/Soil/examples/imex_testing.jl exp mod_picard iters_1 dt_256

julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_1
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_4
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_16
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_64
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_128
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_256

julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_2 dt_1
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_2 dt_4
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_2 dt_16
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_2 dt_64
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_2 dt_128
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_2 dt_256

julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_10 dt_1
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_10 dt_4
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_10 dt_16
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_10 dt_64
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_10 dt_128
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_10 dt_256

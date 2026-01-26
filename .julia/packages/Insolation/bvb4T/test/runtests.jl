using Test

push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Dates
using Statistics
using Roots
using Optim

using Insolation
import Insolation.Parameters as IP
import Insolation.OrbitalData
import ClimaParams as CP

FT = Float32
param_set = IP.InsolationParameters(FT)

@testset "Orbital Params" begin
    include("test_orbit_param.jl")
end
@testset "Types" begin
    include("test_types.jl")
end
@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end
@testset "Insolation" begin
    include("test_insolation.jl")
end
@testset "Equinox Date" begin
    include("test_perihelion.jl")
    include("test_equinox.jl")
end

FT = Float64
param_set = IP.InsolationParameters(FT)

@testset "Orbital Params" begin
    include("test_orbit_param.jl")
end
@testset "Types" begin
    include("test_types.jl")
end
@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end
@testset "Insolation" begin
    include("test_insolation.jl")
end
@testset "Equinox Date" begin
    include("test_perihelion.jl")
    include("test_equinox.jl")
end

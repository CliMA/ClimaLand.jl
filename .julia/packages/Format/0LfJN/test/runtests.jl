using Format
using Test
using Printf
using Random

@testset "cformat"    begin include( "cformat.jl" ) end
@testset "printf"     begin include( "printf.jl" ) end
@testset "fmtspec"    begin include( "fmtspec.jl" ) end
@testset "formatexpr" begin include( "formatexpr.jl" ) end
@testset "fmt"        begin include( "fmt.jl" ) end

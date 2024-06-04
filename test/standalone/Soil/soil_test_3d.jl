using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using Statistics
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: HybridBox, SphericalShell
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP

FT = Float64
@show x = 10 ^ 11 * eps(FT)
# x - 2eps() ≤ x ≤ x + 2eps()
@test 2.220446049205904e-5 ≤ x ≤ 2.220446049294722e-5

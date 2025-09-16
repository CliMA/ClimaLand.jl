# test/standalone/Vegetation/test_owus.jl
using Test
using ClimaLand                     # loads the package (and Canopy)
using ClimaLand.Canopy.OWUSStomata  # <-- bring the submodule into scope

@testset "OWUS smoke" begin
    owus = OWUSStomatalModel(fww=0.6, s_star=0.4, s_w=0.15)

    s_vals = 0.0:0.05:1.0
    E0     = 3.5e-3 / 86400  # m/s
    VPD    = 1500.0          # Pa
    P_air  = 101325.0        # Pa
    T_air  = 298.15          # K

    for s in s_vals
        gref = Ref(0.0)
        stomatal_conductance!(gref, owus; s=s, E0=E0, VPD=VPD, P_air=P_air, T_air=T_air)
        @test gref[] ≥ 0
    end
end
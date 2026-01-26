using Test

import JuliaFormatter
import ClimaDiagnostics

@testset "Formatting" begin
    @test JuliaFormatter.format(
        ClimaDiagnostics;
        verbose = false,
        overwrite = false,
    )
end

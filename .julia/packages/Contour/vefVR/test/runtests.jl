using Contour, Test

@test length(detect_ambiguities(Base, Contour)) == 0

include("verify_vertices.jl")
include("interface.jl")


#issue 59
@inferred collect(())

using Aqua
# Aqua tests
# Intervals brings a bunch of ambiquities unfortunately
Aqua.test_all(Contour)


@static if Base.VERSION >= v"1.7"

    @info "Running JET..."

    using JET
    display(JET.report_package(Contour))
end
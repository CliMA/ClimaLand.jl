using Test
using ClimaCore: Spaces, Geometry, Fields
using ClimaLSM
using ClimaLSM: Domains, condition, SavingAffect, saving_initialize

TestFloatTypes = (Float32, Float64)

mutable struct Integrator{FT}
    t::FT
    p::Array{FT}
end

function upd_integrator(integrator::Integrator{FT}, dt::FT) where {FT}
    integrator.t += dt
    integrator.p .+= FT(1)
end

@testset "callback cond, affect! test" begin
    for FT in TestFloatTypes
        t0 = FT(0)
        dt = FT(10)
        tf = FT(100)
        t_range = collect(t0:dt:tf)

        # specify when to save in the callback (test different dts)
        saveats = (
            collect(t0:dt:tf),
            collect(t0:(2 * dt):tf),
            collect(t0:(tf - t0):tf),
            collect(t0:(0.5 * dt):tf),
        )
        for saveat in saveats
            # note: delete non-multiples of dt in saveat before calling callback
            deleteat!(saveat, findall(x -> (x - t0) % dt != 0, saveat))

            saved_values = (;
                t = Array{FT}(undef, length(saveat)),
                saveval = Array{Array{FT}}(undef, length(saveat)),
            )

            # set up components of callback
            cond = condition(saveat)
            saveiter = 0
            affect! = SavingAffect(saved_values, saveat, saveiter)

            p_init = copy(t_range)
            integrator = Integrator{FT}(t0, p_init)

            # use this to manually save `integrator` at various `t`
            integ_saved = (;
                t = Array{FT}(undef, length(saveat)),
                p = Array{Array{FT}}(undef, length(saveat)),
            )

            # simulate callback behavior
            saveiter = 0
            for t in t_range
                if cond(0, t, 0)
                    affect!(integrator)

                    # save the current state of `integrator`
                    saveiter += 1
                    integ_saved.t[saveiter] = integrator.t
                    integ_saved.p[saveiter] = copy(integrator.p)
                end
                upd_integrator(integrator, dt)
            end

            @test affect!.saved_values.t == integ_saved.t
            @test affect!.saved_values.saveval == integ_saved.p
        end
    end
end

@testset "callback saving_initialize test" begin
    for FT in TestFloatTypes
        t0 = FT(0)
        dt = FT(10)
        tf = FT(100)

        # we're only testing the save at t0 so one saveat range is sufficient
        saveat = collect(t0:dt:tf)
        saved_values = (;
            t = Array{FT}(undef, length(saveat)),
            saveval = Array{Array{FT}}(undef, length(saveat)),
        )

        saveiter = 0
        affect! = SavingAffect(saved_values, saveat, saveiter)
        cb = (; affect! = affect!)

        p_init = collect(t0:dt:tf)
        integrator = Integrator{FT}(t0, p_init)

        # case 1: t not in saveat
        t1 = FT(-1)
        @test saving_initialize(cb, 0, t1, integrator) == false
        @test cb.affect!.saved_values.t != integrator.t
        @test cb.affect!.saved_values.saveval != integrator.p[1]

        # case 2: t in saveat
        t2 = t0
        saving_initialize(cb, 0, t2, integrator)
        @test cb.affect!.saved_values.t[1] == integrator.t[1]
        @test cb.affect!.saved_values.saveval[1] == integrator.p
    end
end

@testset "dss! - 1D surface spaces" begin
    for FT in TestFloatTypes
        # Test for Spaces.PointSpace and Spaces.FiniteDifferenceSpace
        coord = Geometry.ZPoint(FT(1))
        space1 = Spaces.PointSpace(coord)

        model_domain = Domains.Column(; zlim = (FT(-1), FT(0)), nelements = 10)
        space2 = model_domain.space
        spaces = (space1, space2)

        for space in spaces
            field = Fields.coordinate_field(space)
            subfield1 = similar(field)
            parent(subfield1) .= parent(copy(field)) .+ FT(1)

            Y = Fields.FieldVector(
                field = field,
                subfields = (subfield1 = subfield1),
            )
            Y_copy = copy(Y)
            p = (;)
            ClimaLSM.dss!(Y, p, FT(0))

            # On a 1D space, we expect dss! to do nothing
            @test Y == Y_copy
        end
    end
end

@testset "dss! - 2D surface spaces" begin
    for FT in TestFloatTypes
        # Test for Spaces.SpectralElementSpace2D and Spaces.ExtrudedFiniteDifferenceSpace
        domain1 = Domains.Plane(;
            xlim = FT.((0.0, 1.0)),
            ylim = FT.((0.0, 1.0)),
            nelements = (2, 2),
            periodic = (true, true),
            npolynomial = 1,
        )
        domain2 = Domains.HybridBox(;
            xlim = FT.((0.0, 1.0)),
            ylim = FT.((0.0, 1.0)),
            zlim = FT.((-1.0, 0.0)),
            nelements = (5, 2, 5),
            npolynomial = 3,
        )

        spaces = (domain1.space, domain2.space.horizontal_space)

        for space in spaces
            field = Fields.zeros(space)
            subfield1 = Fields.zeros(space)
            subfield2 = Fields.zeros(space)

            # make fields non-constant
            for i in eachindex(parent(field))
                parent(field)[i] = sin(i)
                parent(subfield1)[i] = cos(i)
                parent(subfield2)[i] = sin(i) * cos(i)
            end

            Y = Fields.FieldVector(
                field = field,
                subfields = (subfield1 = subfield1, subfield2 = subfield2),
            )
            Y_copy = copy(Y)
            p = (;)
            ClimaLSM.dss!(Y, p, FT(0))

            # On a 2D space, we expect dss! to change Y
            @test Y.field != Y_copy.field
            @test Y.subfields.subfield1 != Y_copy.subfields.subfield1
            @test Y.subfields.subfield2 != Y_copy.subfields.subfield2
        end
    end
end

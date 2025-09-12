using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore: Spaces, Geometry, Fields
using ClimaLand
using ClimaLand: Domains, SavingAffect, IntervalBasedCallback
using Dates
import ClimaUtilities.TimeManager: ITime, date
## Callback tests
mutable struct Integrator{FT}
    t::Any
    p::Array{FT}
end

function upd_integrator(integrator::Integrator{FT}, dt) where {FT}
    integrator.t += dt
    integrator.p .+= FT(1)
end

FT = Float32

@testset "IntervalBasedCallback" begin
    # with epoch test
    start_date = ITime(0, epoch = DateTime(2010))
    dt = ITime(3600)

    increment_p = (x) -> x.p .+= 1
    # p is only incremented when the cb is triggered
    function step(integrator, cb)
        integrator.t += dt
        cb.condition(nothing, integrator.t, integrator) &&
            cb.affect!(integrator)
        return
    end

    # test update every two hours (every two dt)
    integrator = Integrator(start_date, [0])
    cb = IntervalBasedCallback(Hour(2), start_date, dt, increment_p)
    step(integrator, cb)
    @test integrator.p == [0]
    step(integrator, cb)
    @test integrator.p == [1]
    for i in 1:38
        step(integrator, cb)
    end
    @test integrator.p == [20]

    # test monthly
    integrator = Integrator(start_date, [0])
    cb = IntervalBasedCallback(Month(1), start_date, dt, increment_p)
    for i in 1:(24 * 31)
        step(integrator, cb)
    end
    @test integrator.p == [1]

    # no epoch test
    t0 = ITime(0)
    integrator = Integrator(t0, [0])
    t0, dt = promote(t0, dt)
    cb = IntervalBasedCallback(dt * 2, t0, dt, increment_p)
    step(integrator, cb)
    @test integrator.p == [0]
    step(integrator, cb)
    @test integrator.p == [1]
    for i in 1:38
        step(integrator, cb)
    end
    @test integrator.p == [20]
end



@testset "NonInterpSavingCallback - DateTimes" begin
    start_date = DateTime(2020)
    stop_date = start_date + Hour(4)
    dt = Hour(1)

    saving_cb_every_hour =
        ClimaLand.NonInterpSavingCallback(start_date, stop_date, Hour(1))

    saving_cb_every_two_hour =
        ClimaLand.NonInterpSavingCallback(start_date, stop_date, Hour(2))

    saving_cb_late_start = ClimaLand.NonInterpSavingCallback(
        start_date,
        stop_date,
        Hour(1);
        first_save_date = start_date + Hour(2),
    )


    all_cbs =
        (saving_cb_every_hour, saving_cb_every_two_hour, saving_cb_late_start)
    for t in start_date:dt:stop_date
        # our simulations use ITime internally, so check evaluating at ITimes
        t_itime = ITime(Second(t - start_date).value, epoch = start_date)
        integrator = Integrator(t_itime, [Hour(t)])
        for cb in all_cbs
            if t == start_date
                cb.initialize(nothing, nothing, nothing, integrator)
            else
                cb.condition(nothing, t, integrator) && cb.affect!(integrator)
            end
        end
    end
    @test saving_cb_every_hour.affect!.saved_values.t ==
          collect(start_date:dt:stop_date)
    @test saving_cb_every_two_hour.affect!.saved_values.t ==
          collect(start_date:(2 * dt):stop_date)
    @test saving_cb_late_start.affect!.saved_values.t ==
          collect((start_date + Hour(2)):dt:stop_date)
    # test saved values are as expected
    @test all(
        cb -> all(
            first.(cb.affect!.saved_values.saveval) .==
            Hour.(cb.affect!.saved_values.t),
        ),
        all_cbs,
    )
end

@testset "NonInterpSavingCallback - Floats" begin
    t0 = 0.0
    tf = 4.0
    dt = 1.0

    save_every_dt = ClimaLand.NonInterpSavingCallback(t0, tf, dt)
    save_every_other_dt = ClimaLand.NonInterpSavingCallback(t0, tf, 2 * dt)
    save_every_dt_late_start =
        ClimaLand.NonInterpSavingCallback(t0, tf, dt; t_first_save = 2.0)

    all_cbs = (save_every_dt, save_every_dt_late_start, save_every_other_dt)
    for t in t0:dt:tf
        t_itime = ITime(t)
        integrator = Integrator(t_itime, [t])
        for cb in all_cbs
            if t == t0
                cb.initialize(nothing, nothing, nothing, integrator)
            else
                cb.condition(nothing, t, integrator) && cb.affect!(integrator)
            end
        end
    end

    @test save_every_dt.affect!.saved_values.t == collect(t0:dt:tf)
    @test save_every_other_dt.affect!.saved_values.t == collect(t0:(2 * dt):tf)
    @test save_every_dt_late_start.affect!.saved_values.t ==
          collect((t0 + 2 * dt):dt:tf)
    # test saved values are as expected
    @test all(
        cb -> all(
            first.(cb.affect!.saved_values.saveval) .==
            float.(cb.affect!.saved_values.t),
        ),
        all_cbs,
    )
end


## DSS 0D/1D tests
@testset "dss! - 0D/1D spaces, FT = $FT" begin
    # Test for Spaces.PointSpace and Spaces.FiniteDifferenceSpace
    domain1 = ClimaLand.Domains.Point(; z_sfc = FT(0))
    domain2 =
        ClimaLand.Domains.Column(; zlim = FT.((-1.0, 0.0)), nelements = (5))
    domains = (domain1, domain2, domain2)
    spaces =
        (domain1.space.surface, domain2.space.subsurface, domain2.space.surface)
    for (domain, space) in zip(domains, spaces)
        field = Fields.coordinate_field(space)
        subfield1 = similar(field)
        parent(subfield1) .= parent(copy(field)) .+ FT(1)

        Y = Fields.FieldVector(
            field = field,
            subfields = (subfield1 = subfield1),
        )
        Y_copy = copy(Y)
        p = (;)
        ClimaLand.dss!(Y, p, FT(0))
        @test p == (;)
        @test ClimaLand.add_dss_buffer_to_aux(p, domain) == p
        # On a 1D space, we expect dss! to do nothing
        @test Y.field == Y_copy.field
        @test Y.subfields == Y_copy.subfields
    end
end

@testset "dss! - 2D spaces, FT = $FT" begin
    # Test for Spaces.SpectralElementSpace2D
    domain1 = ClimaLand.Domains.Plane(;
        xlim = FT.((0.0, 1.0)),
        ylim = FT.((0.0, 1.0)),
        nelements = (2, 2),
        periodic = (true, true),
        npolynomial = 1,
    )
    domain2 = ClimaLand.Domains.SphericalSurface(;
        radius = FT(2),
        nelements = 10,
        npolynomial = 1,
    )

    domains = (domain1, domain2)
    for domain in domains
        space = domain.space.surface
        field = Fields.zeros(space)
        subfield1 = Fields.zeros(space)
        subfield2 = Fields.zeros(space)

        # make fields non-constant in order to check the
        # impact of the dss step
        ArrayType = ClimaComms.array_type(ClimaComms.device())
        local_size = size(parent(field))

        # We prepare a new Array or a CuArray (depending on the ClimaComms device) and swap
        # in-place the parent of our fields with the new data
        sin_i = ArrayType(
            reshape([sin(i) for i in eachindex(parent(field))], local_size),
        )
        cos_i = ArrayType(
            reshape([cos(i) for i in eachindex(parent(field))], local_size),
        )
        sin_cos_i = sin_i .* cos_i

        parent(field) .= sin_i
        parent(subfield1) .= cos_i
        parent(subfield2) .= sin_cos_i

        Y = Fields.FieldVector(
            field = field,
            subfields = (subfield1 = subfield1, subfield2 = subfield2),
        )
        Y_copy = copy(Y)
        p = (;)
        p = ClimaLand.add_dss_buffer_to_aux(p, domain)
        @test typeof(p.dss_buffer_2d) ==
              typeof(Spaces.create_dss_buffer(Fields.zeros(space)))
        ClimaLand.dss!(Y, p, FT(0))

        # On a 2D space, we expect dss! to change Y
        # but we changed dss! to do nothing, so replace != with  ==
        @test Y.field == Y_copy.field
        @test Y.subfields.subfield1 == Y_copy.subfields.subfield1
        @test Y.subfields.subfield2 == Y_copy.subfields.subfield2
    end
end

@testset "dss! - 3D spaces, FT = $FT" begin
    domain1 = ClimaLand.Domains.HybridBox(;
        xlim = FT.((0.0, 1.0)),
        ylim = FT.((0.0, 1.0)),
        zlim = FT.((0.0, 1.0)),
        nelements = (2, 2, 2),
        periodic = (true, true),
        npolynomial = 1,
    )
    domain2 = ClimaLand.Domains.SphericalShell(;
        radius = FT(2),
        depth = FT(1.0),
        nelements = (10, 5),
        npolynomial = 1,
    )

    domains = (domain1, domain2)
    for domain in domains
        space = domain.space.subsurface
        field = Fields.zeros(space)
        subfield1 = Fields.zeros(space)
        subfield2 = Fields.zeros(space)

        # make fields non-constant in order to check
        # the impact of the dss step
        ArrayType = ClimaComms.array_type(ClimaComms.device())
        local_size = size(parent(field))

        sin_i = ArrayType(
            reshape([sin(i) for i in eachindex(parent(field))], local_size),
        )
        cos_i = ArrayType(
            reshape([cos(i) for i in eachindex(parent(field))], local_size),
        )
        sin_cos_i = sin_i .* cos_i

        parent(field) .= sin_i
        parent(subfield1) .= cos_i
        parent(subfield2) .= sin_cos_i

        Y = Fields.FieldVector(
            field = field,
            subfields = (subfield1 = subfield1, subfield2 = subfield2),
        )
        Y_copy = copy(Y)
        p = (;)
        p = ClimaLand.add_dss_buffer_to_aux(p, domain)
        @test typeof(p.dss_buffer_3d) ==
              typeof(Spaces.create_dss_buffer(Fields.zeros(space)))
        @test typeof(p.dss_buffer_2d) == typeof(
            Spaces.create_dss_buffer(Fields.zeros(domain.space.surface)),
        )
        ClimaLand.dss!(Y, p, FT(0))

        # On a 3D space, we expect dss! to change Y
        # but we changed dss! to do nothing, so replace != with  ==
        @test Y.field == Y_copy.field
        @test Y.subfields.subfield1 == Y_copy.subfields.subfield1
        @test Y.subfields.subfield2 == Y_copy.subfields.subfield2
    end
end

@testset "dss! - npolynomial = 0, FT = $FT" begin
    # Test for Spaces.SpectralElementSpace2D
    domain1 = ClimaLand.Domains.Plane(;
        xlim = FT.((0.0, 1.0)),
        ylim = FT.((0.0, 1.0)),
        nelements = (2, 2),
        periodic = (true, true),
    )
    domain2 =
        ClimaLand.Domains.SphericalSurface(; radius = FT(2), nelements = 10)
    domain3 = ClimaLand.Domains.HybridBox(;
        xlim = FT.((0.0, 1.0)),
        ylim = FT.((0.0, 1.0)),
        zlim = FT.((0.0, 1.0)),
        nelements = (2, 2, 2),
        periodic = (true, true),
    )
    domain4 = ClimaLand.Domains.SphericalShell(;
        radius = FT(2),
        depth = FT(1.0),
        nelements = (10, 5),
    )

    domains = (domain1, domain2, domain3, domain4)
    for domain in domains
        p = (;)
        p = ClimaLand.add_dss_buffer_to_aux(p, domain)
        @test !haskey(p, :dss_buffer_2d)
        @test !haskey(p, :dss_buffer_3d)
    end
end

@testset "count_nans_state, FT = $FT" begin
    # Test on a 3D spherical domain
    domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(2),
        depth = FT(1.0),
        nelements = (10, 5),
    )

    # Construct some fields
    space = domain.space.subsurface
    var1 = Fields.zeros(space)
    var2 = Fields.zeros(space)
    var3 = Fields.zeros(space)
    fieldvec = Fields.FieldVector(var2 = var2, var3 = var3)

    # Construct a FieldVector containing the fields and a nested FieldVector
    Y = Fields.FieldVector(var1 = var1, fieldvec = fieldvec)

    # Count and log the number of NaNs in the state (test verbose and non-verbose cases)
    @test_logs (:info, "Checking NaNs in var1") (:info, "No NaNs found") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:info, "No NaNs found") (
        :info,
        "Checking NaNs in var3",
    ) (:info, "No NaNs found") ClimaLand.call_count_nans_state(
        Y,
        verbose = true,
    )

    @test_logs (:info, "Checking NaNs in var1") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:info, "Checking NaNs in var3") ClimaLand.call_count_nans_state(
        Y,
    )

    # Add some NaNs to the fields.
    # Note that we only check if NaNs are present in the surface field, since
    # any NaN in the column will propagate to the surface within several timesteps.
    # Note: this code uses `parent` and scalar indexing,
    # which shouldn't be replicated outside of tests
    ClimaComms.allowscalar(ClimaComms.device()) do
        parent(var1)[:, 1, 1, 1, 1] .= NaN
        parent(var2)[:, 1, 1, 1, 1] .= NaN
        parent(var2)[:, 1, 1, 1, 2] .= NaN
    end

    # Count and log the number of NaNs in the state (test verbose and non-verbose cases)
    @test_logs (:info, "Checking NaNs in var1") (:warn, "1 NaNs found") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:warn, "2 NaNs found") (
        :info,
        "Checking NaNs in var3",
    ) (:info, "No NaNs found") ClimaLand.call_count_nans_state(
        Y,
        verbose = true,
    )

    @test_logs (:info, "Checking NaNs in var1") (:warn, "1 NaNs found") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:warn, "2 NaNs found") (
        :info,
        "Checking NaNs in var3",
    ) ClimaLand.call_count_nans_state(Y)

    # Test with a mask
    sfc_space = domain.space.surface

    mask_zeros = Fields.zeros(sfc_space)
    @test_logs (:info, "Checking NaNs in var1") (:info, "No NaNs found") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:info, "No NaNs found") (
        :info,
        "Checking NaNs in var3",
    ) (:info, "No NaNs found") ClimaLand.call_count_nans_state(
        Y,
        mask = mask_zeros,
        verbose = true,
    )
    @test_logs (:info, "Checking NaNs in var1") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:info, "Checking NaNs in var3") ClimaLand.call_count_nans_state(
        Y,
        mask = mask_zeros,
    )

    mask_ones = Fields.ones(sfc_space)
    @test_logs (:info, "Checking NaNs in var1") (:warn, "1 NaNs found") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:warn, "2 NaNs found") (
        :info,
        "Checking NaNs in var3",
    ) (:info, "No NaNs found") ClimaLand.call_count_nans_state(
        Y,
        mask = mask_ones,
        verbose = true,
    )

    @test_logs (:info, "Checking NaNs in var1") (:warn, "1 NaNs found") (
        :info,
        "Checking NaNs in fieldvec",
    ) (:info, "Checking NaNs in var2") (:warn, "2 NaNs found") (
        :info,
        "Checking NaNs in var3",
    ) ClimaLand.call_count_nans_state(Y, mask = mask_ones)
end

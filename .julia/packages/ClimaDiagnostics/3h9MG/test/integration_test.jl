using Test
using Profile
using ProfileCanvas

import SciMLBase
import NCDatasets

import ClimaDiagnostics

import Dates

import ClimaComms
@static if pkgversion(ClimaComms) >= v"0.6"
    ClimaComms.@import_required_backends
end

import LazyBroadcast: lazy

const context = ClimaComms.context()
ClimaComms.init(context)

include("TestTools.jl")

"""
Set up a full test problem

Increasing `more_compute_diagnostics` adds more copies of a compute diagnostic with no output.
Useful to stress allocations.
"""
function setup_integrator(
    output_dir;
    space = SphericalShellSpace(; context),
    context,
    more_compute_diagnostics = 0,
    dict_writer = nothing,
)
    t0 = 0.0
    tf = 10.0
    dt = 1.0
    args, kwargs = create_problem(space; t0, tf, dt)

    @info "Writing output to $output_dir"

    dummy_writer = ClimaDiagnostics.Writers.DummyWriter()
    h5_writer = ClimaDiagnostics.Writers.HDF5Writer(output_dir)
    if space isa ClimaCore.Spaces.PointSpace
        nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
            space,
            output_dir;
            start_date = Dates.DateTime(2015, 2, 2),
        )
    else
        nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
            space,
            output_dir;
            num_points = (10, 5, 3),
            start_date = Dates.DateTime(2015, 2, 2),
        )
    end

    function compute_my_var!(out, u, p, t)
        if isnothing(out)
            return u.my_var
        else
            out .= u.my_var
            return nothing
        end
    end

    function compute_my_var_lazy(u, p, t)
        return lazy.(u.my_var .+ 10.0)
    end

    function compute_my_var_field(u, p, t)
        return u.my_var
    end

    simple_var = ClimaDiagnostics.DiagnosticVariable(;
        compute! = compute_my_var!,
        short_name = "YO",
        long_name = "YO YO",
    )

    simple_var_lazy = ClimaDiagnostics.DiagnosticVariable(;
        compute = compute_my_var_lazy,
        short_name = "YO LAZY",
        long_name = "YO YO LAZY",
    )

    simple_var_field = ClimaDiagnostics.DiagnosticVariable(;
        compute = compute_my_var_field,
        short_name = "YO LAZY FIELD",
        long_name = "YO YO LAZY FIELD",
    )

    average_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = nc_writer,
        reduction_time_func = (+),
        output_schedule_func = ClimaDiagnostics.Schedules.DivisorSchedule(2),
        pre_output_hook! = ClimaDiagnostics.average_pre_output_hook!,
    )
    inst_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = nc_writer,
    )
    inst_diagnostic_lazy = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var_lazy,
        output_writer = nc_writer,
    )
    inst_diagnostic_field = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var_field,
        output_writer = nc_writer,
    )
    average_diagnostic_lazy = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var_lazy,
        output_writer = nc_writer,
        reduction_time_func = (+),
        output_schedule_func = ClimaDiagnostics.Schedules.DivisorSchedule(2),
        pre_output_hook! = ClimaDiagnostics.average_pre_output_hook!,
    )
    average_diagnostic_field = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var_field,
        output_writer = nc_writer,
        reduction_time_func = (+),
        output_schedule_func = ClimaDiagnostics.Schedules.DivisorSchedule(2),
        pre_output_hook! = ClimaDiagnostics.average_pre_output_hook!,
    )
    inst_every3s_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = nc_writer,
        output_schedule_func = ClimaDiagnostics.Schedules.EveryDtSchedule(
            3.0,
            t_last = t0,
        ),
    )
    inst_diagnostic_h5 = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = h5_writer,
    )

    inst_every3s_diagnostic_another = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = nc_writer,
        output_schedule_func = ClimaDiagnostics.Schedules.EveryDtSchedule(
            3.0,
            t_last = t0,
        ),
    )
    scheduled_diagnostics = [
        average_diagnostic,
        inst_diagnostic,
        inst_diagnostic_lazy,
        inst_diagnostic_field,
        average_diagnostic_lazy,
        average_diagnostic_field,
        inst_diagnostic_h5,
        inst_every3s_diagnostic,
        inst_every3s_diagnostic_another,
    ]
    if !isnothing(dict_writer)
        inst_diagnostic_dict = ClimaDiagnostics.ScheduledDiagnostic(
            variable = simple_var,
            output_writer = dict_writer,
        )
        scheduled_diagnostics = [scheduled_diagnostics..., inst_diagnostic_dict]
    end

    @test inst_every3s_diagnostic_another == inst_every3s_diagnostic
    @test !(inst_every3s_diagnostic_another === inst_every3s_diagnostic)

    # Add more weight, useful for stressing allocations
    compute_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dummy_writer,
    )
    scheduled_diagnostics = [
        scheduled_diagnostics...,
        [compute_diagnostic for _ in 1:more_compute_diagnostics]...,
    ]

    return ClimaDiagnostics.IntegratorWithDiagnostics(
        SciMLBase.init(args...; kwargs...),
        scheduled_diagnostics,
    )
end

sphere_space = SphericalShellSpace(; context)
purely_horizontal_space = ClimaCore.Spaces.level(sphere_space, 1)
# list of tuples of (space, space_name, dimensions of written diagnostics without time)
spaces_test_list = [
    (sphere_space, "SphericalShellSpace", (10, 5, 10)),
    (purely_horizontal_space, "purely horizontal space", (10, 5)),
    (
        ClimaCore.Spaces.PointSpace(context, ClimaCore.Geometry.ZPoint(1.0)),
        "PointSpace",
        (),
    ),
]
if context isa ClimaComms.SingletonCommsContext
    spaces_test_list = [
        spaces_test_list...,
        (
            ColumnCenterFiniteDifferenceSpace(10, context),
            "purely vertical space",
            (10,),
        ),
    ]
end
for (space, space_name, written_space_dims) in spaces_test_list
    @testset "A full problem using a $space_name" begin
        mktempdir() do output_dir
            output_dir = ClimaComms.bcast(context, output_dir)
            dict_writer = ClimaDiagnostics.Writers.DictWriter()
            if (space isa ClimaCore.Spaces.PointSpace) &&
               pkgversion(ClimaCore) < v"0.14.27"
                @test_throws "HDF5Writer only supports Fields with PointSpace for ClimaCore >= 0.14.27" setup_integrator(
                    output_dir;
                    context,
                    space,
                    dict_writer,
                )
            else
                integrator =
                    setup_integrator(output_dir; context, space, dict_writer)

                SciMLBase.solve!(integrator)

                if ClimaComms.iamroot(context)
                    NCDatasets.NCDataset(
                        joinpath(output_dir, "YO_1it_inst.nc"),
                    ) do nc
                        @test nc["YO"].attrib["short_name"] == "YO"
                        @test nc["YO"].attrib["long_name"] ==
                              "YO YO, Instantaneous"
                        @test size(nc["YO"]) == (11, written_space_dims...)
                        @test nc["YO"].attrib["start_date"] ==
                              string(Dates.DateTime(2015, 2, 2))
                    end

                    NCDatasets.NCDataset(
                        joinpath(output_dir, "YO_2it_average.nc"),
                    ) do nc
                        @test nc["YO"].attrib["short_name"] == "YO"
                        @test nc["YO"].attrib["long_name"] ==
                              "YO YO, average within every 2 iterations"
                        @test size(nc["YO"]) == (5, written_space_dims...)
                    end

                    NCDatasets.NCDataset(
                        joinpath(output_dir, "YO_3s_inst.nc"),
                    ) do nc
                        @test nc["YO"].attrib["short_name"] == "YO"
                        @test nc["YO"].attrib["long_name"] ==
                              "YO YO, Instantaneous"
                        @test size(nc["YO"]) == (4, written_space_dims...)
                    end
                end
                @test count(
                    occursin.(
                        Ref(r"YO_1it_inst_\d*\.\d\.h5"),
                        readdir(output_dir),
                    ),
                ) == 11
                reader = ClimaCore.InputOutput.HDF5Reader(
                    joinpath(output_dir, "YO_1it_inst_10.0.h5"),
                    context,
                )
                @test parent(
                    ClimaCore.InputOutput.read_field(reader, "YO_1it_inst"),
                ) == parent(integrator.u.my_var)
                close(reader)
                @test length(keys(dict_writer.dict["YO_1it_inst"])) == 11
                @test dict_writer.dict["YO_1it_inst"][integrator.t] ==
                      integrator.u.my_var
            end
        end
    end
end

@testset "Performance" begin
    mktempdir() do output_dir
        output_dir = ClimaComms.bcast(context, output_dir)

        # Flame
        integrator = setup_integrator(output_dir; context)
        prof = Profile.@profile SciMLBase.solve!(integrator)
        ClimaComms.iamroot(context) && (results = Profile.fetch())
        ClimaComms.iamroot(context) &&
            ProfileCanvas.html_file("flame.html", results)

        # Allocations
        integrator = setup_integrator(output_dir; context)
        prof = Profile.Allocs.@profile SciMLBase.solve!(integrator)
        ClimaComms.iamroot(context) && (results = Profile.Allocs.fetch())
        ClimaComms.iamroot(context) &&
            (allocs = ProfileCanvas.view_allocs(results))
        ClimaComms.iamroot(context) &&
            ProfileCanvas.html_file("allocs.html", allocs)
    end
end

import ClimaComms, MPI
using Logging, Test, LoggingExtras

ctx = ClimaComms.context()
(mypid, nprocs) = ClimaComms.init(ctx)

@testset "FileLogger" begin
    ClimaComms.with_tempdir(ctx) do log_dir
        io = IOBuffer()
        logger = ClimaComms.FileLogger(io, ctx, log_dir)
        fname = ClimaComms.iamroot(ctx) ? "output.log" : "logs/rank_$mypid.log"
        with_logger(logger) do
            test_str = "Test message from rank $mypid"
            @info test_str
            log_content = read(joinpath(log_dir, fname), String)
            @test occursin(test_str, log_content)
            # Test writing to IOBuffer
            if ClimaComms.iamroot(ctx)
                @test occursin(test_str, String(take!(io)))
            else
                @test isempty(String(take!(io)))
            end
        end
    end
end

@testset "MPILogger" begin
    io = IOBuffer()
    logger = ClimaComms.MPILogger(io, ctx)
    with_logger(logger) do
        @info "smoke test"
    end
    str = String(take!(io))
    @test contains(str, "[P$mypid]  Info: smoke test\n")

    # Test with file IOStream
    test_filename, io = mktemp()
    logger = ClimaComms.MPILogger(io, ctx)
    with_logger(logger) do
        test_str = "Test message from rank $mypid"
        @info test_str
        flush(io)
        close(io)

        log_content = read(test_filename, String)
        @test occursin(test_str, log_content)
        @test occursin(test_str, log_content)
    end
end

@testset "OnlyRootLogger" begin
    logger = ClimaComms.OnlyRootLogger(ctx)
    @test typeof(ClimaComms.OnlyRootLogger()) == typeof(logger)
    if ClimaComms.iamroot(ctx)
        @test logger isa Logging.ConsoleLogger
    else
        @test logger isa Logging.NullLogger
    end
end

io = IOBuffer()
summary(io, ctx)
summary_str = String(take!(io))
print(summary_str)

@testset "ClimaComms Summary Tests" begin
    if ClimaComms.iamroot(ctx)
        @test contains(summary_str, string(nameof(typeof(ctx))))
        @test contains(
            summary_str,
            string(nameof(typeof(ClimaComms.device(ctx)))),
        )
    end

    if ctx isa ClimaComms.MPICommsContext
        @testset "MPI Context Tests" begin
            ClimaComms.iamroot(ctx) &&
                @test contains(summary_str, "Total Processes: $nprocs")
            @test contains(summary_str, "Rank: $(mypid-1)")
        end
    end
end

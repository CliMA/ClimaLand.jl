@testset "warnings" begin
    @testset "exclusivity warning" begin
        withenv("JULIA_DEBUG" => "all", "SUPPRESS_THREADINGUTILITIES_WARNING" => "0", "JULIA_EXCLUSIVE" => "") do
            @test (@test_logs (:debug, "The JULIA_EXCLUSIVE environment variable is not set to 1. Therefore, the kernel is allowed to move threads to different cores. We recommend that you set the JULIA_EXCLUSIVE environment variable to 1. To suppress this warning, set the SUPPRESS_THREADINGUTILITIES_WARNING environment variable to 1.") match_mode=:any ThreadingUtilities._print_exclusivity_warning()) isa Nothing
        end
        withenv("SUPPRESS_THREADINGUTILITIES_WARNING" => "1", "JULIA_EXCLUSIVE" => "") do
            @test ThreadingUtilities._print_exclusivity_warning() isa Nothing
        end
        withenv("JULIA_DEBUG" => "all", "SUPPRESS_THREADINGUTILITIES_WARNING" => "0", "JULIA_EXCLUSIVE" => "0") do
            @test (@test_logs (:debug, "The JULIA_EXCLUSIVE environment variable is not set to 1. Therefore, the kernel is allowed to move threads to different cores. We recommend that you set the JULIA_EXCLUSIVE environment variable to 1. To suppress this warning, set the SUPPRESS_THREADINGUTILITIES_WARNING environment variable to 1.") match_mode=:any ThreadingUtilities._print_exclusivity_warning()) isa Nothing
        end
        withenv("SUPPRESS_THREADINGUTILITIES_WARNING" => "1", "JULIA_EXCLUSIVE" => "0") do
            @test ThreadingUtilities._print_exclusivity_warning() isa Nothing
        end
    end

    @testset "utility functions" begin
        @test ThreadingUtilities._string_to_bool("true")
        @test ThreadingUtilities._string_to_bool("1")
        @test !ThreadingUtilities._string_to_bool("false")
        @test !ThreadingUtilities._string_to_bool("0")
        @test !ThreadingUtilities._string_to_bool("")
    end
end

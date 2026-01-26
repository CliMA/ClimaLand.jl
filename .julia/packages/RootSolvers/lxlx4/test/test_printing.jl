using Printf

@testset "Solution pretty printing" begin
    # Test that solution objects display correctly with pretty printing
    # This validates the show methods for both solution types

    @testset "CompactSolution" begin
        # Test successful convergence case for CompactSolution
        # This validates the basic display format and convergence status
        sol = find_zero(
            x -> x^2 - 100^2,
            SecantMethod{Float64}(0.0, 1000.0),
            CompactSolution(),
        )
        sol_str = sprint(show, sol)  # Capture the string output of show(sol)

        # Validate the display format for successful convergence
        @test startswith(sol_str, "CompactSolutionResults{Float64}")  # Check type header
        @test contains(sol_str, "converged")                          # Check convergence status


        # test compact printing
        sol_str = sprint(show, sol, context = :compact => true)
        @test count(sol_str, "\n") == 0  # Ensure no newlines in compact mode

        # Test failed convergence case for CompactSolution
        # This validates display format when method doesn't converge
        sol = find_zero(
            x -> x^2 - 100^2,
            SecantMethod{Float64}(0.0, 1e3),
            CompactSolution(),
            RelativeSolutionTolerance(eps(10.0)),  # Very strict tolerance
            2,
        )                                      # Only 2 iterations (too few)
        sol_str = sprint(show, sol)

        # Validate the display format for failed convergence
        @test startswith(sol_str, "CompactSolutionResults{Float64}")  # Check type header
        @test contains(sol_str, "failed to converge")                 # Check failure status
    end

    @testset "VerboseSolution" begin
        # Test successful convergence case for VerboseSolution
        # This validates the detailed display format with iteration history
        sol = find_zero(
            x -> x^2 - 100^2,
            BrentsMethod{Float64}(-10.0, 50000.0),
            VerboseSolution(),
        )
        sol_str = sprint(show, sol)

        # Validate the comprehensive display format for successful convergence
        @test startswith(sol_str, "VerboseSolutionResults{Float64}")  # Check type header
        @test contains(sol_str, "converged")                          # Check convergence status
        @test contains(sol_str, "Root: $(sol.root)")                  # Check root value display
        @test contains(sol_str, "Error: $(sol.err)")                  # Check error value display
        @test contains(sol_str, "Iterations: $(length(sol.root_history)-1)")  # Check iteration count
        @test contains(sol_str, "History")                            # Check history section header
        history_lines = split(sol_str, "History:\n")[2]
        @test count("\n", history_lines) <= 21 # Check that long histories are truncated


        # Test failed convergence case for VerboseSolution
        # This validates detailed display format when method doesn't converge
        sol = find_zero(
            x -> x^2 - 100^2,
            SecantMethod{Float64}(0.0, 1e3),
            VerboseSolution(),
            RelativeSolutionTolerance(eps(10.0)),  # Very strict tolerance
            2,
        )                                      # Only 2 iterations (too few)
        sol_str = sprint(show, sol)

        # Validate the display format for failed convergence
        @test startswith(sol_str, "VerboseSolutionResults{Float64}")  # Check type header
        @test contains(sol_str, "failed to converge")                 # Check failure status
    end
end

@testset "THREADPOOL" begin
  @test isconst(ThreadingUtilities, :THREADPOOL) # test that ThreadingUtilities.THREADPOOL is a constant
  @test ThreadingUtilities.THREADPOOL isa Vector{UInt}
  @test eltype(ThreadingUtilities.THREADPOOL) === UInt
  sys_threads::Int = parse(Bool, get(ENV, "GITHUB_ACTIONS", "false")) ? Threads.nthreads() : (Sys.CPU_THREADS)::Int
  @test length(ThreadingUtilities.THREADPOOL) == (ThreadingUtilities.THREADBUFFERSIZE÷sizeof(UInt)) * (min(Threads.nthreads(),sys_threads) - 1) + (256 ÷ sizeof(UInt)) - 1
end

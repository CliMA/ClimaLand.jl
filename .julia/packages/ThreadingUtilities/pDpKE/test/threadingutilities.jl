struct Copy{P} end
function (::Copy{P})(p::Ptr{UInt}) where {P}
  _, (ptry, ptrx, N) = ThreadingUtilities.load(p, P, 2 * sizeof(UInt))
  N > 0 || throw("This function throws if N == 0 for testing purposes.")
  @simd ivdep for n ∈ 1:N
    unsafe_store!(ptry, unsafe_load(ptrx, n), n)
  end
  ThreadingUtilities._atomic_store!(p, ThreadingUtilities.SPIN)
end
@generated function copy_ptr(::A, ::B) where {A,B}
  c = Copy{Tuple{A,B,Int}}()
  quote
    @cfunction($c, Cvoid, (Ptr{UInt},))
  end
end
function setup_copy!(p, y, x)
  N = length(y)
  @assert length(x) == N
  py = pointer(y)
  px = pointer(x)
  fptr = copy_ptr(py, px)
  offset = ThreadingUtilities.store!(p, fptr, sizeof(UInt))
  ThreadingUtilities.store!(p, (py, px, N), offset)
end

@inline launch_thread_copy!(tid, y, x) = ThreadingUtilities.launch(setup_copy!, tid, y, x)

function test_copy(tid, N = 100_000)
  a = rand(N)
  b = rand(N)
  c = rand(N)
  x = similar(a) .= NaN
  y = similar(b) .= NaN
  z = similar(c) .= NaN
  GC.@preserve a b c x y z begin
    launch_thread_copy!(tid, x, a)
    yield()
    @assert !ThreadingUtilities.wait(tid)
    launch_thread_copy!(tid, y, b)
    yield()
    @assert !ThreadingUtilities.wait(tid)
    launch_thread_copy!(tid, z, c)
    yield()
    @assert !ThreadingUtilities.wait(ThreadingUtilities.taskpointer(tid))
  end
  @test a == x
  @test b == y
  @test c == z
end

@testset "ThreadingUtilities.jl" begin
  for tid ∈ eachindex(ThreadingUtilities.TASKS)
    @test unsafe_load(Ptr{UInt32}(ThreadingUtilities.taskpointer(tid))) == 0x00000001
  end
  @test all(eachindex(ThreadingUtilities.TASKS)) do tid
    ThreadingUtilities.load(
      ThreadingUtilities.taskpointer(tid),
      ThreadingUtilities.ThreadState,
    ) === ThreadingUtilities.WAIT
  end
  @test all(eachindex(ThreadingUtilities.TASKS)) do tid
    ThreadingUtilities._atomic_load(
      reinterpret(Ptr{UInt32}, ThreadingUtilities.taskpointer(tid)),
    ) === reinterpret(UInt32, ThreadingUtilities.WAIT)
  end
  foreach(test_copy, eachindex(ThreadingUtilities.TASKS))

  x = rand(UInt, 3)
  GC.@preserve x begin
    ThreadingUtilities._atomic_store!(pointer(x), zero(UInt))
    @test ThreadingUtilities._atomic_xchg!(pointer(x), ThreadingUtilities.WAIT) ==
          ThreadingUtilities.TASK
    @test ThreadingUtilities._atomic_umax!(pointer(x), ThreadingUtilities.TASK) ==
          ThreadingUtilities.WAIT
    @test ThreadingUtilities._atomic_umax!(pointer(x), ThreadingUtilities.SPIN) ==
          ThreadingUtilities.WAIT
    @test ThreadingUtilities.load(pointer(x), ThreadingUtilities.ThreadState) ==
          ThreadingUtilities.SPIN
  end
  # Make all tasks error
  for tid ∈ eachindex(ThreadingUtilities.TASKS)
    launch_thread_copy!(tid, Float64[], Float64[])
  end
  sleep(1)
  @test all(istaskfailed, ThreadingUtilities.TASKS)
  # Test that `wait` reports the error for each task
  for tid in eachindex(ThreadingUtilities.TASKS)
    @test_throws TaskFailedException ThreadingUtilities.wait(tid)
  end
  # Test that none of the tasks are in the failed state
  @test !any(istaskfailed, ThreadingUtilities.TASKS)
  # test copy on the reinitialized tasks
  foreach(test_copy, eachindex(ThreadingUtilities.TASKS))
end

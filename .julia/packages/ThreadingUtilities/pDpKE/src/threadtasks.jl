struct ThreadTask
  p::Ptr{UInt}
end
Base.pointer(tt::ThreadTask) = tt.p

@inline taskpointer(tid::T) where {T} = THREADPOOLPTR[] + tid * (THREADBUFFERSIZE % T)

@inline function _call(p::Ptr{UInt})
  fptr = load(p + sizeof(UInt), Ptr{Cvoid})
  assume(fptr ≠ C_NULL)
  ccall(fptr, Cvoid, (Ptr{UInt},), p)
end
@inline function launch(f::F, tid::Integer, args::Vararg{Any,K}) where {F,K}
  p = taskpointer(tid)
  f(p, args...)
  # exchange must happen atomically, to prevent it from switching to `WAIT` after reading
  state = _atomic_xchg!(p, TASK)
  state == WAIT && wake_thread!(tid)
  return nothing
end

function (tt::ThreadTask)()
  p = pointer(tt)
  max_wait = one(UInt32) << 20
  wait_counter = max_wait
  GC.@preserve THREADPOOL begin
    while true
      if _atomic_state(p) == TASK
        _call(p)
        wait_counter = zero(UInt32)
        continue
      end
      pause()
      if (wait_counter += one(UInt32)) > max_wait
        wait_counter = zero(UInt32)
        _atomic_cas_cmp!(p, SPIN, WAIT) && Base.wait()
      end
    end
  end
end

function _sleep(p::Ptr{UInt})
  _atomic_store!(p, WAIT)
  Base.wait()
  return nothing
end

function sleep_all_tasks()
  fptr = @cfunction(_sleep, Cvoid, (Ptr{UInt},))
  for tid ∈ eachindex(TASKS)
    p = taskpointer(tid)
    ThreadingUtilities.store!(p, fptr, sizeof(UInt))
    _atomic_cas_cmp!(p, SPIN, TASK)
  end
  for tid ∈ eachindex(TASKS)
    wait(tid)
  end
end

# 1-based tid, pushes into task 2-nthreads()
@noinline function wake_thread!(_tid::T) where {T<:Integer}
  tid = _tid % Int
  tidp1 = tid + one(tid)
  assume(unsigned(length(Base.Workqueues)) > unsigned(tid))
  assume(unsigned(length(TASKS)) > unsigned(tidp1))
  @inbounds push!(Base.Workqueues[tidp1], TASKS[tid])
  ccall(:jl_wakeup_thread, Cvoid, (Int16,), tid % Int16)
end

@noinline function checktask(tid)
  t = TASKS[tid]
  if istaskfailed(t)
    initialize_task(tid)
    throw(TaskFailedException(t))
  end
  yield()
  false
end

# Like checktask, but without throwing errors.
# This is used in Polyester.reset_threads!().
@noinline function reinit_task(tid)
  t = TASKS[tid]
  if istaskfailed(t)
    initialize_task(tid)
  end
  yield()
end

# 1-based tid
@inline tasktid(p::Ptr{UInt}) = (p - THREADPOOLPTR[]) ÷ (THREADBUFFERSIZE)
@inline wait(tid::Integer) = wait(taskpointer(tid), tid)
@inline wait(p::Ptr{UInt}) = wait(p, tasktid(p))
@inline function wait(p::Ptr{UInt}, tid)
  counter = 0x00000000
  while _atomic_state(p) == TASK
    pause()
    if ((counter += 0x00000001) > 0x00010000)
      checktask(tid) && return true
    end
  end
  false
end

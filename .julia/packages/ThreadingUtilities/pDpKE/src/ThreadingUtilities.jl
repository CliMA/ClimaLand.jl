module ThreadingUtilities
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

using ManualMemory: load, store!

"""
  pause()

For use in spin-and-wait loops, like spinlocks.
"""
@inline pause() = ccall(:jl_cpu_pause, Cvoid, ())

if VERSION ≥ v"1.6.0-DEV.674"
  @inline assume(b::Bool)::Cvoid = Base.llvmcall(("    declare void @llvm.assume(i1)\n\n    define void @entry(i8) alwaysinline {\n    top:\n        %b = trunc i8 %0 to i1\ncall void @llvm.assume(i1 %b)\nret void\n    }\n", "entry"), Cvoid, Tuple{Bool}, b)
else
  @inline assume(b::Bool)::Cvoid = Base.llvmcall(("declare void @llvm.assume(i1)", "%b = trunc i8 %0 to i1\ncall void @llvm.assume(i1 %b)\nret void"), Cvoid, Tuple{Bool}, b)
end

@enum ThreadState::UInt32 begin
  TASK = 0   # 0: task available
  WAIT = 1   # 1: waiting
  SPIN = 2   # 2: spinning
end
const TASKS = Task[]
const LINESPACING = 256 # maximum cache-line size among contemporary CPUs.
const THREADBUFFERSIZE = 512
const THREADPOOL = UInt[]
const THREADPOOLPTR =  Ref{Ptr{UInt}}(C_NULL);

include("atomics.jl")
include("threadtasks.jl")
include("warnings.jl")

function initialize_task(tid)
  _atomic_store!(taskpointer(tid), WAIT)
  t = Task(ThreadTask(taskpointer(tid)));
  t.sticky = true # create and pin
  # set to tid, we have tasks 2...nthread, from 1-based ind perspective
  ccall(:jl_set_task_tid, Cvoid, (Any, Cint), t, tid % Cint)
  TASKS[tid] = t
  return nothing
end

function __init__()
  _print_exclusivity_warning()
  sys_threads::Int = parse(Bool, get(ENV, "GITHUB_ACTIONS", "false")) ? Threads.nthreads() : (Sys.CPU_THREADS)::Int
  nt = min(Threads.nthreads(), sys_threads) - 1
  resize!(THREADPOOL, (THREADBUFFERSIZE ÷ sizeof(UInt)) * nt + (LINESPACING ÷ sizeof(UInt)) - 1)
  copyto!(THREADPOOL, zero(UInt))
  # align to LINESPACING boundary, and then subtract THREADBUFFERSIZE to make the pointer 1-indexed
  THREADPOOLPTR[] = reinterpret(Ptr{UInt}, (reinterpret(UInt, pointer(THREADPOOL))+LINESPACING-1) & (-LINESPACING)) - THREADBUFFERSIZE
  resize!(TASKS, nt)
  foreach(initialize_task, 1:nt)
end

end # module

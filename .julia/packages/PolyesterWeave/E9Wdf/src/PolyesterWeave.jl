module PolyesterWeave
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
  @eval Base.Experimental.@max_methods 1
end

using BitTwiddlingConvenienceFunctions: nextpow2
using ThreadingUtilities: _atomic_store!, _atomic_or!, _atomic_xchg!
using Static
using IfElse: ifelse

export request_threads, free_threads!

@static if VERSION ≥ v"1.6.0-DEV.674"
  @inline function assume(b::Bool)
    Base.llvmcall(
      (
        """
      declare void @llvm.assume(i1)

      define void @entry(i8 %byte) alwaysinline {
      top:
        %bit = trunc i8 %byte to i1
        call void @llvm.assume(i1 %bit)
        ret void
      }
  """,
        "entry"
      ),
      Cvoid,
      Tuple{Bool},
      b
    )
  end
else
  @inline assume(b::Bool) = Base.llvmcall(
    (
      "declare void @llvm.assume(i1)",
      "%b = trunc i8 %0 to i1\ncall void @llvm.assume(i1 %b)\nret void"
    ),
    Cvoid,
    Tuple{Bool},
    b
  )
end

const WORKERS = Ref{NTuple{8,UInt64}}(ntuple(((zero ∘ UInt64)), Val(8)))

include("unsignediterator.jl")
include("request.jl")

dynamic_thread_count() = min((Sys.CPU_THREADS)::Int, Threads.nthreads())
worker_mask_init(x) = if x < 0
    0x0000000000000000
  elseif x ≥ 64
    0xffffffffffffffff
  else
    ((0x0000000000000001 << x) - 0x0000000000000001)
  end
# function static_thread_init()
#   nt = num_threads()
#   Base.Cartesian.@ntuple 8 i -> worker_mask_init(nt - (i-1)*64)
# end
function reset_workers!()
  # workers = ntuple(((zero ∘ UInt64)), Val(8))
  nt = dynamic_thread_count() - 1
  WORKERS[] = Base.Cartesian.@ntuple 8 i -> worker_mask_init(nt - (i - 1) * 64)
end
__init__() = reset_workers!()

include("utility.jl")

end

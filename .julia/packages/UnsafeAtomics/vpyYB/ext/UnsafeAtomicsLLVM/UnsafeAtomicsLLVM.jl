baremodule UnsafeAtomicsLLVM

module Internal

using Core: LLVMPtr
using LLVM
using UnsafeAtomics: UnsafeAtomics, Ordering

include("internal.jl")

end  # module Internal

end  # baremodule UnsafeAtomicsLLVM

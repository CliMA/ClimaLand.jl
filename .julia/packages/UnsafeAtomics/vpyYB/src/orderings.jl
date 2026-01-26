struct LLVMOrdering{name} <: Ordering end

const unordered = LLVMOrdering{:unordered}()
const monotonic = LLVMOrdering{:monotonic}()
const acquire = LLVMOrdering{:acquire}()
const release = LLVMOrdering{:release}()
const acq_rel = LLVMOrdering{:acq_rel}()
const seq_cst = LLVMOrdering{:seq_cst}()

const orderings = (unordered, monotonic, acquire, release, acq_rel, seq_cst)

const ConcreteOrdering = Union{map(typeof, orderings)...}

llvm_ordering(::LLVMOrdering{name}) where {name} = name

Base.string(o::LLVMOrdering) = String(llvm_ordering(o))
Base.print(io::IO, o::LLVMOrdering) = print(io, string(o))

Base.show(io::IO, o::ConcreteOrdering) = print(io, UnsafeAtomics, '.', llvm_ordering(o))

base_ordering(::LLVMOrdering{name}) where {name} = name
base_ordering(::LLVMOrdering{:seq_cst}) = :sequentially_consistent
base_ordering(::LLVMOrdering{:acq_rel}) = :acquire_release
struct LLVMSyncScope{name} <: SyncScope end

const none = LLVMSyncScope{:none}()
const singlethread = LLVMSyncScope{:singlethread}()

const syncscopes = (none, singlethread)
const ConcreteSyncScopes = Union{map(typeof, syncscopes)...}

llvm_syncscope(::LLVMSyncScope{name}) where {name} = name

Base.string(s::LLVMSyncScope) = string("syncscope(\"", llvm_syncscope(s), "\")")
Base.string(s::typeof(none)) = ""

Base.print(io::IO, s::LLVMSyncScope) = print(io, string(s))

Base.show(io::IO, o::ConcreteSyncScopes) = print(io, UnsafeAtomics, '.', llvm_syncscope(o))
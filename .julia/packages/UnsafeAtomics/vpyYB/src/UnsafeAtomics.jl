baremodule UnsafeAtomics

abstract type Ordering end
abstract type SyncScope end

function load end
function store! end
function cas! end
function modify! end
function fence end

function add! end
function sub! end
function xchg! end
function and! end
function nand! end
function or! end
function xor! end
function max! end
function min! end

# =>
right(_, x) = x

module Internal

using Base.Sys: WORD_SIZE
using Base: bitcast, llvmcall

using ..UnsafeAtomics: UnsafeAtomics, Ordering, SyncScope, right

include("utils.jl")
include("orderings.jl")
include("syncscopes.jl")
include("core.jl")

end  # module Internal

const unordered = Internal.unordered
const monotonic = Internal.monotonic
const acquire = Internal.acquire
const release = Internal.release
const acq_rel = Internal.acq_rel
const seq_cst = Internal.seq_cst

# Julia names
const acquire_release = acq_rel
const sequentially_consistent = seq_cst

# SyncScope
const none = Internal.none
const singlethread = Internal.singlethread

end  # baremodule UnsafeAtomics

module TestOrderings

using UnsafeAtomics: UnsafeAtomics, Ordering
using Test

const ORDERINGS = [:unordered, :monotonic, :acquire, :release, :acq_rel, :seq_cst]

function check_ordering(name::Symbol)
    ord = getfield(UnsafeAtomics, name)
    @test ord isa Ordering
    @test string(ord) == sprint(print, ord) == string(name)
    @test repr(ord) == "UnsafeAtomics.$name"
end

function test_ordering()
    @testset for name in ORDERINGS
        check_ordering(name)
    end
end

function test_aliases()
    @test UnsafeAtomics.acq_rel === UnsafeAtomics.acquire_release
    @test UnsafeAtomics.seq_cst === UnsafeAtomics.sequentially_consistent
end

end  # module

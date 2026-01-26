using Test
using Markdown: MD

@testset "ReturnTypes" begin

    using CpuId

    # Moved upwards temporarily for better diagnostics
    println(cpuinfo())
    println(cpufeaturetable())
    println(hvinfo())
    println()

    # Can't do real testing on results when target machine is unknown.
    # Thus, let's simply check whether the result types are correct,
    # which also fails if a test throws.

    @test isa( CpuId.cpuid(), NTuple{4, UInt32} )
    @test isa( CpuId.cpuid(0x00), NTuple{4, UInt32} )
    @test isa( CpuId.cpuid(0x00, 0x00), NTuple{4, UInt32} )

    # Test the low-level cpufeature querying function
    @test isa( CpuId.cpufeature(CpuId.SSE)      , Bool )
    @test isa( CpuId.cpufeature(:SSE)           , Bool )
    @test isa( CpuId.cpufeaturedesc(:SSE)       , String )
    @test isa( CpuId.cpufeatures()              , Vector{Symbol} )

    @test_throws UndefVarError CpuId.cpufeature(:UNKNOWNFEATURE)
    @test CpuId.cpufeaturedesc(:UNKNOWNFEATURE) ==
                            "Unknown feature, no description available!"

    # LLVM eliminates calls to hasleaf(...) if the executing machine supports
    # that leaf.  Thus test whether the reverse actually throws...
    function test_nonexisting_leaf()
        leaf = 0x8000_008f
        CpuId.hasleaf(leaf) || CpuId._throw_unsupported_leaf(leaf)
        CpuId.cpuid(leaf)
    end
    @test_throws ErrorException test_nonexisting_leaf()

    @test isa( address_size()         , Integer )
    @test isa( physical_address_size(), Integer )
    @test isa( cachelinesize()        , Integer )
    @test isa( cachesize()            , Tuple )
    for i in 0:5
        @test isa( cachesize(i)       , Integer )
    end
    @test isa( cacheinclusive()       , Tuple )
    for i in 0:5
        @test isa( cacheinclusive(i)  , Integer )
    end
    @test isa( cpuarchitecture()      , Symbol )
    @test isa( cpubrand()             , String )
    @test isa( cpumodel()             , Dict )
    @test isa( cpuvendor()            , Symbol )
    @test isa( CpuId.cpuvendorstring(), String )
    @test isa( hypervised()           , Bool )
    @test isa( hvvendor()             , Symbol )
    @test isa( CpuId.hvvendorstring() , String )
    @test isa( hvversion()            , Dict{Symbol,Any} )
    @test isa( hypervised()           , Bool )
    @test isa( simdbits()             , Integer )
    @test isa( simdbytes()            , Integer )
    @test isa( cpucores()             , Integer )
    @test isa( cputhreads()           , Integer )
    @test isa( has_cpu_frequencies()  , Bool )
    @test isa( cpu_base_frequency()   , Integer )
    @test isa( cpu_bus_frequency()    , Integer )
    @test isa( cpu_max_frequency()    , Integer )
    @test isa( cpuinfo()              , MD )
    @test isa( cpufeaturetable()      , MD )
    @test isa( hvinfo()               , MD )

    @test isa( cpucycle()             , UInt64 )
    @test isa( cpucycle_id()          , Tuple{UInt64,UInt64} )

    @test isa( perf_revision()        , Int )
    @test isa( perf_gen_counters()    , Int )
    @test isa( perf_gen_bits()        , Int )
    @test isa( perf_fix_counters()    , Int )
    @test isa( perf_fix_bits()        , Int )

    # Check if trailing null characters are correctly identified
    # as hypervisor vendor KVM
    @test get( CpuId._cpuid_vendor_id, "KVMKVMKVM\0\0\0", :Unknown) === :KVM

    # If we're on Linux, then also dump /proc/cpuinfo for comparison when on a
    # remote CI, but only print data of the first CPU.
    Sys.islinux() && run(`sed -e '/^$/,$d' /proc/cpuinfo`)

end

print("\n\n-----\nMocking CpuId\n-----\n\n")
flush(stdout) ; flush(stderr)

include("mock.jl")
include("mockdb.jl")

dump_cpuid_table() ; flush(stdout) ; flush(stderr)

print("\n\n-----\nMocking CpuId\n-----\n\n")
flush(stdout) ; flush(stderr)

# Run the known cpuid records
@testset "Mocking" begin
    for i in 1:length(_mockdb)
        # temporarily replace the low-level cpuid function with known records
        eval(quote
            @testset "Mocked #$($i) $(strip(cpubrand()))" begin
                mock_cpuid($i)
                flush(stdout) ; flush(stderr)
                @test isa( cpubrand()       , String )
                @test isa( cpuinfo()        , MD )
                @test isa( cpufeaturetable(), MD )
                @test isa( hvinfo()         , MD )
                println("Tested recorded cpuid table #",$i," for '", strip(cpubrand()), "'")

                for (fns, res) in last(_mockdb[$i])
                    if isa(res, Tuple) && length(res) > 0 && first(res) === :broken
                        @test_broken getfield(CpuId, fns)() == last(res)
                    else
                        @test getfield(CpuId, fns)() == res
                    end
                end

            end
        end)
        flush(stdout) ; flush(stderr)
    end
end

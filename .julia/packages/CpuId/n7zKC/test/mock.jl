#=--- CpuId / test / mock.jl ----------------------------------------------=#

#
# This file contains a number of cpuid dumps from known CPUs.
# These dumps are used to test the analysis function of CpuId.jl.
#

import CpuId
import CpuId.CpuInstructions

"""
Entry for mocking `cpuid`, maps input and output of the `cpuid` instruction
on a known CPU; viz. `(eax, ecx) -> (eax, ebx, ecx, edx)`.
"""
const _mockdb_entry = Pair{ Dict{NTuple{2,UInt32}, NTuple{4,UInt32}}, Dict{Symbol,Any} }


"""
List of cpuid dumps.
Each entry is represented by a dictionary which maps input and output of the
`cpuid` instruction on a known CPU; viz. `(eax, ecx) -> (eax, ebx, ecx, edx)`.
"""
const _mockdb = _mockdb_entry[]


"""
Temporarily overwrite the low-level cpuid function, to return pre-recorded results.
"""
function mock_cpuid(idx::Integer)

    global _mockdb

    !(1 ≤ idx ≤ length(_mockdb)) &&
        error("The cpuid mocking database does only have entries " *
              "with indices 1:$(length(_mockdb))")

    _fake_cpuid(eax=0, ecx=0)::NTuple{4,UInt32} =
        get( first(_mockdb[idx]), (UInt32(eax), UInt32(ecx))
           , (zero(UInt32), zero(UInt32), zero(UInt32), zero(UInt32),) )

    CpuInstructions.eval( :(cpuid(eax=0, ecx=0) = $_fake_cpuid(eax, ecx) ) )
end


"""
Dump a table of cpuid answers for mocking.  Ouput in a copy & paste ready format.
"""
function dump_cpuid_table()

    println("Dump of all cpuid leafs, used for mocking 'cpuid'...\n")

    println("# ", strip(CpuId.cpubrand()), " with",
            CpuId.hypervised() ? " "*string(CpuId.hvvendor()) : "out", " hypervisor" )
    println("push!( _mockdb, (Dict(")

    for minleaf in [0x0000_0000,0x2000_0000,0x4000_0000,0x8000_0000]
        # get the maximum leaf
        maxleaf = first(CpuId.cpuid(minleaf))
        maxleaf < minleaf && continue
        for leaf in minleaf:maxleaf
            println("    (",(leaf, 0x0000_0000), " => ", CpuId.cpuid(leaf), "),")
        end
    end

    leaf = 0x0000_0004
    if CpuId.hasleaf(leaf)
        for subleaf in 0x0000_0001:0x0000_000f
            eax, ebx, ecx, edx = CpuId.cpuid(leaf, subleaf)
            println("    (",(leaf, subleaf), " => ", (eax, ebx, ecx, edx), "),")
            eax & 0x1f == 0 && break
        end
    end

    leaf = 0x0000_000b
    if CpuId.hasleaf(leaf)
        for subleaf in 0x0000_0001:0x0000_000f
            eax, ebx, ecx, edx = CpuId.cpuid(leaf, subleaf)
            println("    (",(leaf, subleaf), " => ", (eax, ebx, ecx, edx), "),")
            ebx & 0xffff == 0x0000 && break
        end
    end

    leaf = 0x8000_001d
    if CpuId.hasleaf(leaf)
        for subleaf in 0x0000_0001:0x0000_000f
            eax, ebx, ecx, edx = CpuId.cpuid(leaf, subleaf)
            println("    (",(leaf, subleaf), " => ", (eax, ebx, ecx, edx), "),")
            eax & 0x1f == 0 && break
        end
    end

    println("  ) => Dict{Symbol,Any}(")
    # print the results of certain identification functions
    println("    :cpuvendor       => :", cpuvendor(),       ",")
    println("    :cpuarchitecture => :", cpuarchitecture(), ",")
    println("    :cpucores        => ", cpucores(),         ",")
    println("    :cputhreads      => ", cputhreads(),       ",")
    println("    :cachesize       => ", cachesize(),        ",")
    println("    :cachelinesize   => ", cachelinesize(),    ",")
    println("    :simdbits        => ", simdbits(),         ",")

    println("  )))\n\nDone.\n")

end

#=--- end of file ---------------------------------------------------------=#

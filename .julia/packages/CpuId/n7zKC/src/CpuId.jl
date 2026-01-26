#=--- CpuId / CpuId.jl ----------------------------------------------------=#

"""
# Module CpuId

Query information about and directly from your CPU.
"""
module CpuId

export cpuvendor, cpubrand, cpumodel, cachesize, cachelinesize,
       simdbytes, simdbits, address_size, physical_address_size,
       cpu_base_frequency, cpu_max_frequency, cpu_bus_frequency,
       has_cpu_frequencies, hypervised, hvvendor, hvversion,
       hvinfo, cpucores, cputhreads, cpucycle, cpucycle_id,
       perf_revision, perf_fix_counters, perf_fix_bits, perf_gen_counters,
       perf_gen_bits, cpuinfo, cpufeature, cpufeatures, cpufeaturedesc,
       cpufeaturetable, cpuarchitecture, cacheinclusive

using Markdown: MD, Table, parse
const MarkdownString = MD     # Rename Markdown constructors
const MarkdownTable = Table   # to avoid deprecation warning in 0.7-beta
const parse_markdown = parse

using Base: @_noinline_meta, @_inline_meta

# Particular feature flag query is also externalized due to largeness of dicts.
include("cpufeature.jl")

# The low-level assembly functions are confined to a helper module.
include("CpuInstructions.jl")

using .CpuInstructions: cpuid, rdtsc, rdtscp


"""
Helper function, tagged noinline to not have detrimental effect on performance.
"""
@noinline _throw_unsupported_leaf(leaf) =
    error("This CPU does not provide information on cpuid leaf 0x$(string(leaf, base=16, pad=8)).")


"""
    hasleaf(leaf::UInt32) ::Bool

Helper function (not exported) to test whether the CPU claims to provide the
given leaf in a `cpuid` instruction call.

Note: It appears LLVM really know its gear: If this function is inlined, and
      just-in-time compiled, then this test is eliminated completly if the
      executing machine does support this feature. Yeah!
"""
function hasleaf(leaf::UInt32) ::Bool
    @_inline_meta
    eax, ebx, ecx, edx = cpuid(leaf & 0xffff_0000)
    eax >= leaf
end


"""
    cpucycle()

Read the CPU's [Time Stamp Counter, TSC](https://en.wikipedia.org/wiki/Time_Stamp_Counter),
directly with a `rdtsc` instruction.  This counter is increased for every CPU
cycle, until reset.  This function has, when inlined, practically no overhead
and is, thus, probably the fasted method to count how many cycles the CPU has
spent working since last read.

Note, the TSC runs at a constant rate if `hasfeature(:TSCINV)==true`;
otherwise, it is tied to the current CPU clock frequency.

Hint: This function is extremely efficient when inlined into your own code.
      Convince yourself by typing `@code_native CpuId.cpucycle()`.
      To use this for benchmarking, simply subtract the results of two calls.
      The result is the actual CPU clock cycles spent, independent of the
      current (and possible non-constant) CPU clock frequency.
"""
function cpucycle end
@eval cpucycle() = $(cpufeature(TSC)) ? rdtsc() : zero(UInt64)


"""
    cpucycle_id()

Read the CPU's [Time Stamp Counter, TSC](https://en.wikipedia.org/wiki/Time_Stamp_Counter),
and executing CPU id directly with a `rdtscp` instruction.  This function is
similar to the `cpucycle()`, but uses an instruction that also allows to
detect if the code has been moved to a different executing CPU.  See also the
comments for `cpucycle()` which equally apply.
"""
function cpucycle_id end
@eval cpucycle_id() = $(cpufeature(RDTSCP)) ? rdtscp() : (zero(UInt64),zero(UInt64))


"""
    cpumodel()

Obtain the CPU model information as a Dict of pairs of
`:Family`, `:Model`, `:Stepping`, and `:CpuType`.
"""
function cpumodel() ::Dict{Symbol, UInt8}
    #  Stepping:    0: 3
    #  Model:       4: 7
    #  Family:      8:11
    #  Processor:  12:13
    #  Ext.Model:  16:19
    #  Ext.Family: 20:27
    eax, ebx, ecx, edx = cpuid(0x01)
    Dict( :Family    => UInt8(((eax & 0x0000_0f00) >> 8)
                             +((eax & 0x0ff0_0000) >> (20-4)))
        , :Model     => UInt8(((eax & 0x0000_00f0) >> 4)
                             +((eax & 0x000f_0000) >> (16-4)))
        , :Stepping  =>  UInt8(eax & 0x0000_000f)
        , :CpuType   =>  UInt8((eax & 0x0000_3000) >> 12))
end


"""
    hypervised()

Check whether the CPU reports to run a hypervisor context, that is,
whether the current process runs in a virtual machine.

A positive answer may indicate that other information reported by the CPU
is fake, such as number of physical and logical cores.  This is because
the hypervisor is free to decide which information to pass.
"""
function hypervised() ::Bool
    # alternative: 0x8000_000a, eax bit 8 set if hv present.
    eax, ebx, ecx, edx = cpuid(0x0000_0001)
    ((ecx >> 31) & one(UInt32)) != zero(UInt32)
end


"""
    hvvendorstring()

Determine the hypervisor vendor string as provided by the cpu by executing a
`cpuid` instruction.  Note, this string has a fixed length of 12 characters.
Use `hvvendor()` if you prefer getting a parsed Julia symbol.  If the CPU is
not running a hypervisor, a string of undefined content will be returned.
"""
function hvvendorstring()
    eax, ebx, ecx, edx = cpuid(0x4000_0000)
    String( reinterpret(UInt8, [ebx, ecx, edx] ) )
end


"""
    hvversion()

Get a dictionary with additional information of the running hypervisor.
The dictionary is empty if no hypervisor is detected, and only tags that are
provided by the hypervisor are included.

Note, the data available is hypervisor vendor dependent.
"""
function hvversion()

    d = Dict{Symbol,Any}()
    !hypervised() && return d

    vendor = d[:vendor] = hvvendor()

    # Signature and version info appear to be vendor specific.

    if vendor == :Microsoft

        # Specs see e.g.
        # https://msdn.microsoft.com/en-us/library/windows/hardware/dn613994(v=vs.85).aspx

        leaf = 0x4000_0001
        if hasleaf(leaf)
            eax, ebx, ecx, edx = cpuid(leaf)
            eax != 0x00 && (d[:signature] = String( reinterpret(UInt8, ( [eax, ] ))))
        end

        leaf = 0x4000_0002
        if hasleaf(leaf)
            eax, ebx, ecx, edx = cpuid(leaf)
            if eax != 0x00 && ebx != 0x00
                d[:version] = VersionNumber(
                                    Int((ebx >> 16) & 0xffff),
                                    Int( ebx        & 0xffff),
                                    Int(eax))
            end
            if ecx != 0x00 && edx != 0x00
                d[:servicepack] = VersionNumber(
                                    Int( ecx ),                     # "servicepack"
                                    Int((edx >> 24) & 0x0000_00ff), # "branch"
                                    Int( edx        & 0x00ff_ffff)) # "number"
            end
        end

        return d

    elseif vendor == :Xen

        # Xen is specified e.g. here
        # https://xenbits.xen.org/docs/unstable/hypercall/x86_64/include,public,arch-x86,cpuid.h.html

        leaf = 0x4000_0001
        if hasleaf(leaf)
            eax, ebx, ecx, edx = cpuid(leaf)
            d[:version] = VersionNumber(Int((eax >> 16) & 0xffff), Int(eax & 0xffff))
        end

        return d

    elseif vendor == :VMware

        # Specs see
        # https://kb.vmware.com/selfservice/microsites/search.do?language=en_US&cmd=displayKC&externalId=1009458
        # VmWare seems to provide frequency information, but nothing else valuable

        leaf = 0x4000_0010
        if hasleaf(leaf)
            eax, ebx, ecx, edx = cpuid(leaf)
            eax != 0x00 && (d[:tscfreq] = Int(eax))
            ebx != 0x00 && (d[:busfreq] = Int(ebx))
        end

        return d

    elseif vendor == :KVM

        return d

    end

    return d
end


"""
    hvinfo() ::MarkdownString

Generate a markdown table of all the detected/available/supported tags of a
running hypervisor.  If there is no hosting hypervisor, an empty markdown
string is returned.
"""
function hvinfo()
    d = hvversion()
    isempty(d) && return MarkdownString()

    md = "| Hypervisor | Value |\n|:------|:------|\n"

    haskey(d, :vendor) &&
        (md *= string( "| Vendor      | '", d[:vendor], "' |\n"))

    haskey(d, :signature) &&
        (md *= string( "| Signature   | '", d[:signature], "' |\n"))

    haskey(d, :version) &&
        (md *= string( "| Version     | ", d[:version], " |\n"))

    haskey(d, :servicepack) &&
        (md *= string( "| Servicepack | ", d[:servicepack], " |\n"))

    haskey(d, :tscfreq) &&
        (md *= string( "| Frequencies | TSC = ", d[:tscfreq] ÷ 1000,
                       " MHz, bus = ", d[:busfreq] ÷ 1000, " MHz |\n"))

    parse_markdown(md)
end


"""
    cpuvendorstring()

Determine the cpu vendor string as provided by the cpu by executing a
`cpuid` instruction.  Note, this string has a fixed length of 12 characters.
Use `cpuvendor()` if you prefer getting a parsed Julia symbol.
"""
function cpuvendorstring()
    eax, ebx, ecx, edx = cpuid(0x00)
    String( reinterpret(UInt8, [ebx, edx, ecx] ) )
end


"Map vendor string of type 'char[12]' provided by `cpuid, eax=0x0` to a Julia symbol."
const _cpuid_vendor_id = Dict(
    "AMDisbetter!" => :AMD,
    "AuthenticAMD" => :AMD,
    "CentaurHauls" => :Centaur,
    "CyrixInstead" => :Cyrix,
    "GenuineIntel" => :Intel,
    "Geode by NSC" => :NSC,
    "NexGenDriven" => :NexGen,
    "RiseRiseRise" => :Rise,
    "SiS SiS SiS " => :SiS,
    "TransmetaCPU" => :Transmeta,
    "GenuineTMx86" => :Transmeta,
    "UMC UMC UMC " => :UMC,
    "VIA VIA VIA " => :VIA,
    "Vortex86 SoC" => :Vortex,
    # Hypervisors
    "KVMKVMKVM\0\0\0" => :KVM,   # KVM actually appends null characters...
    "Microsoft Hv" => :Microsoft,
    "prl hyperv  " => :Parallels,
    "VMwareVMware" => :VMware,
    "XenVMMXenVMM" => :Xen,
)


"""
    cpuvendor()

Determine the cpu vendor as a Julia symbol.  In case the CPU vendor
identification is unknown `:Unknown` is returned (then also consider raising
an issue on Github).
"""
cpuvendor() = get(_cpuid_vendor_id, cpuvendorstring(), :Unknown)


"""
    hvvendor()

Determine the hypervisor vendor as a Julia symbol or `:NoHypervisor` if not
running a hypervisor. In case the hypervisor vendor identification is unknown
`:Unknown` is returned (then also consider raising an issue on Github).
"""
hvvendor() = hypervised() ?
                get(_cpuid_vendor_id, hvvendorstring(), :Unknown) :
                :NoHypervisor


"""
    cpubrand()

Determine the cpu brand as a string as provided by the CPU through executing
the `cpuid` instruction.  This function throws if no CPU brand information is
available form the CPU, which should never be the case on recent hardware.
"""
function cpubrand() ::String

    leaf = 0x8000_0004
    hasleaf(leaf) || _throw_unsupported_leaf(leaf)

    # Extract the information from leaf 0x8000_0002..0x8000_0004
    rstrip( String( reinterpret(UInt8,
                    [cpuid(0x8000_0002)..., cpuid(0x8000_0003)..., cpuid(0x8000_0004)..., 0x0000_0000] )
                  )
          , '\0')
end


"""
    cpuarchitecture()

This function tries to infer the CPU microarchitecture with a call to the
`cpuid` instruction.  For now, only Intel CPUs are suppored according to the
following table.  Others are identified as `:Unknown`.

Table C-1 of Intel's Optimization Reference Manual:

| Family_Model                     | Microarchitecture   |
| :------------------------------- | :------------------ |
| 06_4EH, 06_5EH                   | Skylake             |
| 06_3DH, 06_47H, 06_56H           | Broadwell           |
| 06_3CH, 06_45H, 06_46H, 06_3FH   | Haswell             |
| 06_3AH, 06_3EH                   | Ivy Bridge          |
| 06_2AH, 06_2DH                   | Sandy Bridge        |
| 06_25H, 06_2CH, 06_2FH           | Westmere            |
| 06_1AH, 06_1EH, 06_1FH, 06_2EH   | Nehalem             |
| 06_17H, 06_1DH                   | Enhanced Intel Core |
| 06_0FH                           | Intel Core          |
"""
function cpuarchitecture() ::Symbol

    # See also the C++ library VC++ which has quite good detection.
    # https://github.com/VcDevel/Vc/blob/master/cmake/OptimizeForArchitecture.cmake
    #
    # See also Table 35-1 in Intel's Architectures Software Developer Manual.

    cpumod = cpumodel()
    family = cpumod[:Family]
    model  = cpumod[:Model]

    family == 0x06 &&
        return (model == 0x66) ? :Cannonlake :
               (model == 0x8e || model == 0x9e) ? :Kabylake :
               (model == 0x4e || model == 0x5e || model == 0x55) ? :Skylake :
               (model == 0x3d || model == 0x47 || model == 0x4f || model == 0x56) ? :Broadwell :
               (model == 0x3c || model == 0x3f || model == 0x45 || model == 0x46) ? :Haswell :
               (model == 0x3a || model == 0x3e) ? :IvyBridge :
               (model == 0x2a || model == 0x2d) ? :SandyBridge :
               (model == 0x25 || model == 0x2c || model == 0x2f) ? :Westmere :
               (model == 0x1a || model == 0x1e || model == 0x1f || model == 0x2e) ?  :Nehalem :
               (model == 0x17 || model == 0x1d) ?  :EnhancedIntelCore :
               (model == 0x0f || model == 0x1d) ?  :IntelCore :
               # Atom models
               (model == 0x5c) ? :Goldmont :
               (model == 0x5a) ? :Silvermont :
               (model == 0x1c) ? :Atom :
               # Xeon Phi models
               (model == 0x57) ? :KnightsLanding :
               # Well, this is awkward...
               :UnknownIntel

    # Xeon Phi family 0x07, model 0x01, or Itanium ?
    family == 0x07 && return :Itanium

    # AMD types
    family == 0x0f && return :K8
    family == 0x1f && return :K10
    family == 0x2f && return :Griffin       # not confirmed
    family == 0x3f && return :Llano         # not confirmed
    family == 0x5f && return :Bobcat
    family == 0x6f && return :Bulldozer
    family == 0x7f && return (model == 0x30) ? :Puma : :Jaguar
    family == 0x8f && return :Zen

    return :Unknown
end


"""
    cachelinesize()

Query the CPU about the L1 data cache line size in bytes.  This is typically
64 byte.  Returns zero if cache line size information is not available from
the CPU.
"""
function cachelinesize() ::Int

    # If the extended leaf is not available, same information
    # should be available through leaf 0x01, ebx[15:8]

    leaf = 0x8000_0006
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    (ecx & 0xff) % Int
end


"""
    simdbytes()

Query the CPU on the maximum supported SIMD vector size in bytes, or
`sizeof(Int)` if no SIMD capability is reported by the invoked `cpuid`
instruction.
"""
function simdbytes() ::Int
    @_inline_meta

    simd = 0

    if hasleaf(0x0000_0007)
        eax, ebx, ecx, edx = cpuid(0x0000_0007)
        simd = ebx & (1<<16) != 0 ? 512 ÷ 8 :  # AVX512F instruction set
               ebx & (1<< 5) != 0 ? 256 ÷ 8 : 0   # AVX2
    end

    if simd == 0
        eax, ebx, ecx, edx = cpuid(0x0000_0001)
        simd = ecx & (1<<28) != 0 ? 256 ÷ 8 :     # AVX
            edx & (1<<26) != 0 ? 128 ÷ 8 :     # SSE2
            edx & (1<<25) != 0 ? 128 ÷ 8 :     # SSE
            edx & (1<<23) != 0 ?  64 ÷ 8 :     # MMX
            sizeof(Int)
    end

    simd
end


"""
    simdbits()

Query the CPU on the maximum supported SIMD vector size in bits, or
`sizeof(Int)` in bits if no SIMD capability is reported by the invoked `cpuid`
instruction.
"""
simdbits() = simdbytes() * 8


"""
    cpucores()

Determine the number of physical cores on the current executing CPU by
invoking a `cpuid` instruction.  On systems with multiple CPUs, this only
gives information on the single CPU that is executing the code.
Returns zero if querying this feature is not supported, which may also be due
to a running hypervisor (as observed on hvvendor() == :Microsoft).

Also, this function does not take logical cores (aka hyperthreading) into
account, but determines the true number of physical cores, which typically
also share L3 caches and main memory bandwidth.

See also the Julia global variable `Base.Sys.CPU_THREADS`, which gives the
total count of all logical cores on the machine.
"""
function cpucores() ::Int

    leaf = 0x0000_000b
    hasleaf(leaf) || return zero(UInt32)

    # The number of cores reported by cpuid is actually already the total
    # number of cores at that level, including all of the lower levels.
    # Thus, we need to find the highest level...which is 0x02 == "Core"
    # on ecx[15:08] per the specs, and divide it by the number of
    # 0x01 == "SMT" logical cores.

    sl = zero(UInt32)
    nc = zero(UInt32) # "Core" count
    nl = zero(UInt32) # "SMT" count
    while (true)
        # ebx[15:0] must be non-zero according to manual
        eax, ebx, ecx, edx = cpuid(leaf, sl)
        ebx & 0xffff == 0x0000 && break
        sl += one(UInt32)
        lt = ((ecx >> 8) & 0xff) & 0x03
        # Let's assume for now there's only one valid entry for each level
        lt == 0x01 && ( nl = ebx & 0xffff; continue )
        lt == 0x02 && ( nc = ebx & 0xffff; continue )
        # others are invalid and shouldn't be considered..
    end

    return iszero(nc) ? # no cores detected? then maybe its AMD?
        # AMD
        ((cpuid(0x8000_0008)[3] & 0x00ff)+1) :
        # Intel, we need nonzero values of nc and nl
        (iszero(nl) ? nc : nc ÷ nl)
end


"""
    cpunodes() -> Int

Determine the number of core complexes, aka nodes, on this processor.
This notion is introduced by AMD, where L3 caches are shared among the
cores of a comples
"""
function cpunodes()

    # AMD gives the nodes per processor
    # in extended topology information.
    # CPUID[0x8000_001e][ECX][10:8]
    cpunodes_amd() = 1 + ((cpuid(0x8000_001e)[3] >> 8) & 0b0111)

    return 1

end


"""
    cputhreads_per_core() -> Int

Determine the of threads per hardware core on the currently executing CPU.
A value larger than one indicates simulatenous multithreading being enabled,
aka SMT, aka Hyperthreading.
"""
function cputhreads_per_core()

    cputhreads_per_core_amd() =
        # AMD gives the threads per physical core
        # in extended topology information.
        ((cpuid(0x8000_001e)[2] >> 8) & 0x00ff) + 1

    return 1

end

"""
    cputhreads()

Determine the number of logical cores on the current executing CPU by
invoking a `cpuid` instruction.  On systems with multiple CPUs, this only
gives information on the single CPU that is executing the code.
Returns zero if querying this feature is not supported, which may also be due
to a running hypervisor (as observed on hvvendor() == :Microsoft).

In contrast to `cpucores()`, this function also takes logical cores aka
hyperthreading into account.  For practical purposes, only I/O intensive code
should make use of these total number of cores; memory or computation bound
code will not benefit, but rather experience a detrimental effect.

See also Julia's global variable `Base.Sys.CPU_THREADS`, which gives the total
count of all logical cores on the machine.  Thus, `Base.Sys.CPU_THREADS ÷
CpuId.cputhreads()` should give you the number of CPUs (packages) in your
system.
"""
function cputhreads() ::Int

    function cputhreads_amd()
        # AMD stores the total number of threads
        # aka logical processors directly
        ((cpuid(0x0000_0001)[2] >> 16) & 0x00ff)
    end

    # 1) First try to detect whether we have legacy style core count encoding
    #    This is also correct for AMD, but not for modern Intel.
    #    nc = ((cpuid(0x0000_0001)[2] >> 16) % 8)
    #    if !iszero(nc) return nc
    # 2) Try the modern intel extended information

    leaf = 0x0000_000b
    hasleaf(leaf) || return zero(UInt32)

    # The number of cores reported by cpuid is actually already the total
    # number of cores at that level, including all of the lower levels.
    # Thus, we only need to find the highest level...which is 0x02 == "Core"
    # on ecx[15:08] per the specs.

    sl = zero(UInt32)
    nc = zero(UInt32)
    while (true)
        # ebx[15:0] must be non-zero according to manual
        eax, ebx, ecx, edx = cpuid(leaf, sl)
        ebx & 0xffff == 0x0000 && break
        sl += one(UInt32)
        lt = (ecx >> 8) & 0xff  # level type, 0x01 == "SMT", 0x02 == "Core", other are invalid.
        # Let's assume for now there's only one valid entry for each level
        lt != 0x02 && continue
        nc = ebx & 0xffff
    end

    return iszero(nc) ? # no cores detected? then maybe its AMD?
        # AMD
        ((cpuid(0x0000_0001)[2] >> 16) & 0x00ff) :
        # Intel
        (nc)

end

@deprecate cpucores_total() cputhreads()

"""
    address_size()

Determine the maximum virtual address size supported by this CPU as
reported by the `cpuid` instructions.

This information may be used to determine the number of high bits that can be
used in a pointer for tagging;  viz. `sizeof(Int) - address_size() ÷ 8`,
which gives on most 64 bit Intel machines 2 bytes = 16 bit for other purposes.
"""
function address_size() ::Int

    leaf = 0x8000_0008
    hasleaf(leaf) || _throw_unsupported_leaf(leaf)

    eax, ebx, ecx, edx = cpuid(leaf)
    (eax & 0xff00) >> 8
end


"""
    physical_address_size()

Determine the maximum phyiscal addresses size supported by this CPU as
reported by the `cpuid` instructions.  Prefer to make use of `address_size`
for practical purposes; use this only for diagnostic issues, such as
determining the theoretical maximum memory size.  Also note that this address
size is manipulated by a running hypervisor.
"""
function physical_address_size() ::Int

    leaf = 0x8000_0008
    hasleaf(leaf) || _throw_unsupported_leaf(leaf)

    eax, ebx, ecx, edx = cpuid(leaf)
    eax & 0xff
end


"""
    cachesize()
    cachesize(lvl::Integer)

Obtain information on the CPU's *data* cache sizes.

Determine the data cache size for each cache level as reported by the CPU
using a set of calls to the `cpuid` instruction.  Returns a tuple with the
tuple indices matching the cache levels; sizes are given in bytes.

If given an integer, then the data cache size of the respective cache level
will be returned.  This is significantly faster than the tuple version above.

Note that these are total cache sizes, where some cache levels are typically
shared by multiple cpu cores, the higher cache levels may include lower levels.
To print the cache levels in kbyte, use e.g. `CpuId.cachesize() .÷ 1024`.

This functions throws an error if cache size detection is not supported.
"""
function cachesize end


"""
Helper function that performs the actual computation of the cache size with
register values retrieved from `cpuid` on leaf 0x04.

> Cache size information on leaf 0x04 is computed with
>     size in bytes = (ways+1) * (partitions+1) * (linesize+1) * (sets+1)
> where
>     ways = ebx[22:31], partitions = ebx[12:21], linesize = ebx[0:11]
>     sets = ecx[:]
"""
function __datacachesize(eax::UInt32, ebx::UInt32, ecx::UInt32) ::UInt32
    @_inline_meta
    (1 + (ebx>>22) & 0x03ff ) *    # ways
    (1 + (ebx>>12) & 0x03ff ) *    # partitions
    (1 +  ebx      & 0x0fff ) *    # linesize
    (1 +  ecx )                    # sets
end

"""
Helper function to determine the cache size for a given subleaf `sl` on
Intel or AMD Extended.
"""
function __cachesize_level(leaf::UInt32, sl::UInt32)
    eax, ebx, ecx, edx = cpuid(leaf, sl)
    # if eax is zero in the lowest 5 bits, we've reached the sentinel.
    eax & 0x1f == 0 && return ()
    # could do a sanity check: cache level reported in eax bits 5:7
    # if lowest bit on eax is zero, then its not a data cache
    eax & 0x01 == 0 && return __cachesize_level(leaf, sl + one(UInt32))
    # otherwise this should be a valid data or shared cache level
    (signed(__datacachesize(eax, ebx, ecx)), __cachesize_level(leaf, sl + one(UInt32))...)
end

@noinline function cachesize()

    # TODO: This function fails compilation if inlined.

    # Determine the correct leaf id
    std_leaf = 0x0000_0004
    amd_leaf = 0x8000_001d

    leaf =
        if hasleaf(amd_leaf) && cpufeature(TOPX)
            amd_leaf        # AMD Extended Cache
        elseif hasleaf(std_leaf)
            std_leaf        # Default Intel
        else
            _throw_unsupported_leaf(std_leaf)
        end

    return (__cachesize_level(leaf,zero(UInt32))...,)
end

cachesize(lvl::Integer) = cachesize(UInt32(lvl))

function cachesize(lvl::UInt32) ::Int
    @_inline_meta

    lvl == 0 && return 0

    # Determine the correct leaf id
    std_leaf = 0x0000_0004
    amd_leaf = 0x8000_001d

    leaf =
        if hasleaf(amd_leaf) && cpufeature(TOPX)
            amd_leaf        # AMD Extended Cache
        elseif hasleaf(std_leaf)
            std_leaf        # Default Intel
        else
            _throw_unsupported_leaf(std_leaf)
        end

    # Loop over all subleaves of leaf 0x04 until the cache level bits say
    # we've reached to target level, then check whether it's a data cache or
    # shared cache.

    # Assuming cache enumeration is linear, it is sufficient to start at
    # sub-leaf 'lvl'. This should limit looping to at most two iterations.
    # Though it might fail... Let's see if there are bug reports for weird
    # architectures.
    sl = lvl - one(UInt32)

    # this is a variation of __cachesize_level
    while (true)
        eax, ebx, ecx, edx = cpuid(leaf, sl)
        # still at a valid cache level ?
        eax & 0x1f == 0 && return 0
        # is this a data cache or shared cache level, eax[0:4]?
        # and is it the correct level, eax[5:7]?
        (eax & 0x01 == 0x01) && ((eax >> 5) & 0x07) == lvl &&
            return __datacachesize(eax, ebx, ecx)
        # not yet found, thus continue looking at next subleaf
        sl += one(UInt32)
    end
end


"""
    cacheinclusive()
    cacheinclusive(lvl::Integer)

Obtain information on the CPU's *data* cache inclusiveness. Returns `true`
for a cache that is inclusive of the lower cache levels, and `false` otherwise.

Determine the data cache size for each cache level as reported by the CPU
using a set of calls to the `cpuid` instruction.  Returns a tuple with the
tuple indices matching the cache levels.

If given an integer, then the data cache inclusiveness of the respective cache level
will be returned.  This is significantly faster than the tuple version above.
"""
function cacheinclusive end


function cacheinclusive()

    function cacheinclusive_level(leaf, sl::UInt32)
        eax, ebx, ecx, edx = cpuid(leaf, sl)
        # if eax is zero in the lowest 5 bits, we've reached the sentinel.
        eax & 0x1f == 0 && return ()
        # could do a sanity check: cache level reported in eax bits 5:7
        # if lowest bit on eax is zero, then its not a data cache
        eax & 0x01 == 0 && return cacheinclusive_level(leaf, sl + one(UInt32))
        # otherwise this should be a valid data or shared cache level
        ((edx & 0x02) != 0x00, cacheinclusive_level(leaf, sl + one(UInt32))...)
    end

    # TODO: This is awkwardly slow and requires some rework.
    #       Potential approach: Recurse to the last found cache level, there
    #       allocate a small array, then fill the array when leaving each
    #       recursion level.

    leaf = 0x0000_0004
    hasleaf(leaf) && return (cacheinclusive_level(leaf,zero(UInt32))...,)
    # no cache data available
    ()
end

cacheinclusive(lvl::Integer) = cacheinclusive(UInt32(lvl))

function cacheinclusive(lvl::UInt32) ::Int
    @_inline_meta

    lvl == 0 && return false

    leaf = 0x0000_0004
    hasleaf(leaf) || _throw_unsupported_leaf(leaf)

    # Loop over all subleaves of leaf 0x04 until the cache level bits say
    # we've reached to target level, then check whether it's a data cache or
    # shared cache.

    # Assuming cache enumeration is linear, it is sufficient to start at
    # sub-leaf 'lvl'. This should limit looping to at most two iterations.
    # Though it might fail... Let's see if there are bug reports for weird
    # architectures.
    sl = lvl - one(UInt32)

    while (true)
        eax, ebx, ecx, edx = cpuid(leaf, sl)
        # still at a valid cache level ?
        eax & 0x1f == 0 && return false
        # is this a data cache or shared cache level, eax[0:4]?
        # and is it the correct level, eax[5:7]?
        (eax & 0x01 == 0x01) && ((eax >> 5) & 0x07) == lvl &&
            return (edx & 0x02) != 0x00
        sl += one(UInt32)
    end
end


"""
    has_cpu_frequencies()

Determine whether the CPU provides clock frequency information.  If true, then
`cpu_base_frequency()`, `cpu_max_frequency()` and `cpu_bus_frequency()` should
be expected to return sensible information.
"""
function has_cpu_frequencies() ::Bool

    leaf = 0x0000_0016
    hasleaf(leaf) || return false

    # frequencies are provided if any of the bits in question are non-zero
    eax, ebx, ecx = cpuid(leaf)
    (eax & 0xffff) != zero(UInt32) ||
    (ebx & 0xffff) != zero(UInt32) ||
    (ecx & 0xffff) != zero(UInt32)
end


"""
    cpu_base_frequency()

Determine the CPU nominal base frequency in MHz as reported directly from the
CPU through a `cpuid` instruction call.  Returns zero if the CPU doesn't
provide base frequency information.

The actual cpu frequency might be lower due to throttling, or higher due to
frequency boosting (see `cpu_max_frequency`).
"""
function cpu_base_frequency() ::Int

    leaf = 0x0000_0016
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    eax & 0xffff
end


"""
    cpu_max_frequency()

Determine the maximum CPU frequency in MHz as reported directly from the CPU
through a `cpuid` instrauction call.  The maximum frequency typically refers
to the CPU's boost frequency.  Returns zero if the CPU doesn't provide maximum
frequency information.
"""
function cpu_max_frequency() ::Int

    leaf = 0x0000_0016
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    ebx & 0xffff
end


"""
    cpu_bus_frequency()

Determine the bus CPU frequency in MHz as reported directly from the CPU through
a `cpuid` instrauction call.  Returns zero if the CPU doesn't
provide bus frequency information.
"""
function cpu_bus_frequency() ::Int

    leaf = 0x0000_0016
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    ecx & 0xffff
end


"""
    perf_revision()

Determine the revision number of the performance monitoring unit.

This information is only available if `cpufeature(PDCM) == true`.
"""
function perf_revision() ::Int

    leaf = 0x0000_000a
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    eax & 0xff
end


"""
    perf_fix_counters()

Determine the number of fixed-function performance counters on the executing
machine.

This information is only available if `cpufeature(PDCM) == true`.
"""
function perf_fix_counters() ::Int

    leaf = 0x0000_000a
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    edx & 0x1f
end


"""
    perf_gen_counters()

Determine the number of general purpose counters performance counters on the
executing CPU.  Number of counters is given as per logical processor.

This information is only available if `cpufeature(PDCM) == true`.
"""
function perf_gen_counters() ::Int

    leaf = 0x0000_000a
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    (eax >> 8) & 0xff
end


"""
    perf_fix_bits()

Determine the number of bits fixed-function counters performance counters on
the executing CPU.

This information is only available if `cpufeature(PDCM) == true`.
"""
function perf_fix_bits() ::Int

    leaf = 0x0000_000a
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    (edx >> 5) & 0xff
end


"""
    perf_gen_bits()

Determine the number of bits general purpose counters performance counters on
the executing CPU.

This information is only available if `cpufeature(PDCM) == true`.
"""
function perf_gen_bits() ::Int

    leaf = 0x0000_000a
    hasleaf(leaf) || return 0

    eax, ebx, ecx, edx = cpuid(leaf)
    (eax >> 16) & 0xff
end


"""
    cpuinfo()

Generate a markdown table with the results of all of the CPU querying
functions provided by the module `CpuId`.  Intended to give a quick overview
for diagnostic purposes e.g. in log files.
"""
function cpuinfo()

    unsupported = "Not supported by CPU"
    cachesz = cachesize()
    modelfl = cpumodel()

    address = string(address_size(), " bits virtual, ", physical_address_size(), " bits physical")
    cache   = string("Level ", 1:length(cachesz), " : ", map(x->div(x,1024), cachesz), " kbytes")
    cachels = string(cachelinesize(), " byte cache line size")
    cores = string( CpuId.cpucores(), " physical cores, "
                  , CpuId.cputhreads(), " logical cores (on executing CPU)")
    frequencies = !has_cpu_frequencies() ? unsupported :
                        string(cpu_base_frequency(), " / ",
                               cpu_max_frequency(), " MHz (base/max), ",
                               cpu_bus_frequency(), " MHz bus")
    hyperthreading = (CpuId.cpucores() == CpuId.cputhreads() ?  "No " : "") * "Hyperthreading hardware capability detected"
    hypervisor = hypervised() ? "Yes, $(hvvendor())" : "No"
    model = string("Family: 0x",     string(modelfl[:Family],   base=16, pad=2),
                   ", Model: 0x",    string(modelfl[:Model],    base=16, pad=2),
                   ", Stepping: 0x", string(modelfl[:Stepping], base=16, pad=2),
                   ", Type: 0x",     string(modelfl[:CpuType],  base=16, pad=2))
    simd = string(simdbits(), " bit = ", simdbytes(), " byte max. SIMD vector size" )
    tsc = string("TSC is ", (cpufeature(TSC) ? "" : "not "), "accessible via `rdtsc`")
    tscinv = cpufeature(TSCINV) ? "TSC runs at constant rate (invariant from clock frequency)" :
                                   "TSC increased at every clock cycle (non-invariant TSC)"
    perfmon = !cpufeature(PDCM) ? "Performance Monitoring Counters (PMC) are not supported" :
                                  "Performance Monitoring Counters (PMC) revision $(perf_revision())"
    perfmon2 = !cpufeature(PDCM) ?  [] :
        [ ["", "Available hardware counters per logical core:"]
        , ["", "$(perf_fix_counters()) fixed-function counters of $(perf_fix_bits()) bit width"]
        , ["", "$(perf_gen_counters()) general-purpose counters of $(perf_gen_bits()) bit width"] ]
    ibs = !cpufeature(IBS) ? [] : [["", "CPU supports AMD's Instruction Based Sampling (IBS)"]]

    MarkdownString( MarkdownTable( [
        [ "Cpu Property",       "Value"              ],
        #----------------------------------------------
        [ "Brand",              strip(cpubrand())    ],
        [ "Vendor",             cpuvendor()          ],
        [ "Architecture",       cpuarchitecture()    ],
        [ "Model",              model                ],
        [ "Cores",              cores                ],
        [ "",                   hyperthreading       ],
        [ "Clock Frequencies",  frequencies          ],
        [ "Data Cache",         cache                ],
        [ "",                   cachels              ],
        [ "Address Size",       address              ],
        [ "SIMD",               simd                 ],
        [ "Time Stamp Counter", tsc                  ],
        [ "",                   tscinv               ],
        [ "Perf. Monitoring",   perfmon              ],
                                perfmon2...,
                                ibs...,
        [ "Hypervisor",         hypervisor           ],
       ], [:l, :l] ) )
end


end # module CpuId

#=--- end of file ---------------------------------------------------------=#

#=--- CpuId / cpufeature.jl -----------------------------------------------=#

#
# Functions to query single bit feature flags.
#

"""
Tuple of cpuid leaf in eax, result register and bit, and a descriptive string.

This table is an edited combination of sources from [Wikipedia page on
`cpuid`](https://en.wikipedia.org/wiki/CPUID),
[sandpile.org](www.sandpile.org/x86/cpuid.htm), and of course Intel's 4670
page combined [Architectures Software Devleoper Manual](
http://www.intel.com/content/www/us/en/architecture-and-technology/64-ia-32-architectures-software-developer-manual-325462.html).

Expect this table to be incomplete and improvable.
"""
struct CpuFeature
    leaf::UInt32
    subleaf::UInt32
    indx::UInt16
    shft::UInt16
end

const __EAX, __EBX, __ECX, __EDX = map(UInt32, (1, 2, 3, 4))

const SSE3          = CpuFeature( 0x0000_0001, 0x00, __ECX,  0 )
const PCLMUL        = CpuFeature( 0x0000_0001, 0x00, __ECX,  1 )
const DTES64        = CpuFeature( 0x0000_0001, 0x00, __ECX,  2 )
const MON           = CpuFeature( 0x0000_0001, 0x00, __ECX,  3 )
const DSCPL         = CpuFeature( 0x0000_0001, 0x00, __ECX,  4 )
const VMX           = CpuFeature( 0x0000_0001, 0x00, __ECX,  5 )
const SMX           = CpuFeature( 0x0000_0001, 0x00, __ECX,  6 )
const EST           = CpuFeature( 0x0000_0001, 0x00, __ECX,  7 )
const TM2           = CpuFeature( 0x0000_0001, 0x00, __ECX,  8 )
const SSSE3         = CpuFeature( 0x0000_0001, 0x00, __ECX,  9 )
const CNXT          = CpuFeature( 0x0000_0001, 0x00, __ECX, 10 )
const SDBG          = CpuFeature( 0x0000_0001, 0x00, __ECX, 11 )
const FMA3          = CpuFeature( 0x0000_0001, 0x00, __ECX, 12 )
const CX16          = CpuFeature( 0x0000_0001, 0x00, __ECX, 13 )
const XTPR          = CpuFeature( 0x0000_0001, 0x00, __ECX, 14 )
const PDCM          = CpuFeature( 0x0000_0001, 0x00, __ECX, 15 )
const PCID          = CpuFeature( 0x0000_0001, 0x00, __ECX, 17 )
const DCA           = CpuFeature( 0x0000_0001, 0x00, __ECX, 18 )
const SSE41         = CpuFeature( 0x0000_0001, 0x00, __ECX, 19 )
const SSE42         = CpuFeature( 0x0000_0001, 0x00, __ECX, 20 )
const X2APIC        = CpuFeature( 0x0000_0001, 0x00, __ECX, 21 )
const MOVBE         = CpuFeature( 0x0000_0001, 0x00, __ECX, 22 )
const POPCNT        = CpuFeature( 0x0000_0001, 0x00, __ECX, 23 )
const TSCDL         = CpuFeature( 0x0000_0001, 0x00, __ECX, 24 )
const AES           = CpuFeature( 0x0000_0001, 0x00, __ECX, 25 )
const XSAVE         = CpuFeature( 0x0000_0001, 0x00, __ECX, 26 )
const OSXSV         = CpuFeature( 0x0000_0001, 0x00, __ECX, 27 )
const AVX           = CpuFeature( 0x0000_0001, 0x00, __ECX, 28 )
const F16C          = CpuFeature( 0x0000_0001, 0x00, __ECX, 29 )
const RDRND         = CpuFeature( 0x0000_0001, 0x00, __ECX, 30 )
const HYPVS         = CpuFeature( 0x0000_0001, 0x00, __ECX, 31 )
const FPU           = CpuFeature( 0x0000_0001, 0x00, __EDX,  0 )
const VME           = CpuFeature( 0x0000_0001, 0x00, __EDX,  1 )
const DE            = CpuFeature( 0x0000_0001, 0x00, __EDX,  2 )
const PSE           = CpuFeature( 0x0000_0001, 0x00, __EDX,  3 )
const TSC           = CpuFeature( 0x0000_0001, 0x00, __EDX,  4 )
const MSR           = CpuFeature( 0x0000_0001, 0x00, __EDX,  5 )
const PAE           = CpuFeature( 0x0000_0001, 0x00, __EDX,  6 )
const MCE           = CpuFeature( 0x0000_0001, 0x00, __EDX,  7 )
const CX8           = CpuFeature( 0x0000_0001, 0x00, __EDX,  8 )
const APIC          = CpuFeature( 0x0000_0001, 0x00, __EDX,  9 )
const SEP           = CpuFeature( 0x0000_0001, 0x00, __EDX, 11 )
const MTRR          = CpuFeature( 0x0000_0001, 0x00, __EDX, 12 )
const PGE           = CpuFeature( 0x0000_0001, 0x00, __EDX, 13 )
const MCA           = CpuFeature( 0x0000_0001, 0x00, __EDX, 14 )
const CMOV          = CpuFeature( 0x0000_0001, 0x00, __EDX, 15 )
const PAT           = CpuFeature( 0x0000_0001, 0x00, __EDX, 16 )
const PSE36         = CpuFeature( 0x0000_0001, 0x00, __EDX, 17 )
const PSN           = CpuFeature( 0x0000_0001, 0x00, __EDX, 18 )
const CLFSH         = CpuFeature( 0x0000_0001, 0x00, __EDX, 19 )
const DS            = CpuFeature( 0x0000_0001, 0x00, __EDX, 21 )
const ACPI          = CpuFeature( 0x0000_0001, 0x00, __EDX, 22 )
const MMX           = CpuFeature( 0x0000_0001, 0x00, __EDX, 23 )
const FXSR          = CpuFeature( 0x0000_0001, 0x00, __EDX, 24 )
const SSE           = CpuFeature( 0x0000_0001, 0x00, __EDX, 25 )
const SSE2          = CpuFeature( 0x0000_0001, 0x00, __EDX, 26 )
const SS            = CpuFeature( 0x0000_0001, 0x00, __EDX, 27 )
const HTT           = CpuFeature( 0x0000_0001, 0x00, __EDX, 28 )
const TM            = CpuFeature( 0x0000_0001, 0x00, __EDX, 29 )
const IA64          = CpuFeature( 0x0000_0001, 0x00, __EDX, 30 )
const PBE           = CpuFeature( 0x0000_0001, 0x00, __EDX, 31 )
const FSGS          = CpuFeature( 0x0000_0007, 0x00, __EBX,  0 )
const TSCADJ        = CpuFeature( 0x0000_0007, 0x00, __EBX,  1 )
const SGX           = CpuFeature( 0x0000_0007, 0x00, __EBX,  2 )
const BMI1          = CpuFeature( 0x0000_0007, 0x00, __EBX,  3 )
const HLE           = CpuFeature( 0x0000_0007, 0x00, __EBX,  4 )
const AVX2          = CpuFeature( 0x0000_0007, 0x00, __EBX,  5 )
const SMEP          = CpuFeature( 0x0000_0007, 0x00, __EBX,  7 )
const BMI2          = CpuFeature( 0x0000_0007, 0x00, __EBX,  8 )
const ERMS          = CpuFeature( 0x0000_0007, 0x00, __EBX,  9 )
const INVPCID       = CpuFeature( 0x0000_0007, 0x00, __EBX, 10 )
const RTM           = CpuFeature( 0x0000_0007, 0x00, __EBX, 11 )
const PQM           = CpuFeature( 0x0000_0007, 0x00, __EBX, 12 )
const FPDPR         = CpuFeature( 0x0000_0007, 0x00, __EBX, 13 )
const MPX           = CpuFeature( 0x0000_0007, 0x00, __EBX, 14 )
const PQE           = CpuFeature( 0x0000_0007, 0x00, __EBX, 15 )
const AVX512F       = CpuFeature( 0x0000_0007, 0x00, __EBX, 16 )
const AVX512DQ      = CpuFeature( 0x0000_0007, 0x00, __EBX, 17 )
const RDSEED        = CpuFeature( 0x0000_0007, 0x00, __EBX, 18 )
const ADX           = CpuFeature( 0x0000_0007, 0x00, __EBX, 19 )
const SMAP          = CpuFeature( 0x0000_0007, 0x00, __EBX, 20 )
const AVX512IFMA    = CpuFeature( 0x0000_0007, 0x00, __EBX, 21 )
const PCOMMIT       = CpuFeature( 0x0000_0007, 0x00, __EBX, 22 )
const CLFLUSH       = CpuFeature( 0x0000_0007, 0x00, __EBX, 23 )
const CLWB          = CpuFeature( 0x0000_0007, 0x00, __EBX, 24 )
const IPT           = CpuFeature( 0x0000_0007, 0x00, __EBX, 25 )
const AVX512PF      = CpuFeature( 0x0000_0007, 0x00, __EBX, 26 )
const AVX512ER      = CpuFeature( 0x0000_0007, 0x00, __EBX, 27 )
const AVX512CD      = CpuFeature( 0x0000_0007, 0x00, __EBX, 28 )
const SHA           = CpuFeature( 0x0000_0007, 0x00, __EBX, 29 )
const AVX512BW      = CpuFeature( 0x0000_0007, 0x00, __EBX, 30 )
const AVX512VL      = CpuFeature( 0x0000_0007, 0x00, __EBX, 31 )
const PREFTCHWT1    = CpuFeature( 0x0000_0007, 0x00, __ECX,  0 )
const AVX512VBMI    = CpuFeature( 0x0000_0007, 0x00, __ECX,  1 )
const UMIP          = CpuFeature( 0x0000_0007, 0x00, __ECX,  2 )
const PKU           = CpuFeature( 0x0000_0007, 0x00, __ECX,  3 )
const OSPKE         = CpuFeature( 0x0000_0007, 0x00, __ECX,  4 )
const RDPID         = CpuFeature( 0x0000_0007, 0x00, __ECX, 22 )
const SGXLC         = CpuFeature( 0x0000_0007, 0x00, __ECX, 30 )
const AVX512VNNIW   = CpuFeature( 0x0000_0007, 0x00, __EDX,  2 )
const AVX512FMAPS   = CpuFeature( 0x0000_0007, 0x00, __EDX,  3 )
const XSAVEOPT      = CpuFeature( 0x0000_000d, 0x01, __EAX,  0 )
const XSAVEC        = CpuFeature( 0x0000_000d, 0x01, __EAX,  1 )
const XSAVES        = CpuFeature( 0x0000_000d, 0x01, __EAX,  3 )
const AHF64         = CpuFeature( 0x8000_0001, 0x00, __ECX,  0 )
const CMPLEG        = CpuFeature( 0x8000_0001, 0x00, __ECX,  1 )
const SVM           = CpuFeature( 0x8000_0001, 0x00, __ECX,  2 )
const EXTAPIC       = CpuFeature( 0x8000_0001, 0x00, __ECX,  3 )
const CR8D          = CpuFeature( 0x8000_0001, 0x00, __ECX,  4 )
const LZCNT         = CpuFeature( 0x8000_0001, 0x00, __ECX,  5 )
const SSE4A         = CpuFeature( 0x8000_0001, 0x00, __ECX,  6 )
const SSEMISALIGN   = CpuFeature( 0x8000_0001, 0x00, __ECX,  7 )
const PREFETCHW     = CpuFeature( 0x8000_0001, 0x00, __ECX,  8 )
const OSVW          = CpuFeature( 0x8000_0001, 0x00, __ECX,  9 )
const IBS           = CpuFeature( 0x8000_0001, 0x00, __ECX, 10 )
const XOP           = CpuFeature( 0x8000_0001, 0x00, __ECX, 11 )
const SKINIT        = CpuFeature( 0x8000_0001, 0x00, __ECX, 12 )
const WDT           = CpuFeature( 0x8000_0001, 0x00, __ECX, 13 )
const LWP           = CpuFeature( 0x8000_0001, 0x00, __ECX, 15 )
const FMA4          = CpuFeature( 0x8000_0001, 0x00, __ECX, 16 )
const TCE           = CpuFeature( 0x8000_0001, 0x00, __ECX, 17 )
const NODEID        = CpuFeature( 0x8000_0001, 0x00, __ECX, 19 )
const TBM           = CpuFeature( 0x8000_0001, 0x00, __ECX, 21 )
const TOPX          = CpuFeature( 0x8000_0001, 0x00, __ECX, 22 )
const PCXCORE       = CpuFeature( 0x8000_0001, 0x00, __ECX, 23 )
const PCXNB         = CpuFeature( 0x8000_0001, 0x00, __ECX, 24 )
const DBX           = CpuFeature( 0x8000_0001, 0x00, __ECX, 26 )
const PERFTSC       = CpuFeature( 0x8000_0001, 0x00, __ECX, 27 )
const PCXL2I        = CpuFeature( 0x8000_0001, 0x00, __ECX, 28 )
const MONX          = CpuFeature( 0x8000_0001, 0x00, __ECX, 29 )
const FPU_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  0 )
const VME_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  1 )
const DE_           = CpuFeature( 0x8000_0001, 0x00, __EDX,  2 )
const PSE_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  3 )
const TSC_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  4 )
const MSR_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  5 )
const PAE_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  6 )
const MCE_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  7 )
const CX8_          = CpuFeature( 0x8000_0001, 0x00, __EDX,  8 )
const APIC_         = CpuFeature( 0x8000_0001, 0x00, __EDX,  9 )
const SYSCALL       = CpuFeature( 0x8000_0001, 0x00, __EDX, 11 )
const MTRR_         = CpuFeature( 0x8000_0001, 0x00, __EDX, 12 )
const PGE_          = CpuFeature( 0x8000_0001, 0x00, __EDX, 13 )
const MCA_          = CpuFeature( 0x8000_0001, 0x00, __EDX, 14 )
const CMOV_         = CpuFeature( 0x8000_0001, 0x00, __EDX, 15 )
const PAT_          = CpuFeature( 0x8000_0001, 0x00, __EDX, 16 )
const PSE36_        = CpuFeature( 0x8000_0001, 0x00, __EDX, 17 )
const MP            = CpuFeature( 0x8000_0001, 0x00, __EDX, 19 )
const NX            = CpuFeature( 0x8000_0001, 0x00, __EDX, 20 )
const MMXEXT_       = CpuFeature( 0x8000_0001, 0x00, __EDX, 22 )
const MMX_          = CpuFeature( 0x8000_0001, 0x00, __EDX, 23 )
const FXSR_         = CpuFeature( 0x8000_0001, 0x00, __EDX, 24 )
const FFXSR         = CpuFeature( 0x8000_0001, 0x00, __EDX, 25 )
const PG1G          = CpuFeature( 0x8000_0001, 0x00, __EDX, 26 )
const RDTSCP        = CpuFeature( 0x8000_0001, 0x00, __EDX, 27 )
const LM            = CpuFeature( 0x8000_0001, 0x00, __EDX, 29 )
const AMD3DNOWEXT   = CpuFeature( 0x8000_0001, 0x00, __EDX, 30 )
const AMD3DNOW      = CpuFeature( 0x8000_0001, 0x00, __EDX, 31 )
const TS            = CpuFeature( 0x8000_0007, 0x00, __EDX,  0 )
const FID           = CpuFeature( 0x8000_0007, 0x00, __EDX,  1 )
const VID           = CpuFeature( 0x8000_0007, 0x00, __EDX,  2 )
const TTP           = CpuFeature( 0x8000_0007, 0x00, __EDX,  3 )
const TM_           = CpuFeature( 0x8000_0007, 0x00, __EDX,  4 )
const STC           = CpuFeature( 0x8000_0007, 0x00, __EDX,  5 )
const MUL100        = CpuFeature( 0x8000_0007, 0x00, __EDX,  6 )
const HWPS          = CpuFeature( 0x8000_0007, 0x00, __EDX,  7 )
const TSCINV        = CpuFeature( 0x8000_0007, 0x00, __EDX,  8 )
const CPB           = CpuFeature( 0x8000_0007, 0x00, __EDX,  9 )
const EFRO          = CpuFeature( 0x8000_0007, 0x00, __EDX, 10 )
const PFI           = CpuFeature( 0x8000_0007, 0x00, __EDX, 11 )
const PA            = CpuFeature( 0x8000_0007, 0x00, __EDX, 12 )
const FP128         = CpuFeature( 0x8000_001a, 0x00, __EAX,  0 )
const MOVU          = CpuFeature( 0x8000_001a, 0x00, __EAX,  1 )
const FP256         = CpuFeature( 0x8000_001a, 0x00, __EAX,  2 )


const CpuFeatureDescription = Dict{Symbol, String}(
    :SSE3          => "128bit Streaming SIMD Extensions 3",
    :PCLMUL        => "PCLMULQDQ support",
    :DTES64        => "64bit debug store",
    :MON           => "MONITOR and MWAIT instructions",
    :DSCPL         => "CPL qualified debug store",
    :VMX           => "Virtual machine extensions",
    :SMX           => "Safer mode instructions",
    :EST           => "Enhanced SpeedStep",
    :TM2           => "Thermal monitor 2",
    :SSSE3         => "128bit Supplemental Streaming SIMD Extension 3",
    :CNXT          => "L1 context ID",
    :SDBG          => "Silicon debug interface",
    :FMA3          => "Fused multiply-add using three operands",
    :CX16          => "CMPXCHG16B instruction",
    :XTPR          => "disabling sending of task priority messages",
    :PDCM          => "Perfmon and debug capabilities",
    :PCID          => "Process context identifiers",
    :DCA           => "Direct cache access for DMA writes",
    :SSE41         => "128bit Streaming SIMD Extensions 4.1",
    :SSE42         => "128bit Streaming SIMD Extensions 4.2",
    :X2APIC        => "x2APIC support",
    :MOVBE         => "MOVBE instruction",
    :POPCNT        => "POPCNT instruction",
    :TSCDL         => "APIC one-shot operation using TSC deadline value",
    :AES           => "AES encryption instruction set",
    :XSAVE         => "XSAVE, XRESTOR, XSETBV, XGETBV",
    :OSXSV         => "XSAVE enabled by OS",
    :AVX           => "256bit Advanced Vector Extensions, AVX",
    :F16C          => "half-precision float support",
    :RDRND         => "On-chip random number generator",
    :HYPVS         => "Running on hypervisor",
    :FPU           => "Onboard x87 FPU",
    :VME           => "Virtual 8086 mode enhancements",
    :DE            => "Debugging extensions",
    :PSE           => "Page size extensions",
    :TSC           => "Time stamp counter",
    :MSR           => "Model Specific Registers, RDMSR and WRMSR instructions",
    :PAE           => "Physical address extension",
    :MCE           => "Machine check exception",
    :CX8           => "CMPXCHG8 instruction (64bit compare and exchange)",
    :APIC          => "APIC on-chip (Advanced Programmable Interrupt Controller)",
    :SEP           => "SYSENTER and SYSEXIT instructions",
    :MTRR          => "Memory Type Range Registers",
    :PGE           => "Page global bit",
    :MCA           => "Machine Check Architecture (MSR)",
    :CMOV          => "Conditional move CMOV and FCMOV instructions",
    :PAT           => "Page attribute table",
    :PSE36         => "36bit page size extension",
    :PSN           => "Processor serial number (only available on Pentium III)",
    :CLFSH         => "CLFLUSH instruction (SSE2)",
    :DS            => "Debug store to save trace of executed jumps",
    :ACPI          => "Thermal monitor and software controlled clock facilities (MSR)",
    :MMX           => "64bit Multimedia Streaming Extensions",
    :FXSR          => "FXSAVE, FXRSTOR instructions",
    :SSE           => "128bit Streaming SIMD Extensions 1",
    :SSE2          => "128bit Streaming SIMD Extensions 2",
    :SS            => "Self Snoop",
    :HTT           => "Max APIC IDs reserved field is valid",
    :TM            => "Thermal monitor with automatic thermal control",
    :IA64          => "IA64 processor emulating x86",
    :PBE           => "Pending break enable wakeup support",
    :FSGS          => "Access to base of %fs and %gs",
    :TSCADJ        => "IA32_TSC_ADJUST",
    :SGX           => "Software Guard Extensions",
    :BMI1          => "Bit Manipulation Instruction Set 1",
    :HLE           => "Transactional Synchronization Extensions",
    :AVX2          => "SIMD 256bit Advanced Vector Extensions 2",
    :SMEP          => "Supervisor-Mode Execution Prevention",
    :BMI2          => "Bit Manipulation Instruction Set 2",
    :ERMS          => "Enhanced REP MOVSB/STOSB",
    :INVPCID       => "INVPCID instruction",
    :RTM           => "Transactional Synchronization Extensions",
    :PQM           => "Platform Quality of Service Monitoring",
    :FPDPR         => "FPU CS and FPU DS deprecated",
    :MPX           => "Intel MPX (Memory Protection Extensions)",
    :PQE           => "Platform Quality of Service Enforcement",
    :AVX512F       => "AVX-512 Foundation",
    :AVX512DQ      => "AVX-512 Doubleword and Quadword Instructions",
    :RDSEED        => "RDSEED instruction",
    :ADX           => "Intel ADX (Multi-Precision Add-Carry Instruction Extensions)",
    :SMAP          => "Supervisor Mode Access Prevention",
    :AVX512IFMA    => "AVX-512 Integer Fused Multiply-Add Instructions",
    :PCOMMIT       => "PCOMMIT instruction",
    :CLFLUSH       => "CLFLUSHOPT Instructions",
    :CLWB          => "CLWB instruction",
    :IPT           => "Intel Processor Trace",
    :AVX512PF      => "AVX-512 Prefetch Instructions",
    :AVX512ER      => "AVX-512 Exponential and Reciprocal Instructions",
    :AVX512CD      => "AVX-512 Conflict Detection Instructions",
    :SHA           => "Intel SHA extensions",
    :AVX512BW      => "AVX-512 Byte and Word Instructions",
    :AVX512VL      => "AVX-512 Vector Length Extensions",
    :PREFTCHWT1    => "PREFETCHWT1 instruction",
    :AVX512VBMI    => "AVX-512 Vector Bit Manipulation Instructions",
    :UMIP          => "User-mode Instruction Prevention",
    :PKU           => "Memory Protection Keys for User-mode pages",
    :OSPKE         => "PKU enabled by OS",
    :RDPID         => "Read Processor ID",
    :SGXLC         => "SGX Launch Configuration",
    :AVX512VNNIW   => "AVX-512 Neural Network Instructions",
    :AVX512FMAPS   => "AVX-512 Multiply Accumulation Single precision",
    :AHF64         => "LAHF and SAHF in PM64",
    :CMPLEG        => "HTT or CMP flag",
    :SVM           => "VMRUN, VMCALL, VMLOAD and VMSAVE etc.",
    :EXTAPIC       => "Extended APIC space",
    :CR8D          => "MOV from and to CR8D",
    :LZCNT         => "LZCNT instruction",
    :SSE4A         => "Streaming SIMD extensions 4A",
    :SSEMISALIGN   => "Misaligned SSE",
    :PREFETCHW     => "PREFETCHW instruction",
    :OSVW          => "Operating-system-visible workaround",
    :IBS           => "Instruction Based Sampling (IBS)",
    :XOP           => "XOP",
    :SKINIT        => "SKINIT, STGI, DEV",
    :WDT           => "Watch dog timer",
    :LWP           => "LWP",
    :FMA4          => "4-operand fused multiply add instruction",
    :TCE           => "Translation cache extension",
    :NODEID        => "Node ID (MSR)",
    :TBM           => "TBM",
    :TOPX          => "Topology extensions on 0x8000'001D to 0x8000'001E",
    :PCXCORE       => "Core performance counter extensions",
    :PCXNB         => "NB performance counter extensions",
    :DBX           => "Data breakpoint extensions",
    :PERFTSC       => "Performance TSC",
    :PCXL2I        => "L2I performance counter extensions",
    :MONX          => "MONITORX/MWAITX instructions",
    :FPU_          => "FPU",
    :VME_          => "CR4",
    :DE_           => "MOV from and to DR4 and DR5",
    :PSE_          => "PSE",
    :TSC_          => "RSC, RDTSC, CR4.TSC",
    :MSR_          => "RDMSR and WDMSR",
    :PAE_          => "PDPTE 64bit",
    :MCE_          => "MCE (MSR)",
    :CX8_          => "CMPXCHG8B",
    :APIC_         => "APIC",
    :SYSCALL       => "SYSCALL and SYSRET",
    :MTRR_         => "MTRR (MSR)",
    :PGE_          => "PDE and PTE",
    :MCA_          => "MCA",
    :CMOV_         => "CMOVxx",
    :PAT_          => "FCMOVxx",
    :PSE36_        => "4 MB PDE bits 16..13",
    :MP            => "MP capable",
    :NX            => "NXE",
    :MMXEXT_       => "MMX-SSE (AMD specific)",
    :MMX_          => "64bit Streaming SIMD extensions (MMX)",
    :FXSR_         => "AMD and Cyrix specific",
    :FFXSR         => "FFXSR",
    :PG1G          => "PG1G",
    :RDTSCP        => "RDTSCP instruction",
    :LM            => "AMD64 long mode",
    :AMD3DNOWEXT   => "AMD 3DNow! (extended)",
    :AMD3DNOW      => "AMD 3DNow!",
    :TS            => "Temperature sensor",
    :FID           => "Frequency ID control",
    :VID           => "Voltage ID control",
    :TTP           => "Thermal trip",
    :TM_           => "Thermal monitoring",
    :STC           => "Software thermal control",
    :MUL100        => "100 MHz multiplier steps",
    :HWPS          => "Hardware P-state support",
    :TSCINV        => "Invariant TSC",
    :CPB           => "Core performance boost",
    :EFRO          => "Read-only MPERF/APERF",
    :PFI           => "Processor feedback interface",
    :PA            => "Processor accumulator",
    :FP128         => "1x128 bit instead of 2x 64-bit processing",
    :MOVU          => "prefer unaligned MOV over MOVL/MOVH",
    :FP256         => "1x256 bit instead of 2x128-bit processing",
)


"""
    cpufeature( feature::Symbol ) ::Bool
    cpufeature( feature::CpuFeature ) ::Bool

Query the CPU whether it supports the given feature.  For fast checking
provide directly the `CpuFeature` defined as a global const in `CpuId`.
Explicitly typed `CpuFeature`s got by the same name as the corresponding
symbols.  Valid symbols are available from `keys(CpuId.CpuFeatureDescription)`.
"""
function cpufeature end

function cpufeature(feature::CpuFeature) ::Bool
    @inbounds exx = getindex( cpuid(feature.leaf, feature.subleaf), feature.indx )
    (exx >> feature.shft) & 0x01 != 0x00
end

# Convenience overload to use symbol notatation
cpufeature(feat::Symbol) = cpufeature( getfield(CpuId, feat) )


"""
    cpufeaturedesc( feature::Symbol ) ::String

Get the textual description of a CPU feature flag given as a *symbol*.
"""
cpufeaturedesc(feature::Symbol) = get( CpuFeatureDescription, feature
                                     , "Unknown feature, no description available!" )


"""
    cpufeatures() ::Vector{Symbol}

Get a list of symbols of all cpu supported features.  Might be extensive and
not exactly useful other than for testing purposes.  Also, this implementation
is not efficient since each feature is queried independently.
"""
cpufeatures() = Symbol[f for f in keys(CpuFeatureDescription) if cpufeature(getfield(CpuId, f))] |> sort


"""
    cpufeaturetable() ::MarkdownString

Generate a markdown table of all the detected/available/supported CPU features
along with some textural description.
"""
function cpufeaturetable() ::MarkdownString
    tbl = MarkdownTable([["Cpu Feature", "Description"]], [:l, :l])
    for f in cpufeatures()
        push!(tbl.rows, [string(f), cpufeaturedesc(f)])
    end
    MarkdownString(tbl)
end


#=--- end of file ---------------------------------------------------------=#

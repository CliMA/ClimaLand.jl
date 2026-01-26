# Development Notes

# AMD

Source: [Preliminary Processor Programming Reference (PPR) for AMD Family 17h Models 00h-0Fh Processors](https://support.amd.com/TechDocs/54945_PPR_Family_17h_Models_00h-0Fh.pdf)

Modern AMD uses terminology 'Processor' > 'Node' (aka 'Core Complex') > 'Core' > 'Thread'.
'Logical Processor Count' refers to the sum (product) of all.
This is known as Simultaneous Multithreading, SMT, aka on Intel as Hyperthreading.

Requirement:
Core::X86::Cpuid::FeatureExtIdEcx[TopologyExtensions] != 0

Todo: Check 'CmpLegacy' and 'HTT'

CPUID_Fn80000001_ECX, Core::X86::Cpuid::FeatureExtIdEcx, Bit 22, TopologyExtensions
topology extensions support, 1=Indicates support for Core::X86::Cpuid::CachePropEax0

CPUID_Fn00000001_EBX, Bits 23:16, Core::X86::Cpuid::SizeId[NC]
Specifies the number of threads in the processor as 
Core::X86::Cpuid::SizeId[NC]+1.

CPUID_Fn8000001E_EBX, Bits 15:8, ThreadsPerCore
The number of threads per core is ThreadsPerCore+1.

CPUID_Fn8000001E_ECX, Bits 10:8, NodesPerProcessor
Node per processor
Valid Values:
Value     | Description
000b      | 1 node per processor.
001b      | 2 nodes per processor.
010b      | Reserved.
011b      | 4 nodes per processor.
111b-100b | Reserved.


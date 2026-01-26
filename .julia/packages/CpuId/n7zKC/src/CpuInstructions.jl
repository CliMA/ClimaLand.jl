#=--- CpuId / CpuInstructions.jl ------------------------------------------=#

"""
# Module 'CpuInstructions'

The module 'CpuInstructions' is part of the package 'CpuId', and provides a
selection of wrapped low-level assembly functions to diagnose potential
computational efficiency issues.

Though primarily intended as a helper module to 'CpuId', the functions may be
used directly in other code e.g. for benchmarking purposes.  Just include the
file directly, or copy & paste.
"""
module CpuInstructions

export cpuid, rdtsc, rdtscp

using Base: llvmcall

"""
    cpuid( [leaf], [subleaf]) ::NTuple{4, UInt32}

Invoke the cpu's hardware instruction `cpuid` with the values of the arguments
stored as registers *EAX = leaf*, *ECX = subleaf*, respectively. Returns a
tuple of the response of registers EAX, EBX, ECX, EDX.  Input values may be
given as individual `UInt32` arguments, or converted from any `Integer`.
Unspecified arguments are assumed zero.

This function is primarily intended as a low-level interface to the CPU.

Note: Expected to work on all CPUs that implement the assembly instruction
      `cpuid`, which is at least Intel and AMD.
"""
function cpuid end

# Convenience function allowing passing other than UInt32 values
function cpuid(leaf=0, subleaf=0)
    # for some reason, we need a dedicated local
    # variable of UInt32 for llvmcall to succeed
    l, s = UInt32(leaf), UInt32(subleaf)
    cpuid_llvm(l, s) ::NTuple{4,UInt32}
end

#
#   TODO:
#   The following llvmcall routines fail when being inlined!
#   Hence the @noinline.
#

#
# Test Sys.ARCH for valid CPU architectures at compile time
#
@static if Sys.ARCH in (:x86, :x86_64)

    # Low level cpuid call, taking eax=leaf and ecx=subleaf,
    # returning eax, ebx, ecx, edx as NTuple(4,UInt32)
    @noinline cpuid_llvm(leaf::UInt32, subleaf::UInt32) =
        llvmcall("""
            ; leaf = %0, subleaf = %1, %2 is some label
            ; call 'cpuid' with arguments loaded into registers EAX = leaf, ECX = subleaf
            %3 = tail call { i32, i32, i32, i32 } asm sideeffect "cpuid",
                "={ax},={bx},={cx},={dx},{ax},{cx},~{dirflag},~{fpsr},~{flags}"
                (i32 %0, i32 %1) #2
            ; retrieve the result values and convert to vector [4 x i32]
            %4 = extractvalue { i32, i32, i32, i32 } %3, 0
            %5 = extractvalue { i32, i32, i32, i32 } %3, 1
            %6 = extractvalue { i32, i32, i32, i32 } %3, 2
            %7 = extractvalue { i32, i32, i32, i32 } %3, 3
            ; return the values as a new tuple
            %8  = insertvalue [4 x i32] undef, i32 %4, 0
            %9  = insertvalue [4 x i32]   %8 , i32 %5, 1
            %10 = insertvalue [4 x i32]   %9 , i32 %6, 2
            %11 = insertvalue [4 x i32]  %10 , i32 %7, 3
            ret [4 x i32] %11
        """
        # llvmcall requires actual types, rather than the usual (...) tuple
        , NTuple{4,UInt32}, Tuple{UInt32,UInt32}
        , leaf, subleaf)

    @inline rdtsc() =
        llvmcall("""
            %1 = tail call { i32, i32 } asm sideeffect "rdtsc", "={ax},={dx},~{dirflag},~{fpsr},~{flags}"() #2
            %2 = extractvalue { i32, i32 } %1, 0
            %3 = extractvalue { i32, i32 } %1, 1
            %4 = zext i32 %2 to i64
            %5 = zext i32 %3 to i64
            %6 = shl nuw i64 %5, 32
            %7 = or i64 %6, %4
            ret i64 %7
        """
        , UInt64, Tuple{})


    @noinline rdtscp() =
        llvmcall("""
            %1 = tail call { i32, i32, i32 } asm sideeffect "rdtscp", "={ax},={dx},={cx},~{dirflag},~{fpsr},~{flags}"() #2
            %2 = extractvalue { i32, i32, i32 } %1, 0
            %3 = extractvalue { i32, i32, i32 } %1, 1
            %4 = zext i32 %2 to i64
            %5 = zext i32 %3 to i64
            %6 = shl nuw i64 %5, 32
            %7 = or i64 %6, %4

            %8 = extractvalue { i32, i32, i32 } %1, 2
            %9 = zext i32 %8 to i64

            %10 = insertvalue [2 x i64] undef, i64  %7, 0
            %11 = insertvalue [2 x i64]  %10 , i64  %9, 1
            ret [2 x i64] %11
        """
        , Tuple{UInt64,UInt64}, Tuple{})

else  # Sys.ARCH  other than (:x86, :x86_64)

    #
    # Create fallback functions to avoid failing (pre-)compilations
    #
    @noinline cpuid_llvm(::UInt32, ::UInt32) =
                    (zero(UInt32), zero(UInt32), zero(UInt32), zero(UInt32))
    @inline   rdtsc()  = zero(UInt64)
    @noinline rdtscp() = (zero(UInt64), zero(UInt64))

end  # Sys.ARCH


end # module CpuInstructions

#=--- end of file ---------------------------------------------------------=#

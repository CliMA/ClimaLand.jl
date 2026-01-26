import LLVM
import InteractiveUtils

using UnsafeAtomics: UnsafeAtomics, acquire, release, acq_rel, seq_cst
using UnsafeAtomics.Internal: OP_RMW_TABLE, inttypes
using Test

llvmptr(xs::Array, i) = reinterpret(Core.LLVMPtr{eltype(xs),0}, pointer(xs, i))

function check_default_ordering(T::Type)
    xs = T[rand(T), rand(T)]
    x1 = rand(T)
    x2 = rand(T)
    check_default_ordering(xs, x1, x2)
end

function check_default_ordering(xs::AbstractArray{T}, x1::T, x2::T) where T
    @debug "xs=$(repr(xs)) x1=$(repr(x1)) x2=$(repr(x2))"

    ptr = llvmptr(xs, 1)
    GC.@preserve xs begin
        @test UnsafeAtomics.load(ptr) === xs[1]
        UnsafeAtomics.store!(ptr, x1)
        @test xs[1] === x1
        sizeof(T) == 0 && return # CAS hangs on zero sized data...
        desired = (old = x1, success = true)
        @test UnsafeAtomics.cas!(ptr, x1, x2) === (old = x1, success = true)
        @test xs[1] === x2
        @testset for (op, name) in OP_RMW_TABLE
            xs[1] = x1
            @test UnsafeAtomics.modify!(ptr, op, x2) === (x1 => op(x1, x2))
            @test xs[1] === op(x1, x2)

            rmw = getfield(UnsafeAtomics, Symbol(name, :!))
            xs[1] = x1
            @test rmw(ptr, x2) === x1
            @test xs[1] === op(x1, x2)

            # Check dispatch to LLVM atomic OP instead of CAS loop.
            if (op == +) || (op == -)
                IR = sprint(io->InteractiveUtils.code_llvm(io,
                    UnsafeAtomics.modify!,
                    typeof.((ptr, +, T(1)))))
                @test occursin("atomicrmw", IR)
            end
        end
    end
end

function test_explicit_ordering(T::Type = UInt)
    xs = T[rand(T), rand(T)]
    x1 = rand(T)
    x2 = rand(T)
    test_explicit_ordering(xs, x1, x2)
end

function test_explicit_ordering(xs::AbstractArray{T}, x1::T, x2::T) where T
    @debug "xs=$(repr(xs)) x1=$(repr(x1)) x2=$(repr(x2))"

    ptr = llvmptr(xs, 1)
    GC.@preserve xs begin

        @test UnsafeAtomics.load(ptr, acquire) === xs[1]
        UnsafeAtomics.store!(ptr, x1, release)
        @test xs[1] === x1
        sizeof(T) == 0 && return # CAS hangs on zero sized data...
        desired = (old = x1, success = true)
        @test UnsafeAtomics.cas!(ptr, x1, x2, acq_rel, acquire) === desired
        @test xs[1] === x2
        @testset for (op, name) in OP_RMW_TABLE
            xs[1] = x1
            @test UnsafeAtomics.modify!(ptr, op, x2, acq_rel) === (x1 => op(x1, x2))
            @test xs[1] === op(x1, x2)

            rmw = getfield(UnsafeAtomics, Symbol(name, :!))
            xs[1] = x1
            @test rmw(ptr, x2, acquire) === x1
            @test xs[1] === op(x1, x2)

            # Test syncscopes.
            if (op == +) || (op == -)
                xs[1] = x1
                @test UnsafeAtomics.modify!(ptr, op, x2, seq_cst, UnsafeAtomics.none) ===
                      (x1 => op(x1, x2))
                @test xs[1] === op(x1, x2)

                xs[1] = x1
                @test UnsafeAtomics.modify!(ptr, op, x2, seq_cst, UnsafeAtomics.singlethread) ===
                      (x1 => op(x1, x2))
                @test xs[1] === op(x1, x2)
            end

            # Check dispatch to LLVM atomic OP instead of CAS loop.
            if (op == +) || (op == -)
                IR = sprint(io->InteractiveUtils.code_llvm(io,
                    UnsafeAtomics.modify!,
                    typeof.((ptr, +, T(1), seq_cst, UnsafeAtomics.singlethread))))
                @test occursin("atomicrmw", IR)
            end
        end
    end
end


@testset "UnsafeAtomicsLLVM" begin
    @testset for T in inttypes
        check_default_ordering(T)
        test_explicit_ordering(T)
    end

    @testset "Zero-sized types" begin
        @test sizeof(Nothing) == 0
        check_default_ordering([nothing, nothing], nothing, nothing)
        test_explicit_ordering([nothing, nothing], nothing, nothing)
    end
end


for (ityp,jtyp) ∈ [("i8", UInt8), ("i16", UInt16), ("i32", UInt32), ("i64", UInt64), ("i128", UInt128)]
    @eval begin
        @inline function _atomic_load(ptr::Ptr{$jtyp})
            Base.llvmcall($(
                @static if VERSION >= v"1.12-DEV"
                    # use opaque pointers #53
                    """
                    %v = load atomic $(ityp), ptr %0 acquire, align $(Base.gc_alignment(jtyp))
                    ret $(ityp) %v
                    """
                else
                    """
                    %p = inttoptr i$(8sizeof(Int)) %0 to $(ityp)*
                    %v = load atomic $(ityp), $(ityp)* %p acquire, align $(Base.gc_alignment(jtyp))
                    ret $(ityp) %v
                    """
                end
            ), $jtyp, Tuple{Ptr{$jtyp}}, ptr)
        end
        @inline function _atomic_store!(ptr::Ptr{$jtyp}, x::$jtyp)
            Base.llvmcall($(
                @static if VERSION >= v"1.12-DEV"
                    """
                    store atomic $(ityp) %1, ptr %0 release, align $(Base.gc_alignment(jtyp))
                    ret void
                    """
                else
                    """
                    %p = inttoptr i$(8sizeof(Int)) %0 to $(ityp)*
                    store atomic $(ityp) %1, $(ityp)* %p release, align $(Base.gc_alignment(jtyp))
                    ret void
                    """
                end
            ), Cvoid, Tuple{Ptr{$jtyp}, $jtyp}, ptr, x)
        end
        @inline function _atomic_cas_cmp!(ptr::Ptr{$jtyp}, cmp::$jtyp, newval::$jtyp)
            Base.llvmcall($(
                @static if VERSION >= v"1.12-DEV"
                    """
                    %c = cmpxchg ptr %0, $(ityp) %1, $(ityp) %2 acq_rel acquire
                    %bit = extractvalue { $ityp, i1 } %c, 1
                    %bool = zext i1 %bit to i8
                    ret i8 %bool
                    """
                else
                    """
                    %p = inttoptr i$(8sizeof(Int)) %0 to $(ityp)*
                    %c = cmpxchg $(ityp)* %p, $(ityp) %1, $(ityp) %2 acq_rel acquire
                    %bit = extractvalue { $ityp, i1 } %c, 1
                    %bool = zext i1 %bit to i8
                    ret i8 %bool
                    """
                end
            ), Bool, Tuple{Ptr{$jtyp}, $jtyp, $jtyp}, ptr, cmp, newval)
        end
    end
end

"""
Operations supported by `atomicrmw`.

- Julia v1.11 use [libLLVM-16](https://releases.llvm.org/16.0.0/docs/LangRef.html#atomicrmw-instruction)
"""
const atomicrmw_ops = [
    "xchg", "add", "sub",
    "and", "nand", "or", "xor",
    "max", "min", "umax", "umin",
    # "fadd", "fsub", 
    # "fmax", "fmin",
    # "uinc_wrap", "udec_wrap",  # Need LLVM 16: VERSION >= v"1.11"
]
for op ∈ atomicrmw_ops
    f = Symbol("_atomic_", op, '!')
    for (ityp,jtyp) ∈ [("i8", UInt8), ("i16", UInt16), ("i32", UInt32), ("i64", UInt64), ("i128", UInt128)]
        @eval begin
            # Do inplace `$(op)` for `*ptr` and `x` atomically, return the old value of `*ptr`.
            @inline function $f(ptr::Ptr{$jtyp}, x::$jtyp)
                Base.llvmcall($(
                    @static if VERSION >= v"1.12-DEV"
                        """
                        %v = atomicrmw $(op) ptr %0, $(ityp) %1 acq_rel
                        ret $(ityp) %v
                        """
                    else
                        """
                        %p = inttoptr i$(8sizeof(Int)) %0 to $(ityp)*
                        %v = atomicrmw $op $(ityp)* %p, $(ityp) %1 acq_rel
                        ret $(ityp) %v
                        """
                    end
                ), $jtyp, Tuple{Ptr{$jtyp}, $jtyp}, ptr, x)
            end
        end
    end
    @eval begin
        @inline function $f(ptr::Ptr{UInt}, x::ThreadState)
            reinterpret(ThreadState, $f(reinterpret(Ptr{UInt32}, ptr), reinterpret(UInt32, x)))
        end
    end
end

@inline _atomic_state(ptr::Ptr{UInt}) = reinterpret(ThreadState, _atomic_load(reinterpret(Ptr{UInt32}, ptr)))
@inline _atomic_store!(ptr::Ptr{UInt}, x::ThreadState) = _atomic_store!(reinterpret(Ptr{UInt32}, ptr), reinterpret(UInt32, x))
@inline function _atomic_cas_cmp!(ptr::Ptr{UInt}, cmp::ThreadState, newval::ThreadState)
    _atomic_cas_cmp!(reinterpret(Ptr{UInt32}, ptr), reinterpret(UInt32, cmp), reinterpret(UInt32, newval))
end

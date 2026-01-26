mutable struct SixelAllocator
    ptr::Ptr{C.sixel_allocator_t}
    function SixelAllocator()
        # libsixel checks if the pointer is NULL during initialization so we use a non-zero
        # "invalid" address value to bypass the NULL check
        ptr_ref = Ref{Ptr{C.sixel_allocator_t}}(1)
        status = C.sixel_allocator_new(ptr_ref, C_NULL, C_NULL, C_NULL, C_NULL)
        check_status(status)
        new(ptr_ref[])
    end
end

mutable struct SixelOutput
    ptr::Ptr{C.sixel_output_t}
    function SixelOutput(fn_write, priv; allocator::SixelAllocator=SixelAllocator())
        # libsixel checks if the pointer is NULL during initialization so we use a non-zero
        # "invalid" address value to bypass the NULL check
        ptr_ref = Ref{Ptr{C.sixel_output_t}}(1)
        status = C.sixel_output_new(ptr_ref, fn_write, priv, allocator.ptr)
        check_status(status)
        new(ptr_ref[])
    end
end

mutable struct SixelDither
    ptr::Ptr{C.sixel_dither_t}
    function SixelDither(
            data, width, height, pixelformat, quality_mode,
            method_for_largest=C.LARGE_NORM,
            method_for_rep=C.REP_CENTER_BOX;
            allocator::SixelAllocator=SixelAllocator())
        # libsixel checks if the pointer is NULL during initialization so we use a non-zero
        # "invalid" address value to bypass the NULL check
        ptr_ref = Ref{Ptr{C.sixel_dither_t}}(1)
        ncolors = -1 # libsixel internally sets -1 to length(palette)
        status = C.sixel_dither_new(ptr_ref, ncolors, allocator.ptr)
        check_status(status)
        ptr = ptr_ref[]
        C.sixel_dither_initialize(ptr,
            data,
            Int32(width),
            Int32(height),
            Int32(pixelformat),
            Int32(method_for_largest),
            Int32(method_for_rep),
            Int32(quality_mode)
        )
        new(ptr)
    end
end

check_status(status) = status == 0 || error("libsixel error code: $status")

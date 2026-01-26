"""
    read!(io, arr, comp)

Read in an array `arr` from the [`TiffFile`](@ref) or [`TiffFileStrip`](@ref)
stream `io`, inflating the data using compression method `comp`. `read!` will
dispatch on the value of compression and use the correct compression technique
to read the data.
"""
Base.read!(tfs::Union{TiffFile, TiffFileStrip}, arr::AbstractVector{UInt8}, comp::CompressionType) = read!(tfs.io, arr, Val(comp))

Base.read!(tfs::Union{TiffFile, TiffFileStrip}, arr::AbstractVector{UInt8}, ::Val{COMPRESSION_NONE}) = read!(tfs.io, arr)

function Base.read!(tfs::TiffFileStrip, arr::AbstractVector{UInt8}, ::Val{COMPRESSION_PACKBITS})
    len = length(arr)

    pos = 1
    nbit = Array{Int8}(undef, 1)
    nxt = Array{UInt8}(undef, 1)
    while pos < len
        read!(tfs.io, nbit)
        n = nbit[1]
        if 0 <= n <= 127
            read!(tfs.io, view(arr, pos:pos+n))
            pos += n
        elseif -127 <= n <= -1
            read!(tfs.io, nxt)
            arr[pos:(pos-n)] .= nxt[1]
            pos += -n
        end
        pos += 1
    end
end

function Base.read!(tfs::TiffFileStrip, arr::AbstractVector{UInt8}, ::Val{COMPRESSION_DEFLATE})
    readbytes!(InflateZlibStream(tfs.io), arr)
end

function Base.read!(tfs::TiffFileStrip, arr::AbstractVector{UInt8}, ::Val{COMPRESSION_ADOBE_DEFLATE})
    readbytes!(InflateZlibStream(tfs.io), arr)
end

function lzw_decode!(io, arr)
    CLEAR_CODE::Int = 256 + 1
    EOI_CODE::Int = 257 + 1
    TABLE_ENTRY_LENGTH_BITS::Int = 16
    TABLE_ENTRY_OFFSET_BITS::Int = 8 * sizeof(Int) - TABLE_ENTRY_LENGTH_BITS

    output_size = length(arr)

    GC.@preserve arr begin
        out_pointer::Ptr{UInt8} = pointer(arr)
        out_position::Int = 0 # current position in out

        table_size::Int = output_size * 2 + 258
        table_pointer::Ptr{UInt8} = reinterpret(Ptr{UInt8}, Libc.malloc(table_size)) # table of strings
        table_offsets_pointer::Ptr{Int} = reinterpret(Ptr{Int}, Libc.malloc(sizeof(Int) * 4097)) # offsets into table

        @inline create_table_entry(length, offset) = Base.shl_int(length, TABLE_ENTRY_OFFSET_BITS) | offset
        @inline table_entry_length(table_entry) = Base.lshr_int(table_entry, TABLE_ENTRY_OFFSET_BITS)
        @inline table_entry_offset(table_entry) = table_entry & (Base.shl_int(1, TABLE_ENTRY_OFFSET_BITS) - 1)

        try
            # InitializeTable();
            foreach(i -> unsafe_store!(table_pointer + i, UInt8(i)), 0:255)
            foreach(i -> unsafe_store!(table_offsets_pointer, create_table_entry(1, i), i+1), 0:259) # length is stored in upper 16 bits

            code = -1

            buffer::Int=0 # buffer for reading in codes
            bitcount::Int=0 # number of valid bits in buffer
            codesize::Int=9 # current number of bits per code
            input::Vector{UInt8} = Vector{UInt8}(undef, bytesavailable(io))
            read!(io, input)
            function getcode(buffer, code, bitcount, codesize, i)
                old_code::Int = code

                # make sure we have enough bits in the buffer
                if bitcount < codesize
                    buffer = Base.shl_int(buffer, 8) | input[i+=1]
                    bitcount += 8
                end

                # one more time (since the max code size is 12 bits, only need to check twice)
                if bitcount < codesize
                    buffer = Base.shl_int(buffer, 8) | input[i+=1]
                    bitcount += 8
                end

                code = Base.lshr_int(buffer, bitcount - codesize) & (Base.shl_int(1, codesize) - 1)
                bitcount -= codesize
                # code + 1 because this is Julia
                (buffer, code + 1, old_code, bitcount, codesize, i)
            end

            @inline check_table_overflow(start, length) = start + length > table_size && error("LZW: table buffer overflow")
            @inline check_output_overflow(start, length) = start + length > output_size && error("LZW: output buffer overflow")

            # annotated with excerpts from the LZW pseudocode in the TIFF 6.0 spec
            # https://developer.adobe.com/content/dam/udp/en/open/standards/tiff/TIFF6.pdf
            table_count::Int = 258 # number of (valid) table entries; 256 one-byte codes + CLEAR_CODE + EOI_CODE
            next_table_offset::Int = 258
            input_pos::Int = 0 # current position in input
            while true
                # GetNextCode()
                (buffer, code, old_code, bitcount, codesize, input_pos) = getcode(buffer, code, bitcount, codesize, input_pos)
                if code == EOI_CODE || out_position >= output_size
                    break
                elseif code == CLEAR_CODE # reset table
                    # InitializeTable();
                    table_count = 258
                    next_table_offset = 258
                    codesize = 9
                    # Code = GetNextCode();
                    (buffer, code, old_code, bitcount, codesize, input_pos) = getcode(buffer, code, bitcount, codesize, input_pos)
                    if code == EOI_CODE
                        break
                    end
                    # WriteString(StringFromCode(Code))
                    r = unsafe_load(table_offsets_pointer, code)
                    len = table_entry_length(r)

                    check_output_overflow(out_position, len)

                    memcpy(out_pointer + out_position, table_pointer + table_entry_offset(r), len)
                    out_position += len
                else
                    if code <= table_count
                        # WriteString(StringFromCode(Code));
                        if code <= 256
                            # this is redundant with the check above, but it makes
                            # the code easier to reason about and less bug prone
                            check_output_overflow(out_position, 1)

                            unsafe_store!(out_pointer + out_position, code - 1)
                            out_position += 1
                        else
                            r = unsafe_load(table_offsets_pointer, code)
                            len = table_entry_length(r)

                            check_output_overflow(out_position, len)

                            memcpy(out_pointer + out_position, table_pointer + table_entry_offset(r), len)
                            out_position += len
                        end

                        # AddStringToTable(StringFromCode(OldCode) + FirstChar(StringFromCode(Code)));
                        table_count += 1
                        len = 1
                        if old_code <= 256
                            check_table_overflow(next_table_offset, 2) # this byte + the next one

                            unsafe_store!(table_pointer + next_table_offset, UInt8(old_code - 1))
                        else
                            r = unsafe_load(table_offsets_pointer, old_code)
                            len = table_entry_length(r)

                            check_table_overflow(next_table_offset, len + 1) # these bytes + the next one

                            memcpy(table_pointer + next_table_offset, table_pointer + table_entry_offset(r), len)
                        end

                        if code <= 256
                            unsafe_store!(table_pointer + next_table_offset + len, UInt8(code - 1))
                        else
                            r = unsafe_load(table_offsets_pointer, code)
                            memcpy(table_pointer + next_table_offset + len, table_pointer + table_entry_offset(r), 1)
                        end
                        unsafe_store!(table_offsets_pointer, create_table_entry(len + 1, next_table_offset), table_count)
                        next_table_offset += len + 1
                    else
                        # WriteString(StringFromCode(OldCode) + FirstChar(StringFromCode(OldCode)));
                        r = unsafe_load(table_offsets_pointer, old_code)
                        len = table_entry_length(r)

                        check_output_overflow(out_position, len + 1)

                        memcpy(out_pointer + out_position, table_pointer + table_entry_offset(r), len)
                        unsafe_store!(out_pointer + out_position + len, unsafe_load(table_pointer + table_entry_offset(r)))
                        out_position += len + 1

                        check_table_overflow(next_table_offset, len + 1)

                        # AddStringToTable(StringFromCode(OldCode) + FirstChar(StringFromCode(OldCode)));
                        table_count += 1
                        memcpy(table_pointer + next_table_offset, table_pointer + table_entry_offset(r), len)
                        memcpy(table_pointer + next_table_offset + len, table_pointer + table_entry_offset(r), 1)
                        unsafe_store!(table_offsets_pointer, create_table_entry(len + 1, next_table_offset), table_count)
                        next_table_offset += len + 1
                    end
                end

                if out_position < output_size
                    if table_count == 511
                        codesize = 10
                    elseif table_count == 1023
                        codesize = 11
                    elseif table_count == 2047
                        codesize = 12
                    end
                end
            end

            out_position != output_size && @warn "LZW: expected $output_size bytes, got $out_position bytes"
            out_position == output_size && code != EOI_CODE && @warn "LZW: missing EOI code"
        catch e
            error("LZW: $e")
            rethrow()
        finally
            Libc.free(table_pointer)
            Libc.free(table_offsets_pointer)
        end
    end
end

function Base.read!(tfs::TiffFileStrip{S}, arr::AbstractVector{UInt8}, ::Val{COMPRESSION_LZW}) where S
    lzw_decode!(tfs, arr)
end

function Base.read!(::Union{TiffFile, TiffFileStrip}, ::AbstractVector{UInt8}, ::Val{C}) where C
    error("Compression $C is not implemented. Please open an issue against TiffImages.jl.")
end

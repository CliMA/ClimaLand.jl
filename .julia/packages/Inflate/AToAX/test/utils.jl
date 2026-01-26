function all_zeros_deflate(min_length)
    n = max(0, cld(min_length - 249, 1032))
    deflate_start = [0xed, 0xc1, 0x01, 0x01, 0x00, 0x00, 0x00, 0x82,
                     0x20, 0xff, 0xaf, 0x6e, 0x48, 0x40]
    deflate_end = [0xaf, 0x06]
    return vcat(deflate_start, zeros(UInt8, n), deflate_end), 249 + n * 1032
end

function minimum_gzip_header()
    return [0x1f, 0x8b,              # Gzip ID bytes.
            0x08,                    # Compression method (deflate).
            0x00,                    # Flags
            0x00, 0x00, 0x00, 0x00,  # MTIME
            0x00,                    # Extra flags
            0xff]                    # OS (unknown)
end

function all_zeros_crc(n)
    c = Inflate.init_crc()
    for i in 1:n
        c = Inflate.update_crc(c, 0x00)
    end
    return Inflate.finish_crc(c)
end

function all_zeros_gzip(min_length)
    deflate, n = all_zeros_deflate(min_length)
    crc = reinterpret(UInt8, [all_zeros_crc(n)])
    len = reinterpret(UInt8, [n % UInt32])
    return vcat(minimum_gzip_header(), deflate, crc, len)
end

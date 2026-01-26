import .Iterators.partition

"""
    $(TYPEDEF)

An image file directory is a sorted collection of the tags representing this
plane in the TIFF file. They behave like dictionaries except that tags aren't
required to be unique, so given an IFD called `ifd`, we can add new tags as
follows:

```jldoctest; setup = :(ifd = TiffImages.IFD(UInt32))
julia> ifd[TiffImages.IMAGEDESCRIPTION] = "Some details";

julia> ifd[TiffImages.IMAGEWIDTH] = 512;

julia> ifd
IFD, with tags:
	Tag(IMAGEWIDTH, 512)
	Tag(IMAGEDESCRIPTION, "Some details")
```

!!! note
    Tags are not required to be unique! See [`TiffImages.Iterable`](@ref) for
    how to work with duplicate tags.
"""
struct IFD{O <: Unsigned}
    tags::DefaultDict{UInt16, Vector{Tag}, DataType}
end

IFD(::Type{O}) where {O <: Unsigned} = IFD{O}(DefaultDict{UInt16, Vector{Tag}}(Vector{Tag}))
IFD(::Type{O}, tags) where {O <: Unsigned} = IFD{O}(tags)

"""
A wrapper to force getindex to return the underlying array instead of only the
first element. Usually the first element is sufficient, but sometimes access to
the array is needed (to add duplicate entries or access them).

```jldoctest; setup = :(ifd = TiffImages.IFD(UInt32))
julia> using TiffImages: Iterable

julia> ifd[TiffImages.IMAGEDESCRIPTION] = "test"
"test"

julia> ifd[Iterable(TiffImages.IMAGEDESCRIPTION)] # since wrapped with Iterable, returns array
1-element Vector{TiffImages.Tag}:
 Tag(IMAGEDESCRIPTION, "test")

julia> ifd[Iterable(TiffImages.IMAGEDESCRIPTION)] = "test2" # since wrapped with Iterable, it appends
"test2"

julia> ifd
IFD, with tags:
	Tag(IMAGEDESCRIPTION, "test")
	Tag(IMAGEDESCRIPTION, "test2")

```
"""
struct Iterable{T}
    key::T
end

Base.length(ifd::IFD) = sum(map(length, values(ifd.tags)))
Base.keys(ifd::IFD) = keys(ifd.tags)
Base.values(ifd::IFD) = values(ifd.tags)
Base.iterate(ifd::IFD) = iterate(ifd.tags)
Base.iterate(ifd::IFD, n::Int) = iterate(ifd.tags, n)

Base.getindex(ifd::IFD, key::Iterable{TiffTag}) = getindex(ifd, Iterable(UInt16(key.key)))
Base.getindex(ifd::IFD, key::Iterable{UInt16}) = getindex(ifd.tags, key.key)
Base.getindex(ifd::IFD, key::TiffTag) = getindex(ifd, UInt16(key))
Base.getindex(ifd::IFD, key::UInt16) = first(getindex(ifd.tags, key))

getdata(f::F, ifd::IFD, key::TiffTag, default) where F = getdata(f, ifd, UInt16(key), default)
function getdata(f::F, ifd::IFD, key::UInt16, default) where F
    val = get(ifd.tags, key, nothing)
    val === nothing && return default
    return f(first(val).data)
end
getdata(ifd::IFD, key, default) = getdata(identity, ifd, key, default)

Base.in(key::TiffTag, v::IFD) = in(UInt16(key), v)
Base.in(key::UInt16, v::IFD) = in(key, keys(v))
Base.delete!(ifd::IFD, key::TiffTag) = delete!(ifd, UInt16(key))
Base.delete!(ifd::IFD, key::UInt16) = delete!(ifd.tags, key)

Base.similar(::IFD{O}) where {O <: Unsigned} = IFD(O)
Base.merge(ifd::IFD{O}, other::IFD) where {O <: Unsigned} = IFD(O, DefaultDict(Vector{Int}, merge(ifd.tags, other.tags)))

Base.setindex!(ifd::IFD, value::Tag, key::UInt16) = setindex!(ifd.tags, [value], key)
Base.setindex!(ifd::IFD, value::Tag, key::TiffTag) = setindex!(ifd, value, UInt16(key))

Base.setindex!(ifd::IFD, value, key::TiffTag) = setindex!(ifd, value, UInt16(key))
Base.setindex!(ifd::IFD, value, key::UInt16) = setindex!(ifd, Tag(key, value), key)

Base.setindex!(ifd::IFD, value, key::Iterable{TiffTag}) = setindex!(ifd, value, Iterable(UInt16(key.key)))
Base.setindex!(ifd::IFD, value, key::Iterable{UInt16}) = setindex!(ifd, Tag(key.key, value), key)
Base.setindex!(ifd::IFD, value::Tag, key::Iterable{UInt16}) = push!(ifd.tags[key.key], value)

function isloaded(ifd::IFD)
    for tags in values(ifd.tags)
        for tag in tags
            (!isloaded(tag)) && return false
        end
    end
    true
end

"""
    sizeof(ifd)

Number of bytes that an IFD will use on disk.
"""
function Base.sizeof(ifd::IFD{O}) where {O}
    sz = O == UInt32 ? 2 : sizeof(O)
    for tags in values(ifd)
        for tag in tags
            # tag, data, length, and data
            sz += 2 + 2 + sizeof(O) + sizeof(O)
            if sizeof(tag) > sizeof(O) # if data in tag is larger than slot
                sz += sizeof(tag) # we have add all the additional bytes
            end
        end
    end
    sz += sizeof(O) # slot for subsequent IFD locations
end
"""
    $SIGNATURES

Checks if the data in this IFD is contiguous on disk. Striped data can be read
faster as one contiguous chunk if possible.
"""
function iscontiguous(ifd::IFD)
    if !(ROWSPERSTRIP in ifd) || nrows(ifd) <= ifd[ROWSPERSTRIP].data
        return true
    else
        return all(diff(ifd[STRIPOFFSETS].data) .== ifd[STRIPBYTECOUNTS].data[1:end-1])
    end
end

"""
    $(SIGNATURES)

Updates an [`TiffImages.IFD`](@ref) by replacing all instances of the
placeholder type [`TiffImages.RemoteData`](@ref) with the actual data from the
file `tf`.
"""
function load!(tf::TiffFile, ifd::IFD)
    for key in sort(collect(keys(ifd)))
        tags = ifd[Iterable(key)]
        for (idx, tag) in enumerate(tags)
            tags[idx] = load(tf, tag)
        end
    end
end

function Base.show(io::IO, ifd::IFD)
    print(io, "IFD, with tags: ")
    for key in sort(collect(keys(ifd)))
        tags = ifd[Iterable(key)]
        for tag in tags
            print(io, "\n\t", tag)
        end
    end
end

function Base.read(tf::TiffFile{O}, ::Type{IFD}) where O <: Unsigned
    # Regular TIFF's use 16bits instead of 32 bits for entry data
    N = O == UInt32 ? read(tf, UInt16) : read(tf, O)

    entries = IFD(O)

    for i in 1:N
        tag = read(tf, Tag)
        push!(entries[Iterable(tag.tag)], tag)
    end

    next_ifd = Int(read(tf, O))
    entries, next_ifd
end

function Base.iterate(file::TiffFile{O}) where {O}
    seek(file.io, file.first_offset)
    iterate(file, (read(file, IFD)))
end

"""
    iterate(file, state) -> IFD, Int

Advances the iterator to the next IFD.

**Output**
- `Vector{Int}`: Offsets within file for all strips corresponding to the current
   IFD
- `Int`: Offset of the next IFD
"""
function Base.iterate(file::TiffFile, state::Tuple{Union{<: IFD, Nothing}, Int})
    curr_ifd, next_ifd_offset = state
    # if current element doesn't exist, exit
    (curr_ifd === nothing) && return nothing
    (next_ifd_offset <= 0) && return (curr_ifd, (nothing, 0))

    seek(file.io, next_ifd_offset)
    next_ifd, next_ifd_offset = read(file, IFD)

    return (curr_ifd, (next_ifd, next_ifd_offset))
end

"""
    $(TYPEDEF)

A strip is a contiguous block of separately-encoded image data. A TIFF
file will typically have multiple strips, each decoding to multiple rows of
pixels in the image

$(FIELDS)
"""
struct TiffFileStrip{O}
    """Strip data"""
    io::IOBuffer

    """The IFD corresponding to this strip"""
    ifd::IFD{O}
end

Base.read!(tfs::TiffFileStrip, arr::AbstractArray) = read!(tfs.io, arr)
Base.bytesavailable(tfs::TiffFileStrip) = bytesavailable(tfs.io)

function Base.read!(target::AbstractArray{T, N}, tf::TiffFile{O, S}, ifd::IFD{O}) where {T, N, O, S}
    offsets = TILEOFFSETS in ifd ? ifd[TILEOFFSETS].data : ifd[STRIPOFFSETS].data
    compression = getdata(CompressionType, ifd, COMPRESSION, COMPRESSION_NONE)

    rtype = rawtype(ifd)
    spp = ispalette(ifd) ? 1 : nsamples(ifd)
    rows = nrows(ifd)
    cols = ncols(ifd)
    bps = bitspersample(ifd)
    samples = reinterpret(rtype, target)

    function info()
        compression_name = "?"
        if compression == COMPRESSION_LZW
            compression_name = "LZW"
        elseif compression == COMPRESSION_DEFLATE || compression == COMPRESSION_ADOBE_DEFLATE
            compression_name = "Zip"
        elseif compression == COMPRESSION_PACKBITS
            compression_name = "PackBits"
        end

        cmprssn = compression == COMPRESSION_NONE ? "uncompressed" : "$compression_name-compressed"
        chunks = (istiled(ifd) ? "tile" : "strip") * (length(offsets) == 1 ? "" : "s")
        "reading $(cmprssn), $(isplanar(ifd) ? "planar" : "chunky") image with $(length(offsets)) $(chunks)"
    end

    @debug info()

    # number of input bytes in each strip or tile
    encoded_bytes = TILEBYTECOUNTS in ifd ? ifd[TILEBYTECOUNTS].data : ifd[STRIPBYTECOUNTS].data

    # number of samples (pixels * channels) in each strip or tile
    if istiled(ifd)
        @debug "tile size: $(tilecols(ifd)) x $(tilerows(ifd))"
        # tiled images are always padded to tile boundaries
        strip_samples = map(_ -> tilecols(ifd) * tilerows(ifd) * (isplanar(ifd) ? 1 : spp), encoded_bytes)
    else
        nstrips = length(encoded_bytes)
        rowsperstrip = getdata(Int, ifd, ROWSPERSTRIP, rows)

        if isplanar(ifd)
            # planar files will have separate strips or tiles for each channel
            temp = fill(rowsperstrip * cols, cld(rows, rowsperstrip))
            temp[end] = (rows - (rowsperstrip * (length(temp) - 1))) * cols
            strip_samples = repeat(temp, spp)
        else
            strip_samples = fill(rowsperstrip * cols * spp, nstrips)
            strip_samples[end] = (rows - (rowsperstrip * (nstrips - 1))) * cols * spp
        end

        @assert sum(strip_samples) == rows * cols * spp
    end

    parallel_enabled = something(tryparse(Bool, get(ENV, "JULIA_IMAGES_PARALLEL", "1")), false)
    do_parallel = parallel_enabled && rows * cols > 250_000 # pixels

    comp = Val(compression)
    if is_complicated(ifd)
        tasks::Vector{Task} = []
        start = 1
        for (offset::Int, len::Int, bytes::Int) in zip(offsets, strip_samples, encoded_bytes)
            @debug "reading strip with $len samples from $bytes encoded bytes"

            seek(tf, offset)
            arr = view(samples, start:(start+len-1))
            data = Vector{UInt8}(undef, bytes)
            read!(tf, data)
            tfs = TiffFileStrip{O}(IOBuffer(data), ifd)

            function go(tfs, arr, comp)
                cls = istiled(ifd) ? tilecols(ifd) : cols
                cls = isplanar(ifd) ? cls : cls * spp # number of samples (not pixels) per column
                rws = fld(length(arr), cls)
                sz = uncompressed_size(ifd, cls, rws)
                read!(tfs, view(reinterpret(UInt8, vec(arr)), 1:sz), comp)
                if is_irregular_bps(ifd)
                    arr .= unpack_integers(arr, rws, cls, bps)
                end
                reverse_prediction!(tfs.ifd, arr)
            end

            if do_parallel
                push!(tasks, Threads.@spawn go(tfs, arr, comp))
            else
                go(tfs, arr, comp)
            end

            start += len
        end

        for task in tasks
            wait(task)
        end
    else
        @debug "fast path for uncomplicated images"

        seek(tf, first(offsets))
        read!(tf, target)
    end

    if isplanar(ifd)
        samplesv = vec(samples)
        deplane!(samplesv, spp)
    end
end

function Base.write(tf::TiffFile{O}, ifd::IFD{O}) where {O <: Unsigned}
    if !isloaded(ifd)
        error("Cannot write unloaded IFDs. Use `load!` to populate tags with remote data")
    end

    N = length(ifd)
    O == UInt32 ? write(tf, UInt16(N)) : write(tf, UInt64(N))

    # keep track of which tags are too large to fit in the IFD slot and need a
    # remote location for their data
    remotedata = Vector{Pair{Tag, Vector{Int}}}()
    sorted_keys = sort(collect(keys(ifd)))
    for k in sorted_keys
        tags = ifd[Iterable(k)]
        for tag in tags
            pos = position(tf.io)
            if !write(tf, tag)
                push!(remotedata, tag => [pos])
            end
        end
    end

    # end position, write a zero by default, but this should be updated if any
    # more IFDs are written
    ifd_end_pos = position(tf.io)
    write(tf, O(0))

    for (tag, poses) in remotedata
        tag = tag::Tag
        data_pos = position(tf.io)
        data = tag.data
        # add NUL terminator to the end of Strings that don't have it already
        if eltype(tag) === String
            data = data::SubString{String}
            if !endswith(data, '\0')
                data *= '\0'
            end
            write(tf, data)  # compile-time dispatch
        else
            write(tf, data)  # run-time dispatch
        end
        push!(poses, data_pos)
    end

    for (tag, poses) in remotedata
        tag = tag::Tag
        orig_pos, data_pos = poses
        seek(tf, orig_pos)
        write(tf, tag, data_pos)
    end

    seek(tf, ifd_end_pos)

    return ifd_end_pos
end

# reverse any pre-processing that might have been applied to sample
# values prior to compression
# https://www.awaresystems.be/imaging/tiff/tifftags/predictor.html
function reverse_prediction!(ifd::IFD, arr::AbstractArray{T,N}) where {T, N}
    pred::Int = predictor(ifd)
    # for planar data, each "pixel" in the strip is actually a single channel
    spp::Int = isplanar(ifd) ? 1 : nsamples(ifd)
    columns = istiled(ifd) ? tilecols(ifd) : ncols(ifd)
    rows = istiled(ifd) ? tilerows(ifd) : cld(length(arr), columns * spp)

    GC.@preserve arr begin
        if pred == 2
            @debug "reversing horizontal differencing"

            temp::Ptr{T} = pointer(arr)
            for row in 1:rows
                start = (row - 1) * columns * spp
                for plane in 1:spp
                    previous::T = unsafe_load(temp, start + plane)
                    for i in (spp + plane):spp:(columns - 1) * spp + plane
                        current = unsafe_load(temp, start + i) + previous
                        unsafe_store!(temp, current, start + i)
                        previous = current
                    end
                end
            end
        elseif pred == 3
            @debug "reversing float differencing"

            columns = columns * spp * sizeof(T)

            buffer::Vector{UInt8} = Vector{UInt8}(undef, columns)

            valn = Val(sizeof(T))

            temp2::Ptr{UInt8} = pointer(reinterpret(UInt8, arr))

            for row in 1:rows
                start = (row - 1) * columns
                for plane in 1:spp
                    prev::UInt8 = unsafe_load(temp2, start + plane)
                    for i in (spp + plane):spp:(columns - 1) + plane
                        curr = unsafe_load(temp2, start + i) + prev
                        unsafe_store!(temp2, curr, start + i)
                        prev = curr
                    end
                end
                vw = view(reinterpret(UInt8, arr), start+1:start+columns)
                deplane!(buffer, vw, valn)
            end
        end
    end
end

# {AAA...BBB...CCC...} => {ABCABCABC...}
function deplane!(arr::AbstractVector{T}, n::Integer) where T
    out = Vector{T}(undef, length(arr))
    deplane!(out, arr, Val(n))
end

const is_mac_or_windows_x64 = (Sys.iswindows() || Sys.isapple()) && Sys.ARCH == :x86_64

# {AAA...BBB...CCC...} => {ABCABCABC...}
function deplane!(buffer::AbstractVector{T}, arr::AbstractVector{T}, n::Val{N}) where {T, N}
    @assert length(buffer) == length(arr)
    @assert length(arr) % N == 0

    GC.@preserve arr buffer begin
        if Int(pointer(arr)) & 0x3f > 0 || length(arr) < 64 || is_mac_or_windows_x64
            # small or not 64-byte aligned
            temp = deplane_slow(arr, N)
            GC.@preserve temp begin
                memcpy(pointer(arr), pointer(temp), sizeof(temp))
            end
        else
            deplane_simd!(buffer, arr, n)
            memcpy(pointer(arr), pointer(buffer), sizeof(buffer))
        end
    end
end

# {AAA...BBB...CCC...} => {ABCABCABC...}
function deplane_slow(arr::AbstractVector{T}, n) where T
    @debug "rearranging planar data"
    reshape(arr, fld(length(arr), n), n)'[:]
end

# {AAA...BBB...CCC...} => {ABCABCABC...}
@generated function deplane_simd!(out::Vector{T}, arr::AbstractVector{T}, ::Val{N}) where {T, N}
    width = cld(sizeof(T) * N, 64) * 64
    count = fld(width, sizeof(T) * N) # pixels per iteration

    sym(x) = Symbol("q", x)

    # vload {AAA...}
    # vload {BBB...}
    # ...
    loads = map(x -> :($(sym(x + 1)) = vload(Vec{$count * $(max(1, x)), T}, ptrA + (index + num_pixels * $x - 1) * $(sizeof(T)))), 0:N-1)

    # shuffle1 = {0, X, 1, X+1, ...}
    # shuffle2 = {0, 1, X, 2, 3, X+1, ...}
    # ...
    perms = []
    for k in 1:N-1
        left = count * k # number of elements in the first shuffle vector (see below)
        # take `k` elements from the first vector, followed by one from the second vector, repeated `count` times
        shuffle::Vector{Int} = mapreduce(x -> vcat(x...), vcat, zip(collect.(partition(0:left-1,k)), left:left+count-1))
        push!(perms, :($(Symbol("shuffle", k)) = $(Val(Tuple(shuffle)))))
    end

    # shufflevector {AAA...}, {BBB...}, shuffle1 => {ABABAB...}
    # shufflevector {ABABAB...}, {CCC...}, shuffle2 => {ABCABCABC...}
    # ...
    shuffles::Vector{Expr} = []
    ll, rr = sym(1), 1
    for x in 1:N-1
        input1 = ll
        ll = sym(x + N)
        input2 = sym(rr += 1)
        push!(shuffles, :($ll = shufflevector($input1, $input2, $(Symbol("shuffle", x)))))
    end
    push!(shuffles, :(final = $ll))

    # assignments for each channel in the final loop
    finish = map(x -> :(out[start + i * $N + $x] = arr[iterations * $count + i + num_pixels * $x + 1]), 0:N-1)

    quote
        @debug "rearranging planar data (SIMD)"

        GC.@preserve arr begin
            ptrA = pointer(arr)
            num_pixels = fld(length(arr), N)
            iterations = fld(num_pixels, $count) - 1
            out_index = 1 # output index

            $(perms...)

            @inbounds for index in 1:$count:iterations*$count
                $(loads...)
                $(shuffles...)

                vstore(final, out, out_index)

                out_index += $count * N
            end

            remaining = num_pixels - iterations * $count
            start = iterations * $count * N + 1

            @inbounds for i in 0:remaining-1
                $(finish...)
            end
        end
    end
end

unpack_integers(v::AbstractVector, n::Integer) = unpack_integers(v, 1, length(v), n)

unpack_integers(v::AbstractVector, c::Integer, n::Integer) = unpack_integers(v, fld(length(v), c), c, n)

# unpack `r` rows of `c` n-bit integers to useful word-sized integers
unpack_integers(v::AbstractVector, r, c, n::Integer) = unpack_integers(v, r, c, Val(n))

# for SIMD, must have (c % M) == 0 && (M * N) % (K * 8) == 0, where M is
# the vector width used by the SIMD algorithm; M = 32 doesn't work for
# {27, 29, 31} (K = 8), so we step up to the next power of two
for N in (27, 29, 31)
    @eval unpack_integers(v::AbstractVector, r, c, n::Val{$N}) = c % 64 == 0 ? unpack_integers_simd(v, n) : unpack_integers_slow(v, r, c, $N)
end

unpack_integers(v::AbstractVector, r, c, n::Val{N}) where N = c % 32 == 0 ? unpack_integers_simd(v, n) : unpack_integers_slow(v, r, c, N)

# {AAA, ABB, BBC, CCC} => {AAAA, BBBB, CCCC}
function unpack_integers_slow(v::AbstractVector{T}, rows::Integer, columns::Integer, n::Integer) where T
    @debug "recoding from $n bits per sample"

    vb::Vector{UInt8} = reinterpret(UInt8, vec(v))
    out::Vector{T} = Vector{T}(undef, length(v))
    i = 0 # input index
    j = 0 # output index
    # encoding is done per row, so decoding is also done per row
    for _ in 1:rows
        buffer::UInt64 = 0
        available = 0 # number of valid bits available in buffer
        for _ in 1:columns
            while available < n
                buffer = (buffer << 8) | vb[i+=1]
                available += 8
            end
            val = (buffer >> (available - n)) & ((UInt64(1) << n) - 1)
            out[j+=1] = T(val)
            available -= n
        end
    end
    out
end

# these are "nice" n, for which all packed n-bit integers are encoded within no more
# than sizeof(T) bytes, where T is the smallest integer type big enough to store n bits
const nice_n = filter(x -> !(max(nextpow(2, x), 8) - x in [0,1,2,3,5]), 1:31)

# {AAA, ABB, BBC, CCC} => {AAAA, BBBB, CCCC}
@generated function unpack_integers_simd(A::AbstractVector{T}, n::Val{N}) where {T, N}
    nice = N in nice_n
    TT = nice ? T : widen(T)
    width = max(32, sizeof(TT) * 8) # SIMD vector width
    m = fld(width * N, sizeof(TT) * 8) # input bytes per vector
    count = fld(width, sizeof(TT)) # number of codes per vector

    lp(v::AbstractVector,n) = vcat(zeros(Int, n - length(v)), v)

    function shuffle(N)
        bitrange = 0:m * 8 - 1
        extents = extrema.(partition(bitrange, N))
        Val(Tuple(mapreduce(y -> lp(collect(UnitRange(map(x -> fld.(x, 8), y)...)), sizeof(TT)), vcat, extents)))
    end

    function shift(N)
        len = N * count
        Vec(7 .- last.(extrema.(Iterators.partition(0:len-1,N))) .% 8...)
    end

    mask(N) = Vec(fill(TT(2^N - 1), count)...)

    function main_block(i)
        sym = Symbol
        quote
            $(sym("a", i)) = vload(Vec{$m, UInt8}, in_ptr)
            # arrange bytes so that each TT-sized lane contains a single code (+ extra bits)
            $(sym("b", i)) = shufflevector($(sym("a", i)), shuffle)
            # shift out extra low-order bits
            $(sym("c", i)) = _bswap(reinterpret(Vec{$count, $TT}, $(sym("b", i)))) >> shift
            # mask out extra high-order bits
            $(sym("d", i)) = $(sym("c", i)) & mask
            in_ptr += $m
        end
    end

    function nice_block()
        quote
            $(main_block(1))

            vstore(d1, out_ptr)
            out_ptr += $width
        end
    end

    function non_nice_block()
        quote
            $(main_block(1)) # gives d1
            $(main_block(2)) # gives d2

            # d1 and d2 come from the main_blocks
            t1 = reinterpret(Vec{$count * 2, T}, d1)
            t2 = reinterpret(Vec{$count * 2, T}, d2)

            # shufflevector {AXBXCX...}, {DXEXFX...} => {ABC...DEF...}
            f = shufflevector(t1, t2, shuffle2)

            vstore(f, out_ptr)
            out_ptr += $width
        end
    end

    quote
        @debug "recoding from $N bits per sample (SIMD, $($nice ? "nice" : "not nice"))"

        # decoded integers per iteration
        num = fld($width, sizeof(T))

        @assert length(A) % num == 0

        out = Vector{T}(undef, length(A))
        GC.@preserve A out begin
            in_ptr::Ptr{UInt8} = reinterpret(Ptr{UInt8}, pointer(A))
            out_ptr::Ptr{T} = pointer(out)

            shuffle = $(shuffle(N))
            shift = $(shift(N))
            mask = $(mask(N))

            shuffle2 = $(Val(Tuple(0:2:count*4-1)))

            iterations = fld(length(A), num)
            for i in 1:iterations
                $(nice ? nice_block() : non_nice_block())
            end

            out
        end
    end
end
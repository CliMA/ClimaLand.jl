"""
    $(TYPEDEF)

A type to represent lazily-loaded TIFF data, returned by `TiffImages.load(filepath; lazyio=true)`.
Useful for opening and operating on images too large to store in memory,
and for incrementally writing new TIFF files.

This works by buffering individual slices. This allows broad format support, including compressed
TIFFs, but with mixed performance depending on your specific access (indexing) patterns.
See discussion in the package documentation, and [`MmappedTIFF`](@ref) for an alternative
with different strengths and weaknesses.

```jldoctest
julia> using TiffImages, ColorTypes

julia> img = empty(LazyBufferedTIFF, Gray{Float32}, joinpath(mktempdir(), "test.tif"))
32-bit LazyBufferedTIFF{Gray{Float32}} 0×0×0 (writable)
    Current file size on disk:   8 bytes
    Addressable space remaining: 4.000 GiB

julia> close(img) # when you're done, close the stream
```

$(FIELDS)
"""
mutable struct LazyBufferedTIFF{T <: ColorOrTuple, O <: Unsigned, AA <: AbstractArray} <: AbstractDenseTIFF{T, 3}

    """
    Pointer to keep track of the backing file
    """
    file::TiffFile{O}
    """
    The associated tags for each slice in this array
    """
    ifds::Vector{IFD{O}}

    dims::NTuple{3, Int}

    """
    An internal cache to fill reading from disk
    """
    working_cache::AA

    """
    An internal cache to fill with fully transformed data
    """
    cache::AA

    """
    The index of the currently loaded slice
    """
    cache_index::Int

    """
    Position of last loaded IFD, updated whenever a slice is appended
    """
    last_ifd_offset::O

    """
    A flag tracking whether this file is editable
    """
    readonly::Bool

    function LazyBufferedTIFF(file::TiffFile{O}, ifds::Vector{IFD{O}}, dims, cache::AA, cache_index::Int, last_ifd_offset::O, readonly::Bool) where {O, AA <: AbstractArray}
        new{eltype(cache), O, typeof(cache)}(file, ifds, dims, cache, cache, cache_index, last_ifd_offset, readonly)
    end
end

function LazyBufferedTIFF(file::TiffFile{O}, ifds::Vector{IFD{O}}) where {O}
    ifd = ifds[1]
    dims = (nrows(ifd), ncols(ifd), length(ifds))
    cache = getcache(ifd)
    LazyBufferedTIFF(file, ifds, dims, cache, -1, zero(O), true)
end

"""
    empty(LazyBufferedTIFF, T, filepath; bigtiff)

Create a new memory-mapped file ready with element type `T` for appending future
slices. The `bigtiff` flag, if true, allows 64-bit offsets for data larger than
~4GB.

```jldoctest; setup=:(rm("test.tif", force=true))
julia> using ColorTypes, FixedPointNumbers # for Gray{N0f8} type

julia> img = empty(LazyBufferedTIFF, Gray{N0f8}, "test.tif"); # make memory-mapped image

julia> push!(img, rand(Gray{N0f8}, 100, 100));

julia> push!(img, rand(Gray{N0f8}, 100, 100));

julia> size(img)
(100, 100, 2)
```
"""
function Base.empty(::Type{LazyBufferedTIFF}, ::Type{T}, filepath; bigtiff=false) where {T <: ColorOrTuple}
    if isfile(filepath)
        error("This file already exists, please use `TiffImages.load` to open")
    end
    LazyBufferedTIFF(T, getstream(format"TIFF", open(filepath, "w+"), filepath); bigtiff = bigtiff)
end

function memmap(t::Type{T}, filepath; bigtiff=false) where {T <: ColorOrTuple}
    Base.depwarn("`memmap` is deprecated, please use empty(LazyBufferedTIFF, $t, $filepath; bigtiff = $bigtiff)", :memmap, force = true)
    empty(LazyBufferedTIFF, t, filepath; bigtiff = bigtiff)
end

function LazyBufferedTIFF(::Type{T}, io::Stream; bigtiff = false) where {T}
    O = bigtiff ? UInt64 : UInt32
    tf = TiffFile{O}(io)

    last_ifd_offset = write(tf) # write out header

    LazyBufferedTIFF(tf, IFD{O}[], (0, 0, 0), Array{T}(undef, 1, 1), -1, O(last_ifd_offset), false)
end

Base.size(A::LazyBufferedTIFF) = A.dims
offset(::LazyBufferedTIFF{T, O, AA}) where {T, O, AA} = O

function Base.getindex(A::LazyBufferedTIFF{T, O, AA}, i1::Int, i2::Int, i::Int) where {T, O, AA}
    (size(A) == (0, 0, 0)) && error("This image has not been initialized, please `push!` data into it first")
    # check the loaded cache is already the correct slice
    if A.cache_index != i
        ifd = ifds(A)[i]

        # if the file isn't open, lets open a handle and update it
        if !isopen(A.file.io)
            path = A.file.filepath
            A.file.io = getstream(format"TIFF", open(path), path)
        end

        read!(A.working_cache, A.file, ifd)

        A.cache = transform(A.working_cache, ifd)

        A.cache_index = i
    end

    return A.cache[i1, i2]
end

function Base.setindex!(A::LazyBufferedTIFF, I...)
    error("Unable to mutate inplace since this array is on disk. Convert to a mutable in-memory version by running "*
          "`copy(arr)`. \n\n𝗡𝗼𝘁𝗲: For large files this can be quite expensive. A future PR will add "*
          "support for writing inplace to disk. See `push!` for appending to an array.")
end

"""
    push!(img::LazyBufferedTIFF, slice::AbstractMatrix)

Push a 2D slice to a memory-mapped file. The slice must be the same `eltype` as the
target `img` and the `img` must be not be readonly.
"""
function Base.push!(A::LazyBufferedTIFF{T, O, AA}, data::AbstractMatrix{T}) where {T, O, AA}
    ifd = _constructifd(data, offset(A.file))

    push!(A, DenseTaggedImage(data, ifd))
end

"""
    push!(img::LazyBufferedTIFF, tiff)

Push a 2D slice and associated IFD to a memory-mapped file. The slice must be
the same `eltype` as the target `img`. This allows for writing custom tags to
the file

```jldoctest; setup = :(using TiffImages, ColorTypes, FixedPointNumbers)
julia> slice = TiffImages.DenseTaggedImage(Gray.(zeros(UInt8, 10, 10)));

julia> ifds(slice)[TiffImages.IMAGEDESCRIPTION] = "Custom info";

julia> temp_file = joinpath(mktempdir(), "tmp.tif");

julia> lazy_tiff = empty(LazyBufferedTIFF, Gray{N0f8}, temp_file; bigtiff = true);

julia> push!(lazy_tiff, slice);

julia> close(lazy_tiff) # done writing, close open stream

julia> ifds(TiffImages.load(temp_file))[TiffImages.IMAGEDESCRIPTION] # read from disk
Tag(IMAGEDESCRIPTION, "Custom info")

```
"""
function Base.push!(A::LazyBufferedTIFF{T, O, AA}, tiff::D) where {T, O, AA, D <: AbstractTIFF{T, 2}}
    (A.readonly) && error("This image is read only")
    data = tiff.data

    # merge a minimal IFD with the correct offset with any user provided tags
    ifd = merge(_constructifd(data, offset(A.file)), ifds(tiff))

    if size(A) == (0, 0, 0) # if this is the initial slice pushed, initialize the size
        A.dims = (size(data)..., 0)
    end

    @assert size(data) == (A.dims[1], A.dims[2]) "Pushed slices must have dimensions: $((A.dims[1], A.dims[2]))"

    pagecache = Vector{UInt8}(undef, size(data, 2) * sizeof(T) * size(data, 1))

    if A.last_ifd_offset + sizeof(pagecache) + sizeof(ifd) > typemax(O)
        @info "No more room in file @ N = $(length(ifds(A)))"
    end
    seekend(A.file.io)
    prev_ifd_slice = _writeslice(pagecache, A.file, data, ifd, A.last_ifd_offset)

    # update data
    A.dims = (A.dims[1], A.dims[2], A.dims[3] + 1)
    push!(ifds(A), ifd)
    A.last_ifd_offset = prev_ifd_slice
    A.cache = data
    A.cache_index = A.dims[3]
    A
end

Base.push!(A::DenseTaggedImage{T, N, O, AA}, data) where {T, N, O, AA <: LazyBufferedTIFF} = push!(A.data, data)

function Base.show(io::IO, ::MIME"text/plain", A::LazyBufferedTIFF{T, O, AA}) where {T, O, AA}
    printstyled(io, O == UInt32 ? "32-bit" : "64-bit"; color = :cyan)
    print(io, " LazyBufferedTIFF{$(T)} ")
    printstyled(io, "$(size(A, 1))×$(size(A, 2))×$(size(A, 3))"; bold=true)
    if A.readonly
        printstyled(io, " (readonly)"; color=:red)
    else
        printstyled(io, " (writable)"; color=:green)
    end
    println(io)
    ondisk = sizeof(A.file)
    ondisk += (size(A) == (0, 0, 0)) ? 0 : sum(sizeof.(ifds(A))) + sizeof(T) * reduce(*, size(A))
    println(io, "    Current file size on disk:   $(Base.format_bytes(ondisk))")
    println(io, "    Addressable space remaining: $(Base.format_bytes(typemax(O) - ondisk))")
end

Base.close(tiff::L) where {L <: LazyBufferedTIFF} = close(tiff.file.io)
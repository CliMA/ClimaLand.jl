module TiffImages

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@compiler_options"))
    @eval Base.Experimental.@compiler_options optimize=1
end

using ColorTypes
using DocStringExtensions
using FileIO
using FixedPointNumbers
using IndirectArrays
using OffsetArrays
using DataStructures
using PkgVersion
using ProgressMeter
using PrecompileTools: @setup_workload, @compile_workload
using Base.Iterators
using Inflate
using UUIDs
using Mmap
using SIMD

const PKGVERSION = @PkgVersion.Version 0

include("enum.jl")
include(joinpath("types", "widepixel.jl"))
include("utils.jl")
include("files.jl")
include("tags.jl")
include("ifds.jl")
include("compression.jl")
include("layout.jl")
include(joinpath("types", "common.jl"))
include(joinpath("types", "dense.jl"))
include(joinpath("types", "strided.jl"))
include(joinpath("types", "lazy.jl"))
include(joinpath("types", "mmapped.jl"))
include("load.jl")

export memmap, LazyBufferedTIFF, ifds, color, nchannels, channel

@deprecate TiffFile(::Type{O}) where O<:Unsigned TiffFile{O}()

@setup_workload begin
    files = filter(x -> endswith(x, ".tif"), readdir(joinpath(@__DIR__, "precomp-tifs", "images"), join = true))
    @compile_workload begin
        for file in files
            TiffImages.load(file)
            TiffImages.load(file; mmap=true)
            TiffImages.load(file; lazyio=true)
        end
    end
end

end # module

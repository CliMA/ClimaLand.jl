using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Inflate
using BenchmarkTools
using CodecZlib
using Graphs
using StatsBase
using Random
using Printf

# Generate an incompressible vector of length `n`.
incompressible(n) = rand(UInt8, n)

# Generate a Huffman compressible vector of length `n`.
function huffman_compressible(n)
    w = 0.844 .^ (0:255)
    return sample(StatsBase.shuffle(UInt8.(0:255)), Weights(w), n)
end

# Generate a runlength compressible vector of length `n`.
function runlength(n)
    x = UInt8[]
    while length(x) < n
        m = min(sample(1:32), n - length(x))
        append!(x, fill(rand(UInt8), m))
    end
    return x
end

# Generate a Graphs save file of length approximately `n`.
# Note: the balance between nodes and edges needs to be adjusted for
# savefiles smaller than about 620 bytes.
function graph(n)
    buf = IOBuffer()
    m = ceil(Int, n / (2 * (log10(n) - 1)))
    local x
    for k = 1:4
        g = erdos_renyi(m ÷ 10, m, is_directed = true)
        savegraph(buf, g, "", LGFormat())
        x = take!(buf)
        m = m * n ÷ length(x)
    end
    return x
end

sizes = Dict(:small => 1_000, :medium => 100_000, :large => 10_000_000)
generators = Dict(:incompressible => incompressible,
                  :huffman => huffman_compressible,
                  :runlength => runlength,
                  :graph => graph)
compressors = Dict(:deflate => DeflateCompressor,
                   :zlib => ZlibCompressor,
                   :gzip => GzipCompressor)
zlib_decompressors = Dict(:deflate => DeflateDecompressor,
                          :zlib => ZlibDecompressor,
                          :gzip => GzipDecompressor)
inflate_decompressors = Dict(:deflate => inflate,
                             :zlib => inflate_zlib,
                             :gzip => inflate_gzip)
zlib_decompressor_streams = Dict(:deflate => DeflateDecompressorStream,
                                 :zlib => ZlibDecompressorStream,
                                 :gzip => GzipDecompressorStream)
inflate_decompressor_streams = Dict(:deflate => InflateStream,
                                    :zlib => InflateZlibStream,
                                    :gzip => InflateGzipStream)

function run_tests()
    results = Dict{Any, Float64}()
    for data_size in keys(sizes)
        for data_type in keys(generators)
            Random.seed!(13)
            x = generators[data_type](sizes[data_size])
            for format in keys(compressors)
                compressed = transcode(compressors[format], x)
                zlib_decompressor = zlib_decompressors[format]
                inflate_decompressor = inflate_decompressors[format]
                GC.gc()
                t1 = @belapsed transcode($zlib_decompressor, $compressed)
                GC.gc()
                t2 = @belapsed $inflate_decompressor($compressed)
                results[(data_size, data_type, format, :zlib, :in_memory)] = t1
                results[(data_size, data_type, format, :inflate, :in_memory)] = t2
                zlib_decompressor_stream = zlib_decompressor_streams[format]
                inflate_decompressor_stream = inflate_decompressor_streams[format]
                GC.gc()
                t3 = @belapsed read($zlib_decompressor_stream(stream)) setup=(stream = IOBuffer($compressed)) evals=1
                GC.gc()
                t4 = @belapsed read($inflate_decompressor_stream(stream)) setup=(stream = IOBuffer($compressed)) evals=1
                results[(data_size, data_type, format, :zlib, :streaming)] = t3
                results[(data_size, data_type, format, :inflate, :streaming)] = t4
            end
        end
    end
    return results
end

function print_results(results, mode)
    mode_string = replace(string(mode), "_" => " ")
    println("Results for $(mode_string):")
    for data_type in [:incompressible, :huffman, :runlength, :graph]
        for data_size in [:small, :medium, :large]
            for format in [:deflate, :zlib, :gzip]
                t1 = results[(data_size, data_type, format, :zlib, mode)]
                t2 = results[(data_size, data_type, format, :inflate, mode)]
                @printf("%6s %14s %7s: %4.1f %8.1f %8.1f \n",
                        data_size, data_type, format, t2 / t1, 1e6*t1, 1e6*t2)
            end
        end
    end
end

function print_markdown_row(row)
    println(join(vcat("", row, ""), " | "))
end

function print_markdown_table(results, mode)
    mode_string = Dict(:in_memory => "In Memory", :streaming => "Streaming")[mode]
    data_types = [:incompressible, :huffman, :runlength, :graph]
    # TODO: Find a proper way to retrieve the version number from Pkg.
    deps = Pkg.API.Context().env.manifest.deps
    version = only(filter(x -> x.name == "Inflate",
                          collect(values(deps)))).version
    print_markdown_row(vcat(version, mode_string, fill("", 4)))
    print_markdown_row(fill("-", 6))
    print_markdown_row(vcat("", "", string.(data_types)))
    for data_size in [:small, :medium, :large]
        row = [string(data_size)]
        for format in [:deflate, :zlib, :gzip]
            push!(row, string(format))
            for data_type in data_types
                t1 = results[(data_size, data_type, format, :zlib, mode)]
                t2 = results[(data_size, data_type, format, :inflate, mode)]
                push!(row, @sprintf("%4.1f", t2 / t1))
            end
            print_markdown_row(row)
            row = [""]
        end
    end
end

results = run_tests()
print_results(results, :in_memory)
print_results(results, :streaming)
print_markdown_table(results, :in_memory)
print_markdown_table(results, :streaming)

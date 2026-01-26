using BenchmarkTools

using QOI
using PNGFiles
using TestImages
using FixedPointNumbers
using ColorTypes

struct TrialData
    load_time::Float64
    save_time::Float64
    size::Float64
end

struct BenchmarkResult
    file::String
    dim::Tuple{Int, Int}
    size::Int # bytes
    png_data::TrialData
    qoi_data::TrialData
end

function save_png_qoi(img)
    t_png = tempname()
    PNGFiles.save(t_png, img)
    t_qoi = tempname()
    QOI.qoi_encode(t_qoi, img)
    return t_png, t_qoi
end

function run_benchmarks(png_images)
    results = BenchmarkResult[]

    minimum_trial(trial) = minimum(trial.times)

    for (i, png_name) in enumerate(png_images)
        @info "Benchmarking $png_name, $i / $(length(png_images))"
        png_image = testimage(png_name)

        # QOI only supports 24-bit colors so convert to that
        T = eltype(png_image)
        if T <: TransparentColor
            png_image = RGBA{N0f8}.(png_image)
        else
            png_image = RGB{N0f8}.(png_image)
        end

        # Load
        f_png, f_qoi = save_png_qoi(png_image)
        b_png_load = @benchmark PNGFiles.load($f_png)
        b_qoi_load = @benchmark QOI.qoi_decode($f_qoi)

        # Save
        t_png = tempname()
        b_png_save = @benchmark PNGFiles.save($t_png, $png_image)
        t_qoi = tempname()
        b_qoi_save = @benchmark QOI.qoi_encode($t_qoi, $png_image)

        # Sizes
        orig_size = sizeof(png_image)
        png_size = filesize(t_png)
        qoi_size = filesize(t_qoi)

        # Sanity check implementation
        @assert QOI.qoi_decode(f_qoi) == PNGFiles.load(f_png)
        @assert QOI.qoi_decode(t_qoi) == PNGFiles.load(t_png)

        png_data = TrialData(minimum_trial(b_png_load), minimum_trial(b_png_save), png_size)
        qoi_data = TrialData(minimum_trial(b_qoi_load), minimum_trial(b_qoi_save), qoi_size)

        result = BenchmarkResult(png_name, size(png_image), orig_size, png_data, qoi_data)
        push!(results, result)
    end

    return results
end

function render_table(io::IO, results::Vector{BenchmarkResult})
    for result in results
        println(io, "`", result.file, "` ", result.dim[1], "x", result.dim[2], " ", round(result.size / 1024; digits=2), "kB")
        println(io)

        mega_pixels = result.size / (1024^2)

        png_dec_ms = round(result.png_data.load_time / 1e6; digits=2)
        qoi_dec_ms = round(result.qoi_data.load_time / 1e6; digits=2)
        png_enc_ms = round(result.png_data.save_time / 1e6; digits=2)
        qoi_enc_ms = round(result.qoi_data.save_time / 1e6; digits=2)
        png_dec_mpps = round(mega_pixels / (result.png_data.load_time/1e9); digits=2)
        qoi_dec_mpps = round(mega_pixels / (result.qoi_data.load_time/1e9); digits=2)
        png_enc_mpps = round(mega_pixels / (result.png_data.save_time/1e9); digits=2)
        qoi_enc_mpps = round(mega_pixels / (result.qoi_data.save_time/1e9); digits=2)
        png_size_kb = round(result.png_data.size / 1024; digits=2)
        qoi_size_kb = round(result.qoi_data.size / 1024; digits=2)
        png_rate = round(result.png_data.size / result.size * 100; digits=1)
        qoi_rate = round(result.qoi_data.size / result.size * 100; digits=1)

        print(io, """
        |           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
        |:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
        |  PNGFiles  |$png_dec_ms|$png_enc_ms|$png_dec_mpps|$png_enc_mpps|$png_size_kb |$png_rate%|
        |  QOI       |$qoi_dec_ms|$qoi_enc_ms|$qoi_dec_mpps|$qoi_enc_mpps|$qoi_size_kb |$qoi_rate%|

        ----------------------

        """)
    end
end

png_images = filter(endswith(".png"), TestImages.remotefiles)

results = run_benchmarks(png_images)
render_table(stdout, results)

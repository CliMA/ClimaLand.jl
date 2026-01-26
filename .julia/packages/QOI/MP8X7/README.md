# QOI.jl - Implementation of the QOI (Quite OK Image) format 

[![Build Status](https://github.com/KristofferC/QOI.jl/workflows/CI/badge.svg)](https://github.com/KristofferC/QOI.jl/actions?query=workflows/CI)[![codecov](https://codecov.io/gh/KristofferC/QOI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/KristofferC/QOI.jl)


This Julia package contains a decoder and encoder for the [QOI image format](https://qoiformat.org/). The QOI format is very simple and can be faster than using PNG (see the [benchmarks](#benchmarks)) at the cost of a slightly worse compression ratio.

The code here is based on the reference C implementation given in https://github.com/phoboslab/qoi.

## FileIO API

This is the simplest API and likely the one that most should use. Simply, use the `load`/`save` API from FileIO.jl to load and save QOI images.

```jl
using FileIO
image = load("test.qoi")
save("test2.qoi", image)
```

## Basic API

- `QOI.qoi_decode(f::Union{String, IO})` - Read an image in the QOI format from the file/IO `f` and return a matrix with `RGB` or `RGBA` colorants from [ColorTypes.jl](https://github.com/JuliaGraphics/ColorTypes.jl)
- `QOI.qoi_encode(f::Union{String, IO}, image::AbstractMatrix{<:Colorant})` - Write the `image` to the the file/IO `f`. 

## Advanced API

The QOI format is read in row-major order.
This means that a transpose is required to create the matrices returned in the basic API.
To avoid this, the following more advanced APIs exist:

- `QOI.qoi_decode_raw(v::AbstractVector{UIt8}})` - Takes a vector of the bytes of an image in QOI format and returns the uncompressed vector of bytes from decoding.  

- `QOI.qoi_encode_raw(image::AbstractVecOrMat{UInt8}, header::QOI.QOIHeader)` - Returns the bytes from compressing the bytes in `image`. The required header is defined as
  ```jl
  struct QOIHeader
    width::UInt32
    height::UInt32
    channels::QOI.QOIChannel # @enum
    colorspace::QOI.QOIColorSpace # @enum
  end
  ```

  where `channels` can be either `QOI.QOI_RGB` or `QOI.QOI_RGBA` and `colorspace` can be either `QOI.QOI_SRGB` or `QOI.QOI_LINEAR`

- `QOI.qoi_encode_raw!(data::AbstractVecOrMat{UInt8}, image::AbstractVecOrMat{UInt8}, header::QOIHeader)` - Same as above except writing the compressed data into `data`.  

## Benchmarks

The benchmarks here compares the speed of encoding/decoding images with QOI.jl and PNGFiles.jl (which uses `libpng`)
It is supposed to mimic https://qoiformat.org/benchmark/.
The images used in the benchmark are taken from `TestImages.jl`, specifically, all the PNG images in https://testimages.juliaimages.org/stable/imagelist/.
The benchmarks were run on Linux on a 12th Gen Intel(R) Core(TM) i9-12900K CPU.
They can be repeated by running the `benchmark/runbenchmarks.jl` file with the associated environment.

<details><summary>Expand to see benchmark results</summary>


`autumn_leaves.png` 105x140 57.42kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |0.24|0.41|233.52|137.4|22.18 |38.6%|
|  QOI       |0.04|0.13|1418.35|425.17|28.25 |49.2%|

----------------------

`barbara_color.png` 576x720 1215.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |9.05|14.07|131.12|84.35|812.84 |66.9%|
|  QOI       |3.52|5.03|337.19|235.74|945.75 |77.8%|

----------------------

`chelsea.png` 300x451 396.39kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |2.24|4.4|172.82|88.02|218.84 |55.2%|
|  QOI       |0.64|1.55|605.92|250.3|233.27 |58.8%|

----------------------

`coffee.png` 400x600 703.12kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |5.35|11.6|128.28|59.19|443.18 |63.0%|
|  QOI       |1.98|5.23|346.76|131.24|493.3 |70.2%|

----------------------

`fabio_color_256.png` 256x256 192.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |1.56|2.11|119.84|88.84|118.0 |61.5%|
|  QOI       |0.54|0.86|344.22|216.96|154.86 |80.7%|

----------------------

`fabio_color_512.png` 512x512 768.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |5.42|8.0|138.44|93.76|327.45 |42.6%|
|  QOI       |2.04|3.13|366.9|239.9|463.9 |60.4%|

----------------------

`fabio_gray_256.png` 256x256 192.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |1.08|1.65|174.16|113.7|108.65 |56.6%|
|  QOI       |0.38|0.78|491.28|239.6|87.99 |45.8%|

----------------------

`fabio_gray_512.png` 512x512 768.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |3.68|5.84|203.88|128.53|259.41 |33.8%|
|  QOI       |1.55|2.99|483.08|251.05|298.29 |38.8%|

----------------------

`lena_gray_16bit.png` 256x256 192.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |1.12|1.7|167.17|109.98|100.71 |52.5%|
|  QOI       |0.4|0.78|468.04|241.2|82.01 |42.7%|

----------------------

`lighthouse.png` 512x768 1152.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |8.4|13.02|133.94|86.39|692.56 |60.1%|
|  QOI       |3.41|8.14|329.92|138.25|638.75 |55.4%|

----------------------

`monarch_color.png` 512x768 1152.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |8.36|13.34|134.6|84.3|613.23 |53.2%|
|  QOI       |3.64|5.04|309.24|223.43|715.79 |62.1%|

----------------------

`monarch_color_256.png` 256x256 192.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |1.24|2.3|150.65|81.55|129.52 |67.5%|
|  QOI       |0.34|1.46|557.65|128.6|144.39 |75.2%|

----------------------

`mountainstream.png` 512x768 1152.0kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |8.8|13.46|127.78|83.6|892.34 |77.5%|
|  QOI       |3.05|4.41|369.44|255.16|824.83 |71.6%|

----------------------

`toucan.png` 150x162 94.92kB

|           | decode ms | encode ms | decode mpps | encode mpps | size kb | rate |
|:-----------|----------:|----------:|------------:|------------:|--------:|-----:|
|  PNGFiles  |0.34|0.58|274.49|158.87|25.84 |27.2%|
|  QOI       |0.06|0.18|1513.03|514.79|43.28 |45.6%|

</details>

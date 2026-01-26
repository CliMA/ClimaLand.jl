# Release Notes

All notable changes to this package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## Unreleased

## [v1.0.0](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v1.0.0) - 2025-08-29

### The API is now stable
No breaking changes.

## [v0.6.0](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.6.0) - 2025-08-26

### BREAKING `can_concatenate` is now a decoder method instead of a `Codec` method [#73](https://github.com/JuliaIO/ChunkCodecs.jl/pull/73)

### BREAKING the return type of `try_encode`, `try_decode`, and `try_resize_decode!` changed to a new `MaybeSize` type [#72](https://github.com/JuliaIO/ChunkCodecs.jl/pull/72)

## [v0.5.3](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.5.3) - 2025-08-09

- Added support for Julia 1.6 [#68](https://github.com/JuliaIO/ChunkCodecs.jl/pull/68)

## [v0.5.2](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.5.2) - 2025-07-28

- Allow values in `decoded_size_range` to saturate `encode_bound` to `typemax(Int64)`. [#64](https://github.com/JuliaIO/ChunkCodecs.jl/pull/64)
- Added support for Julia 1.9 [#63](https://github.com/JuliaIO/ChunkCodecs.jl/pull/63)

## [v0.5.1](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.5.1) - 2025-07-06

- Added the `is_lossless` trait. [#59](https://github.com/JuliaIO/ChunkCodecs.jl/pull/59)

## [v0.5.0](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.5.0) - 2025-05-23

### BREAKING the resizing behavior in `try_resize_decode!` changed [#45](https://github.com/JuliaIO/ChunkCodecs.jl/pull/45)

`try_resize_decode!` is a `Base.readbytes!` style function where the passed in `dst` vector might be resized if needed. Before this change there was a final `resize!` to shrink `dst` to just contain output if it was resized. The main user facing function `decode` already does a final shrinking `resize!`, so requiring it in `try_resize_decode!` was redundant.

`dst` may now be longer than the returned number of bytes, even if `dst` was grown with `resize!`.

`try_resize_decode!` will now only grow `dst`, never shrink it.

If `max_size` is less than `length(dst)`, instead of throwing an error, `try_resize_decode!` will now act as if `max_size == length(dst)`.

- Added the `grow_dst!` helper function. [#45](https://github.com/JuliaIO/ChunkCodecs.jl/pull/45)

## [v0.4.2](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.4.2) - 2025-04-07

- Added the `decode!` function. [#41](https://github.com/JuliaIO/ChunkCodecs.jl/pull/41)

## [v0.4.1](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.4.1) - 2025-04-06

- Added support for Julia 1.10 [#33](https://github.com/JuliaIO/ChunkCodecs.jl/pull/33)

## [v0.4.0](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.4.0) - 2025-02-24

### BREAKING the `codec` function is replaced with a `.codec` property [#11](https://github.com/JuliaIO/ChunkCodecs.jl/pull/11)

## [v0.3.0](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.3.0) - 2025-01-03

### BREAKING `encode_bound` is required to be monotonically increasing [#7](https://github.com/JuliaIO/ChunkCodecs.jl/pull/7)

### Added

- `ShuffleCodec` and HDF5 compatibility test [#6](https://github.com/JuliaIO/ChunkCodecs.jl/pull/6)

## [v0.2.0](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.2.0) - 2024-12-29

### BREAKING `try_resize_decode!`'s signature is changed. [#5](https://github.com/JuliaIO/ChunkCodecs.jl/pull/5)

`max_size` is a required positional argument instead of an optional keyword argument.

### Fixed

- When using a `Codec` as a decoder the `max_size` option was ignored. [#5](https://github.com/JuliaIO/ChunkCodecs.jl/pull/5)

## [v0.1.1](https://github.com/JuliaIO/ChunkCodecs.jl/tree/ChunkCodecCore-v0.1.1) - 2024-12-20

### Added

- Initial release

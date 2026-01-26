# Wrapping headers

This directory contains a script that can be used to automatically generate wrappers from C headers provided by the SuiteSparse libraries.
This is done using Clang.jl.

# Setup

You need to create an `include` folder that contains the header files `amd.h`, `camd.h`, `ccolamd.h`, `colamd.h` and `SuiteSparse_config.h`.

# Usage

Either run `julia wrapper.jl` directly, or include it and call the `main()` function.
Be sure to activate the project environment in this folder (`julia --project`), which will install `Clang.jl` and `JuliaFormatter.jl`.
You can also call `main(library)` if you want to generate the wrapper for a specific SuiteSparse `library`.
The possible values for `library` are:
- `"all"` (default);
- `"amd"`;
- `"camd"`;
- `"colamd"`;
- `"ccolamd"`.

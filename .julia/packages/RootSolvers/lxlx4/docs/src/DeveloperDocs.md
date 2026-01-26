# Developer Documentation

This section contains documentation for internal methods and developer-focused functionality that may be useful for advanced development and extending RootSolvers.jl.

## Internal Methods

These functions are used internally by the solvers but are exported and may be useful for advanced development.

```@docs
method_args
value_deriv
default_tol
```

## Core Types

### RootSolvingMethod

```@docs
RootSolvingMethod
```

### SolutionType

```@docs
SolutionType
```

## Extending RootSolvers.jl

### Adding New Root-Finding Methods

To add a new root-finding method, you need to:

1. **Define the method struct**:
```julia
struct MyNewMethod{FT} <: RootSolvingMethod{FT}
    x0::FT
    # Add other fields as needed
end
```

2. **Implement `method_args`**:
```julia
method_args(method::MyNewMethod) = (method.x0,)
```

3. **Implement the main `find_zero` method**:
```julia
function find_zero(
    f::F,
    ::MyNewMethod,
    x0::FT,
    soltype::SolutionType,
    tol::AbstractTolerance,
    maxiters::Int,
) where {F <: Function, FT}
    # Your implementation here
    return _find_zero_my_method(f, x0, soltype, tol, maxiters)
end
```

4. **Implement the core algorithm**:
```julia
function _find_zero_my_method(f, x0, soltype, tol, maxiters)
    # Your root-finding algorithm implementation
    # Return a SolutionResults object
end
```

### Adding New Tolerance Types

To add a new tolerance type:

1. **Define the tolerance struct**:
```julia
struct MyTolerance{FT} <: AbstractTolerance{FT}
    tol::FT
end
```

2. **Implement the callable interface**:
```julia
(tol::MyTolerance)(x1, x2, y) = # your convergence criterion
```

### Adding New Solution Types

To add a new solution type:

1. **Define the solution type**:
```julia
struct MySolution <: SolutionType end
```

2. **Define the results struct**:
```julia
struct MySolutionResults{FT} <: AbstractSolutionResults{FT}
    root::FT
    converged::Bool
    # Add other fields as needed
end
```

3. **Implement the constructor**:
```julia
SolutionResults(soltype::MySolution, args...) = MySolutionResults(args...)
```

4. **Implement history functions**:
```julia
init_history(::MySolution, x::FT) where {FT <: Real} = # your initialization
push_history!(history, x, ::MySolution) = # your push logic
```

## Performance Considerations

### GPU Compatibility

When extending RootSolvers.jl for GPU compatibility:

- Use `ifelse` instead of `if-else` blocks where possible
- Avoid dynamic dispatch in hot loops
- Ensure all functions are type-stable
- Use `CompactSolution` for GPU operations (avoid `VerboseSolution`)

### Memory Management

- `CompactSolution` is memory-efficient and GPU-compatible
- `VerboseSolution` stores iteration history and is CPU-only
- Consider memory usage when implementing new solution types

### Type Stability

- Ensure all functions return consistent types
- Use `Base.Fix1` for function composition instead of anonymous functions
- Avoid type instability in hot loops

## Testing Guidelines

### Test Structure

The tests are located in the `test/` directory:

- **`test/runtests.jl`**: Main test suite with comprehensive tests for all methods (including GPU tests)
- **`test/test_helper.jl`**: Helper functions and test utilities
- **`test/test_printing.jl`**: Tests for solution printing and formatting

### What Tests Cover

#### `test/runtests.jl`
- **All root-finding methods**: Secant, Regula Falsi, Brent's, Newton's (AD and manual)
- **All tolerance types**: Solution, Residual, Relative, and combined tolerances
- **All solution types**: Compact and Verbose solutions
- **Edge cases**: Non-finite inputs, convergence failures, high-multiplicity roots
- **Broadcasting**: Array and GPU compatibility
- **Type stability**: Different floating-point types (Float32, Float64)

#### GPU Tests (integrated in `test/runtests.jl`)
- **GPU kernel tests**: CUDA array compatibility
- **Broadcasting on GPU**: Parallel root-finding on GPU arrays

#### `test/test_helper.jl`
- **Test utilities**: Helper functions for generating test problems
- **Method type definitions**: Test-specific method types
- **Problem generators**: Functions to create test cases

#### `test/test_printing.jl`
- **Solution formatting**: Compact and verbose solution display
- **Color output**: Terminal color coding for convergence status
- **History display**: Iteration history formatting

### Running Tests

#### Basic Test Suite
```bash
# From the project root
julia --project=. -e "using Pkg; Pkg.test()"

# Or using the test script directly
julia --project=. test/runtests.jl
```

#### Specific Test Files
```bash
# Run the main test suite (CPU only)
julia --project=. test/runtests.jl

# Run tests with GPU arrays (requires CUDA.jl and compatible GPU)
julia --project=. test/runtests.jl CuArray

# Run printing tests
julia --project=. test/test_printing.jl
```

#### GPU Tests
```bash
# Run all tests including GPU tests (only run when CUDA.jl and compatible GPU is available)
julia --project=. test/runtests.jl CuArray
```

### Test Coverage

The test suite aims for comprehensive coverage:

1. **Unit tests**: Test individual functions and methods
2. **GPU tests**: Test CUDA compatibility and performance (integrated in main test suite)
3. **Edge case tests**: Handle difficult convergence scenarios
4. **Type stability tests**: Ensure GPU compatibility
5. **Broadcasting tests**: Test array and GPU operations

### Continuous Integration

Tests are automatically run on:
- **GitHub Actions**: Multiple Julia versions and operating systems
- **Code coverage**: Tracked via CodeCov
- **GPU tests**: Run on CUDA-compatible runners
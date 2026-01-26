# StructUtils

[![CI](https://github.com/JuliaServices/StructUtils.jl/workflows/CI/badge.svg)](https://github.com/JuliaServices/StructUtils.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaServices/StructUtils.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaServices/StructUtils.jl)
[![deps](https://juliahub.com/docs/StructUtils/deps.svg)](https://juliahub.com/ui/Packages/StructUtils/HHBkp?t=2)
[![version](https://juliahub.com/docs/StructUtils/version.svg)](https://juliahub.com/ui/Packages/StructUtils/HHBkp)
[![pkgeval](https://juliahub.com/docs/StructUtils/pkgeval.svg)](https://juliahub.com/ui/Packages/StructUtils/HHBkp)

## Installation

The package is registered in the [`General`](https://github.com/JuliaRegistries/General) registry and so can be installed at the REPL with `] add StructUtils`.

## Documentation

- [**STABLE**](https://juliaservices.github.io/StructUtils.jl/stable) &mdash; **most recently tagged version of the documentation.**
- [**LATEST**](https://juliaservices.github.io/StructUtils.jl/dev) &mdash; *in-development version of the documentation.*

## Basic Usage

StructUtils.jl provides macros to enhance struct definitions and tools to convert between different data representations.

### Struct Definition Macros

```julia
using StructUtils, Dates

# @noarg - No-argument constructor for mutable structs
@noarg mutable struct User
    id::Int
    name::String
    created_at::DateTime = now()
    active::Bool = true
end

# Create with empty constructor and set values later
user = User()
user.id = 42
user.name = "Alice"

# @defaults - Default values for trailing fields
@defaults struct Point
    x::Float64
    y::Float64
    z::Float64 = 0.0  # Optional with default value
end

# Create with only required fields (trailing defaults are automatically added)
point = Point(1.0, 2.0)  # Point(1.0, 2.0, 0.0)

# @kwarg - Keyword constructor
@kwarg struct Config
    port::Int = 8080
    host::String = "localhost"
    debug::Bool = false
    timeout::Int = 30
end

# Create with keyword arguments (any/all can be omitted if they have defaults)
config = Config(port=9000, debug=true)  # Config(9000, "localhost", true, 30)

# @tags - Field metadata for custom handling
@tags struct Person
    id::Int &(name="person_id",)
    first_name::String &(name="firstName",)
    birth_date::Date &(dateformat="yyyy-mm-dd",)
    internal_note::String = "" &(ignore=true,)
end

# Default values and field tags can be provided in any of the struct macros: @noarg, @kwarg, @defaults, or @tags
@defaults struct Document
    id::Int
    title::String
    created::DateTime &(dateformat="yyyy-mm-dd HH:MM:SS",)
    status::String = "draft" &(json=(name="documentStatus",),)
end
```

### Programmatic Object Construction

The `StructUtils.make` function allows converting between different data representations:

```julia
# Convert a Dict to a struct
data = Dict(:port => 9000, :host => "example.com", :debug => true)
config = StructUtils.make(Config, data)  # Config(9000, "example.com", true, 30)

# Convert a struct to a Dict
dict = StructUtils.make(Dict{Symbol,Any}, config)  
# Dict(:port => 9000, :host => "example.com", :debug => true, :timeout => 30)

# JSON-like nested structures
person_data = Dict(
    "person_id" => 123,
    "firstName" => "Jane",
    "birth_date" => "1990-01-15"
)
person = StructUtils.make(Person, person_data)  # Person(123, "Jane", Date("1990-01-15"), "")

# Convert between array and Vector
StructUtils.make(Vector{Int}, [1, 2, 3])  # [1, 2, 3]

# Create a nested struct with custom field tags
author = Person(1, "John", Date(1980, 5, 10), "VIP customer")
doc_data = Dict(
    "id" => 42,
    "title" => "Annual Report",
    "created" => "2023-04-15 09:30:00",
    "documentStatus" => "published",
    "author" => StructUtils.make(Dict{String,Any}, author)
)

# Type with field references to other types
@tags struct Article
    id::Int
    title::String
    author::Person &(lower=x -> StructUtils.make(Dict{String,Any}, x),)
end

article = StructUtils.make(Article, doc_data)
# Article(42, "Annual Report", Person(1, "John", Date(1980, 5, 10), ""))
```

## Contributing and Questions

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[ci-img]: https://github.com/JuliaServices/StructUtils.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/JuliaServices/StructUtils.jl/actions?query=workflow%3ACI+branch%3Amaster
[codecov-img]: https://codecov.io/gh/JuliaServices/StructUtils.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaServices/StructUtils.jl
[issues-url]: https://github.com/JuliaServices/StructUtils.jl/issues

# StructUtils.jl Documentation

StructUtils.jl provides flexible tools for working with Julia structs, making it easier to build, manipulate, and convert between different data structures. It offers macros for defining struct behaviors and a powerful `make` function for programmatic construction of objects from various data sources.

```@contents
Pages = ["index.md"]
Depth = 3
```

## Installation

The package is registered in the [`General`](https://github.com/JuliaRegistries/General) registry and can be installed at the REPL with:

```julia
] add StructUtils
```

## Quick Start

StructUtils.jl offers several key features:

1. **Struct definition macros** - Enhance struct definitions with special behaviors:
   ```julia
   using StructUtils
   
   # Define a struct with default values
   @defaults struct Config
       port::Int = 8080
       host::String = "localhost"
       debug::Bool = false
   end
   
   # Only need to provide non-default values
   config = Config(9000)  # Config(9000, "localhost", false)
   ```

2. **Programmatic object construction** - Convert between different data representations:
   ```julia
   # Convert a Dict to our Config struct
   dict = Dict(:port => 9000, :host => "example.com")
   config = StructUtils.make(Config, dict)  # Config(9000, "example.com", false)
   
   # Convert a Config back to a Dict
   dict_again = StructUtils.make(Dict{Symbol,Any}, config)  # Dict(:port => 9000, :host => "example.com", :debug => false)
   ```

## Core Concepts

StructUtils.jl is built around several key concepts:

1. **Struct Styles** - Define customization points for handling different struct types
2. **Field Tags** - Add metadata to struct fields for controlling serialization and deserialization
3. **Struct Macros** - Enhance struct definitions with special behaviors
4. **The `make` function** - Programmatically construct objects from various data sources

Let's explore each of these concepts in detail.

## Struct Styles

At the core of StructUtils.jl is the concept of a `StructStyle`:

```julia
abstract type StructStyle end
struct DefaultStyle <: StructStyle end
```

Struct styles provide a way to customize how structs are handled. The `DefaultStyle` is used by default, but you can create custom styles to override behavior for specific types, especially those you don't own:

```julia
struct MyCustomStyle <: StructStyle end

# Override behavior for a type you don't own
StructUtils.lift(::MyCustomStyle, ::Type{UUID}, x::AbstractString) = UUID(x)
```

This approach allows library authors to provide custom serialization/deserialization behavior for their types without modifying the original package.

## Struct Definition Macros

StructUtils.jl provides several macros to enhance struct definitions:

### `@nonstruct` - Non-struct-like Types

The `@nonstruct` macro marks a struct as not struct-like for StructUtils purposes. This is useful for custom types that should be treated as primitives rather than structs during serialization/deserialization:

```julia
@nonstruct struct Email
    value::String
end

# Define conversion methods
StructUtils.lift(::Type{Email}, x::String) = Email(x)
StructUtils.lower(x::Email) = x.value

# Now Email will be serialized as a string, not an object
user = User(email=Email("alice@example.com"))
dict = StructUtils.make(Dict{String,Any}, user)  # email field is a string
```

**Note**: `@nonstruct` does not support field defaults, tags, or other StructUtils macros since you're explicitly opting out of struct-like behavior.

### `@noarg` - No-argument Constructor

The `@noarg` macro creates a no-argument constructor for mutable structs and allows setting default values:

```julia
@noarg mutable struct User
    id::Int
    name::String
    created_at::DateTime = now()
    active::Bool = true
end

# Now you can create a User without arguments
user = User()  # Fields will be undefined except those with defaults
user.id = 1
user.name = "Alice"
```

### `@defaults` - Default Values

The `@defaults` macro creates an additional constructor that allows omitting arguments with default values:

```julia
@defaults struct Point
    x::Float64
    y::Float64
    z::Float64 = 0.0  # Default value
end

# You can omit the z argument
point = Point(1.0, 2.0)  # Point(1.0, 2.0, 0.0)
```

### `@kwarg` - Keyword Constructor

Similar to Base Julia's `@kwdef`, but with enhanced capabilities:

```julia
@kwarg struct HttpConfig
    port::Int = 8080
    host::String = "localhost"
    timeout::Int = 30
end

# Create with keyword arguments
config = HttpConfig(port=9000)  # HttpConfig(9000, "localhost", 30)
```

### `@tags` - Field Metadata

The `@tags` macro allows attaching metadata to struct fields using the `&(...)` syntax:

```julia
@tags struct Person
    id::Int &(json=(name="person_id",),)
    first_name::String &(json=(name="firstName",),)
    birth_date::Date &(dateformat="yyyy-mm-dd",)
    internal_note::String &(json=(ignore=true,),)
end
```

Each field can have tags that control how it's handled by different libraries. For example:

- `name` - Use a different name when serializing/deserializing
- `dateformat` - Specify a format for date parsing/formatting
- `ignore` - Skip this field during serialization/deserialization
- `lift`/`lower` - Custom functions to convert values during serialization/deserialization

## Field Tags Syntax

Field tags use a special syntax: `&(namespace=(key=value,),)`. The namespace (like `json`) allows different libraries to use their own tags without conflicts.

Common field tags include:

- `name`: Alternative name to use when matching source keys
- `ignore`: Skip this field during serialization/deserialization (boolean)
- `dateformat`: Format string or `DateFormat` object for date fields
- `lift`: Function to convert source values to field type
- `lower`: Function to convert field values to serialization format
- `choosetype`: Function to determine concrete type for abstract fields

Example of a field with multiple tags:

```julia
@tags struct Document
    id::Int &(json=(name="doc_id",), db=(column="document_id",))
    created::DateTime &(json=(dateformat="yyyy-mm-dd",), db=(column="creation_date",))
    data::Any &(json=(ignore=true,), db=(ignore=true,))
end
```

## The `make` Function

The core functionality of StructUtils.jl is in the `make` function, which creates an object of a specific type from a source object:

```julia
StructUtils.make(T, source) -> T
StructUtils.make(T, source, style) -> T
```

This function can convert between many kinds of objects:

```julia
# Convert a Dict to a struct
user = StructUtils.make(User, Dict("id" => 1, "name" => "Alice"))

# Convert a struct to a Dict
dict = StructUtils.make(Dict{String,Any}, user)

# Convert a struct to a NamedTuple
nt = StructUtils.make(NamedTuple, user)

# Convert a JSON object to a struct (with JSON.jl)
user = JSON.parse(json_string, User)  # Uses StructUtils.make under the hood
```

### How `make` Works

The `make` function follows these steps:

1. **Type Analysis**: Determine if the target type is:
   - Dictionary-like (`AbstractDict`, `Vector{Pair}`)
   - Array-like (`AbstractArray`, `Tuple`, `Set`)
   - No-arg constructible (`@noarg` or overridden `noarg` function)
   - Regular struct (default constructor)
   - Primitive type (requiring a `lift` function)

2. **Object Construction**:
   - For dictionary-like types: Create an empty dictionary and add key-value pairs
   - For array-like types: Create an empty array and push values
   - For no-arg types: Create an empty instance and set fields
   - For regular structs: Collect field values and call the constructor
   - For primitive types: Use `lift` to convert the source value

3. **Field Mapping**:
   - Match source keys to target fields, respecting field tags
   - Convert values to appropriate field types
   - Handle missing values, defaults, and special types

## Implementing StructUtils Interfaces

To make your types work well with StructUtils.jl, you can implement several interfaces:

### Type Classification

These functions determine how your type is handled by `make`:

```julia
# For dictionary-like types
StructUtils.dictlike(::Type{MyDict}) = true

# For array-like types
StructUtils.arraylike(::Type{MyArray}) = true

# For types with empty constructors
StructUtils.noarg(::Type{MyType}) = true

# For types with keyword constructors
StructUtils.kwdef(::Type{MyType}) = true
```

### Value Conversion

These functions control how values are converted during serialization/deserialization:

```julia
# Convert a source value to your type
StructUtils.lift(::Type{MyType}, x) = MyType(x)

# Convert a key to your type (for dictionary keys)
StructUtils.liftkey(::Type{MyType}, x::String) = MyType(parse(Int, x))

# Convert your type to a serializable form
StructUtils.lower(x::MyType) = string(x)

# Convert your type to a serializable key
StructUtils.lowerkey(x::MyType) = string(x)
```

### Field Metadata

These functions control field behavior:

```julia
# Define default values for fields
StructUtils.fielddefaults(::StructUtils.StructStyle, ::Type{MyType}) = (field1=1, field2="default")

# Define tags for fields
StructUtils.fieldtags(::StructUtils.StructStyle, ::Type{MyType}) = (field1=(name="f1",), field2=(ignore=true,))

# Define a namespace for field tags
StructUtils.fieldtagkey(::MyStyle) = :mylib
```

## Advanced Features

### Type Selection for Abstract Types

When working with abstract types, you need a way to determine the concrete type to construct. The `@choosetype` macro helps with this:

```julia
abstract type Vehicle end
struct Car <: Vehicle; make::String; model::String; end
struct Truck <: Vehicle; make::String; model::String; payload::Float64; end

# Define how to choose concrete types based on source data
StructUtils.@choosetype Vehicle x -> x["type"] == "car" ? Car : Truck

# Now make can create the right type
car = StructUtils.make(Vehicle, Dict("type" => "car", "make" => "Toyota", "model" => "Corolla"))
```

### The Selectors Module

StructUtils includes a `Selectors` module that provides a powerful way to query objects:

```julia
using StructUtils.Selectors

# Create a nested structure
data = Dict("users" => List([
    Dict("id" => 1, "name" => "Alice"),
    Dict("id" => 2, "name" => "Bob")
]))

# Query with selectors
users = data["users"]  # Get the users array
names = users[:].name  # Get all user names
```

The selector syntax supports various operations:
- `x["key"]` / `x.key` - Select by key
- `x[:]` - Select all values
- `x[~, "key"]` - Recursively select all values with key
- `x[:, (k,v) -> Bool]` - Filter by predicate

### Non-Struct-Like Structs with `@nonstruct`

What if you have a custom struct that you want behave more like a primitive type rather than a struct?

The `@nonstruct` macro is perfect for this use case. By marking your struct as non-struct-like, you tell StructUtils to treat it as a primitive type that should be converted directly using `lift` and `lower` methods rather than constructing it from field values.

Here's an example of a custom email type that should be serialized as a JSON string:

```julia
@nonstruct struct Email
    value::String
    
    function Email(value::String)
        # Validate email format
        if !occursin(r"^[^@]+@[^@]+\.[^@]+$", value)
            throw(ArgumentError("Invalid email format: $value"))
        end
        new(value)
    end
end

# Define how to convert from various sources to Email
StructUtils.lift(::Type{Email}, x::String) = Email(x)

# Define how to convert Email to a serializable format
StructUtils.lower(x::Email) = x.value

# Now you can use Email in your structs and it will be serialized as a string
@defaults struct User
    id::Int = 1
    name::String = "default"
    email::Email
end

# Create a user with an email
user = User(email=Email("alice@example.com"))

# Convert to Dict - email will be a string, not an object
dict = StructUtils.make(Dict{String,Any}, user)
# Result: Dict("id" => 1, "name" => "default", "email" => "alice@example.com")

# Convert back from Dict
user_again = StructUtils.make(User, dict)
```

Another example - a custom numeric type that represents a percentage:

```julia
@nonstruct struct Percent <: Number
    value::Float64
    
    function Percent(value::Real)
        if value < 0 || value > 100
            throw(ArgumentError("Percentage must be between 0 and 100"))
        end
        new(Float64(value))
    end
end

# Convert from various numeric types
StructUtils.lift(::Type{Percent}, x::Number) = Percent(x)
StructUtils.lift(::Type{Percent}, x::String) = Percent(parse(Float64, x))

# Convert to a simple number for serialization
StructUtils.lower(x::Percent) = x.value

# Use in a struct
@defaults struct Product
    name::String = "default"
    discount::Percent = Percent(0.0)
end

# Create and serialize
product = Product(discount=Percent(15.5))
dict = StructUtils.make(Dict{String,Any}, product)
# Result: Dict("name" => "default", "discount" => 15.5)
```

The key points about `@nonstruct`:

1. **No field defaults or tags**: Since you're opting out of struct-like behavior, field defaults and tags are not supported.

2. **Requires `lift` and `lower` methods**: You must define how to convert to/from your type.

3. **Fields are private**: The struct's fields are considered implementation details for the `make` process.

4. **Perfect for wrapper types**: Great for types that wrap primitives but need custom validation or behavior.

## Complex Example

Let's put everything together in a complex example, similar to the FrankenStruct example in the JSON.jl documentation:

```julia
using Dates, StructUtils

# Abstract type for polymorphism
abstract type AbstractMonster end

struct Dracula <: AbstractMonster
    num_victims::Int
end

struct Werewolf <: AbstractMonster
    witching_hour::DateTime
end

# Type chooser for AbstractMonster
StructUtils.@choosetype AbstractMonster x -> 
    x isa Dict && haskey(x, "monster_type") && x["monster_type"] == "vampire" ? 
    Dracula : Werewolf

# Custom numeric type with special parsing
struct Percent <: Number
    value::Float64
end

# Custom value lifting
StructUtils.lift(::Type{Percent}, x::Number) = Percent(Float64(x))
StructUtils.liftkey(::Type{Percent}, x::String) = Percent(parse(Float64, x))

# Our complex struct with various field types and defaults
@defaults struct FrankenStruct
    id::Int = 0
    name::String = "Jim"
    address::Union{Nothing, String} = nothing
    rate::Union{Missing, Float64} = missing
    type::Symbol = :a &(json=(name="franken_type",),)
    notsure::Any = nothing
    monster::AbstractMonster = Dracula(0)
    percent::Percent = Percent(0.0)
    birthdate::Date = Date(0) &(dateformat="yyyy/mm/dd",)
    percentages::Dict{Percent, Int} = Dict{Percent, Int}()
    matrix::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
end

# Create a FrankenStruct from a nested dictionary
source = Dict(
    "id" => 1,
    "address" => "123 Main St",
    "franken_type" => "b",
    "monster" => Dict("monster_type" => "vampire", "num_victims" => 10),
    "percent" => 0.1,
    "birthdate" => "2023/10/01",
    "percentages" => Dict("0.1" => 1, "0.2" => 2),
    "matrix" => [[1.0, 2.0], [3.0, 4.0]]
)

# Create a FrankenStruct from the source
frankenstein = StructUtils.make(FrankenStruct, source)

# Convert it back to a dictionary
dict_again = StructUtils.make(Dict{String,Any}, frankenstein)
```

In this example:

1. We define a polymorphic type hierarchy with `AbstractMonster`
2. We implement custom type selection using `@choosetype`
3. We define a custom numeric type `Percent` with special parsing
4. We create a complex struct with various field types and tags
5. We use `make` to create an instance from a nested dictionary

## Summary

StructUtils.jl provides a comprehensive suite of tools for working with Julia structs:

1. **Struct definition macros** enhance structs with special behaviors:
   - `@noarg` - No-argument constructor for mutable structs
   - `@defaults` - Default values for struct fields
   - `@kwarg` - Keyword constructor
   - `@tags` - Field metadata

2. **Field tags** provide metadata for fields:
   - `name` - Alternative name
   - `ignore` - Skip during serialization/deserialization
   - `dateformat` - Format for date fields
   - `lift`/`lower` - Custom conversion functions

3. **The `make` function** converts between different data representations:
   - Dict → Struct
   - Struct → Dict
   - Struct → NamedTuple
   - Array → Vector
   - etc.

4. **Custom interfaces** allow for specialized behavior:
   - Type classification (`dictlike`, `arraylike`, etc.)
   - Value conversion (`lift`, `lower`, etc.)
   - Field metadata (`fielddefaults`, `fieldtags`, etc.)

StructUtils.jl integrates well with other packages like JSON.jl for seamless serialization and deserialization of complex Julia types.


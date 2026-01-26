module CompositionsBase

# Use README as the docstring of the module:
@doc let path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    replace(read(path, String), r"^```julia"m => "```jldoctest README")
end CompositionsBase

if !isdefined(Base, Symbol("@var_str"))
    macro var_str(x)
        return Symbol(x)
    end
end

const compose = ∘

# Since `⨟` may not be rendered correctly in all environments, let's
# use ASCII version as the definition and then `⨟` as an alias.  This
# is not symmetric with how `compose` is defined but appropriately
# "opposite."
opcompose(fs...) = ∘(reverse(fs)...)
const var"⨟" = opcompose

"""
    g ⨟ f
    opcompose(g, f)

The opposite composition operator defined as

    g ⨟ f = f ∘ g
    ⨟(f) = f
    ⨟(fs...) = ∘(reverse(fs)...)
"""
(var"⨟", opcompose)

export ∘, compose, opcompose
@eval export $(Symbol("⨟"))  # for Julia 1.0

if isdefined(Base, :ComposedFunction)
    export decompose, deopcompose
    """
        decompose(comp)

    Destructure a composition into its pieces:
    ```jldoctest
    julia> using CompositionsBase

    julia> decompose(sin ∘ cos)
    (sin, cos)

    julia> decompose(sin∘cos∘identity)
    (sin, cos, identity)

    julia> decompose(sin)
    (sin,)
    ```

    Requires at least Julia version 1.6. See also [deopcompose](@ref).
    """
    decompose(f) = (f,)
    decompose(o::Base.ComposedFunction) = (decompose(o.outer)..., decompose(o.inner)...)

    """
        deopcompose(comp)

    Destructure a composition into its pieces:
    ```jldoctest
    julia> using CompositionsBase

    julia> deopcompose(sin ⨟cos⨟tan)
    (sin, cos, tan)

    julia> decompose(sin⨟cos⨟tan)
    (tan, cos, sin)
    ```

    Requires at least Julia version 1.6. See also [decompose](@ref).
    """
    deopcompose(f) = (f,)
    deopcompose(o::Base.ComposedFunction) = (deopcompose(o.inner)..., deopcompose(o.outer)...)
end

end # module

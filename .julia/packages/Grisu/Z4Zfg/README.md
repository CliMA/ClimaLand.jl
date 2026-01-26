# Grisu

[![Build Status](https://travis-ci.com/JuliaAttic/Grisu.jl.svg?branch=master)](https://travis-ci.com/JuliaAttic/Grisu.jl)

The (internal) Grisu module was removed in Julia 1.6. However, some packages
relies on this module. To keep this working, the Grisu module was filtered out
as a normal package that can be depended on.

Use it as follows, add a dependency on Grisu and use this instead of normally
loading it:

```julia
if isdefined(Base, :Grisu)
    import Base.Grisu
else
    import Grisu
end
```


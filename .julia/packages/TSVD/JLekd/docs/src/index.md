# The TSVD documentation

## Functions
```@meta
CurrentModule = TSVD
DocTestSetup = quote
    using MatrixDepot, TSVD
    try
        matrixdepot("LPnetlib/lp_osa_30")
    catch
    	nothing
    end
end
```

```@docs
tsvd
```

```@meta
DocTestSetup = nothing
```

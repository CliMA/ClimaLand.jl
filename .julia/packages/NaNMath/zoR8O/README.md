# NaNMath

[![CI](https://github.com/JuliaMath/NaNMath.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaMath/NaNMath.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaMath/NaNMath.jl/graph/badge.svg?token=uoFvfoAd4T)](https://codecov.io/gh/JuliaMath/NaNMath.jl)
[![deps](https://juliahub.com/docs/General/NaNMath/stable/deps.svg)](https://juliahub.com/ui/Packages/General/NaNMath?t=2)

Implementations of basic math functions which return ``NaN`` instead of throwing a ``DomainError``.

Example:
```julia
import NaNMath
NaNMath.log(-100) # NaN
NaNMath.pow(-1.5,2.3) # NaN
```

In addition this package provides functions that aggregate arrays and ignore elements that are NaN.
The following functions are implemented:

```
sum
maximum
minimum
extrema
mean
median
var
std
min
max
```

Example:
```julia
using NaNMath; nm=NaNMath
nm.sum([1., 2., NaN]) # result: 3.0
```


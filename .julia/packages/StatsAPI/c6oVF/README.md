# StatsAPI.jl

[![CI](https://github.com/JuliaStats/StatsAPI.jl/workflows/CI/badge.svg)](https://github.com/JuliaStats/StatsAPI.jl/actions?query=workflow%3ACI)
[![deps](https://juliahub.com/docs/StatsAPI/deps.svg)](https://juliahub.com/ui/Packages/StatsAPI/3a8mN?t=2)
[![version](https://juliahub.com/docs/StatsAPI/version.svg)](https://juliahub.com/ui/Packages/StatsAPI/3a8mN)
[![pkgeval](https://juliahub.com/docs/StatsAPI/pkgeval.svg)](https://juliahub.com/ui/Packages/StatsAPI/3a8mN)

### Purpose
This package provides a namespace for statistics-related generic function definitions to solve the optional dependency problem; packages wishing to share and/or extend functions can avoid depending directly on each other by moving the function definition to StatsAPI.jl and each package taking a dependency on it. As such, it is paramount for StatsAPI.jl to be as minimal as possible, defining only generic function stubs and very little else. PRs proposing external dependencies or involved definitions will not be accepted.

### Adding New Functions
When a function is proposed to be defined in StatsAPI.jl, it must include clear documentation of its purpose, convention, and API, as well as specify which package will "own" any generic fallback definitions. Functions will not be exported from StatsAPI.jl, but are left to extending packages to choose whether it is exported from their package or not.

### Extending StatsAPI.jl Functions
Packages wishing to extend a function defined in StatsAPI.jl should first take a dependency on StatsAPI.jl, then define their own method, taking care to properly follow the correct API, and including their package's type in the signature. They *should not* attempt to "pirate" methods with types not owned by the package, or provide their own generic fallbacks, this includes defining methods for `Base` types. If additional methods are desired, discussion should be had with the maintainers of the "owning" package of the generic function. It is recommended to specify version compatibility like `StatsAPI = "1"` in your package, so additional patch release updates to StatsAPI.jl are supported. In the very unlikely event that StatsAPI.jl releases a breaking release 2.0, your package will need to update any broken code and tag a new release with version compatibility like `StatsAPI = "2"`.

### StatsAPI.jl Function Users
*Users* of functions defined in StatsAPI.jl should be aware that taking a dependency on StatsAPI.jl itself is not encouraged, since it exists to coordinate inter-package function sharing, and not to provide any functionality. To actually *use* functions defined in StatsAPI.jl, please review which package "owns" the function, and take a dependency on it, to ensure the default definitions are available, with additional methods being available as non-owning packages are loaded which extend the function defined in StatsAPI.jl.

### StatsAPI.jl Functions
As this package is developer-focused, please see the [source code](https://github.com/JuliaStats/StatsAPI.jl/blob/master/src/StatsAPI.jl) directly for additional information on current functions defined and their documentation.

# ShaderAbstractions

[![Build Status](https://travis-ci.com/SimonDanisch/ShaderAbstractions.jl.svg?branch=master)](https://travis-ci.com/SimonDanisch/ShaderAbstractions.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/SimonDanisch/ShaderAbstractions.jl?svg=true)](https://ci.appveyor.com/project/SimonDanisch/ShaderAbstractions-jl)
[![Codecov](https://codecov.io/gh/SimonDanisch/ShaderAbstractions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SimonDanisch/ShaderAbstractions.jl)

Abstractions to construct shader objects for any WebGL/OpenGL/EGS/Vulkan context!
You construct the objects, and then render them in any backend.

```Julia
using ShaderAbstractions, LinearAlgebra
using ShaderAbstractions: VertexArray
using Test, GeometryBasics
import GeometryBasics

struct WebGL <: ShaderAbstractions.AbstractContext end

m = GLNormalMesh(Sphere(Point3f(0), 1f0))

mvao = VertexArray(m)
instances = VertexArray(positions = rand(GeometryBasics.Point{3, Float32}, 100))

x = ShaderAbstractions.InstancedProgram(
    WebGL(),
    "void main(){}\n", "void main(){}\n",
    mvao,
    instances,
    Dict(
        :model => Mat4f(I),
        :view => Mat4f(I),
        :projection => Mat4f(I),
    )
)
```

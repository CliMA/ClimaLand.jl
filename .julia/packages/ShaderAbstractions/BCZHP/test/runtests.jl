using ShaderAbstractions, LinearAlgebra
using ShaderAbstractions: VertexArray
using Test, GeometryBasics
import GeometryBasics

struct WebGL <: ShaderAbstractions.AbstractContext end

m = normal_mesh(Sphere(Point3f(0), 1f0))

mvao = VertexArray(m)
instances = VertexArray(positions = rand(Point{3, Float32}, 100))

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

@test x.program.fragment_source == replace(read(joinpath(@__DIR__, "test.frag"), String), "\r\n" => "\n")
@test x.program.vertex_source   == replace(read(joinpath(@__DIR__, "test.vert"), String), "\r\n" => "\n")

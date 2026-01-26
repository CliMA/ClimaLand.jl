abstract type ShaderStage end
struct Vertex <: ShaderStage end
struct Geometry <: ShaderStage end
struct Fragment <: ShaderStage end

struct Program
    context::AbstractContext
    vertexarray::VertexArray
    uniforms::Dict{Symbol, Any}
    vertex_source::String
    fragment_source::String
end

struct InstancedProgram
    program::Program
    per_instance::VertexArray
end

function getter_function(io, T, t_str, name)
    println(io, t_str, " get_$(name)(){return $name;}")
end

function getter_function(io, ::AbstractSampler, t_str, name)
end


function input_element(context::AbstractContext, stage::Vertex, io::IO, element::AbstractVector{T}, name::Symbol, uniforms) where {T}
    t_str = type_string(context, T)
    println(io, "in ", t_str, " $name;")
    getkey = Symbol(string(name, "_", "getter"))
    if haskey(uniforms, getkey)
        println(io, uniforms[getkey])
    else
        getter_function(io, T, t_str, name)
    end
end

function input_block(context::AbstractContext, io, vertex_attributes, uniforms)
    for (name, element) in buffers(vertex_attributes)
        input_element(context, Vertex(), io, element, name, uniforms)
    end
end

function InstancedProgram(
        context::AbstractContext,
        vertshader, fragshader,
        instance::Union{GeometryBasics.AbstractMesh, VertexArray},
        per_instance::VertexArray,
        uniforms::Dict{Symbol};
    )
    instance_attributes = sprint() do io
        println(io, "\n// Per instance attributes: ")
        input_block(context, io, per_instance, uniforms)
        println(io)
    end

    p = Program(
        context,
        instance_attributes * vertshader,
        fragshader,
        instance, uniforms
    )
    return InstancedProgram(p, per_instance)
end

function vertex_header(context::AbstractContext)
    return """
    #version 300 es
    precision mediump int;
    precision mediump float;
    precision mediump sampler2D;
    precision mediump sampler3D;
    precision mediump isampler2D;
    precision mediump isampler3D;
    precision mediump usampler2D;
    precision mediump usampler3D;
    """
end

function fragment_header(context::AbstractContext)
    return """
    #version 300 es
    precision mediump int;
    precision mediump float;
    precision mediump sampler2D;
    precision mediump sampler3D;
    precision mediump isampler2D;
    precision mediump isampler3D;
    precision mediump usampler2D;
    precision mediump usampler3D;

    out vec4 fragment_color;
    """
end

function Program(
        context::AbstractContext,
        vertshader, fragshader,
        mesh::Union{VertexArray, GeometryBasics.AbstractMesh}, 
        uniforms::Dict{Symbol}
    )
    converted_uniforms = Dict{Symbol, Any}()
    uniform_block = sprint() do io
        println(io, "\n// Uniforms: ")
        for (name, v) in uniforms
            endswith(string(name), "_getter") && continue
            vc = convert_uniform(context, v)
            t_str = try
               type_string(context, vc)
            catch e
                @error("Type $(typeof(vc)) isn't supported for uniform: $(name)")
                rethrow(e)
            end
            println(io, "uniform ", t_str, " $name;")
            getkey = Symbol(string(name, "_", "getter"))
            if haskey(uniforms, getkey)
                println(io, uniforms[getkey])
            else
                getter_function(io, vc, t_str, name)
            end
            converted_uniforms[name] = vc
        end
        println(io)
    end
    va = VertexArray(mesh)
    src = sprint() do io
        println(io, "// Instance inputs: ")
        input_block(context, io, va, uniforms)
        println(io, uniform_block)
        println(io)
        println(io, vertshader)
    end
    Program(
        context, 
        va,
        converted_uniforms,
        vertex_header(context) * src,
        fragment_header(context) * uniform_block * fragshader
    )
end

import ClimaCore

"""
    NodeProbe(node, from; capacity = 3000)

Records, at every step from the date `from`, the values at the surface node of linear
index `node` (in the parent arrays) of every real-valued field of the state `Y` and of
the cache `p`, keeping the last `capacity` steps. `probe_step!` returns `true` once a
value of `Y` at the node is not finite; `write_probe` writes the records to a CSV file.
"""
mutable struct NodeProbe
    node::Int
    from::DateTime
    capacity::Int
    names::Vector{String}
    rows::Vector{Vector{Float64}}
    dates::Vector{DateTime}
    n_state::Int
end

NodeProbe(node, from; capacity = 3000) =
    NodeProbe(node, from, capacity, String[], Vector{Float64}[], DateTime[], 0)

function probe_collect!(names, values, x, prefix, node)
    if x isa ClimaCore.Fields.Field
        T = eltype(x)
        components = T <: Real ? nothing : fieldnames(T)
        isnothing(components) || all(S -> S <: Real, fieldtypes(T)) || return nothing
        a = parent(x)
        column = Array(view(reshape(a, :, size(a, ndims(a))), :, node))
        nf = isnothing(components) ? 1 : length(components)
        nv = length(column) ÷ nf
        for (k, v) in enumerate(column)
            level, f = mod1(k, nv), cld(k, nv)
            name = isnothing(components) ? prefix : "$prefix.$(components[f])"
            push!(names, nv == 1 ? name : "$name[$level]")
            push!(values, Float64(v))
        end
    elseif x isa ClimaCore.Fields.FieldVector || x isa NamedTuple
        for name in propertynames(x)
            probe_collect!(names, values, getproperty(x, name), "$prefix.$name", node)
        end
    end
    return nothing
end

function probe_step!(probe::NodeProbe, simulation)
    t = date(simulation._integrator.t)
    t < probe.from && return false
    names = String[]
    values = Float64[]
    probe_collect!(names, values, simulation._integrator.u, "Y", probe.node)
    n_state = length(values)
    probe_collect!(names, values, simulation._integrator.p, "p", probe.node)
    if isempty(probe.names)
        probe.names = names
        probe.n_state = n_state
    end
    push!(probe.rows, values)
    push!(probe.dates, t)
    if length(probe.rows) > probe.capacity
        popfirst!(probe.rows)
        popfirst!(probe.dates)
    end
    return !all(isfinite, view(values, 1:n_state))
end

function write_probe(probe::NodeProbe, path)
    isempty(probe.rows) && return nothing
    open(path, "w") do io
        println(io, join(["date"; probe.names], ","))
        for (t, row) in zip(probe.dates, probe.rows)
            length(row) == length(probe.names) || continue
            println(io, join([string(t); string.(row)], ","))
        end
    end
    @info "Probe written" path rows = length(probe.rows)
end

import OrderedCollections: OrderedDict

"""The `DictWriter` is a writer that does not write to disk, but to memory (in a
dictionary).

This is particularly useful for testing and debugging. This is not type stable
(the underlying dictionary does not know in advance what types might be used).
"""
struct DictWriter{T <: AbstractDict} <: AbstractWriter
    """Underlying dictionary. Keys are the short names of the diagnostics, values are
    dictionaries that map the time to the value."""
    dict::T
end

"""
    DictWriter()

A simple in-memory writer. Useful for interactive work and debugging.

You can retrieve values using the typical dictionary interface and using as keys
the names of the stored diagnostics.

Example
========

Assuming we have a diagnostic with short output name "mydiag" stored in `dictW`.
`dictW["mydiag"]` will be a dictionary with keys the timesteps when the data was
saved. The values are the diagnostic output (typically a `ClimaCore` `Field`).
"""
function DictWriter()
    return DictWriter(Dict())
end

"""
    write_field!(writer::DictWriter, field, diagnostic, u, p, t)

Add an entry to the `writer` at time `t` for the current `diagnostic` with value
`field`.

`DictWriter` is backed by a dictionary. Most typically, the keys of this
dictionary are either strings, the `output_short_name` of the diagnostic. If the
`output_short_name` is not available, use the diagnostic itself. The values of
this dictionary is another dictionary that maps the time `t` to the `field` at
that value.

`DictWriter` implements a basic read-only dictionary interface to access the
times and values.
"""
function write_field!(writer::DictWriter, field, diagnostic, u, p, t)
    key_name =
        diagnostic isa ScheduledDiagnostic ? output_short_name(diagnostic) :
        diagnostic
    diagnostic_dict = get!(writer.dict, key_name, OrderedDict())
    diagnostic_dict[t] = copy(field)
    return nothing
end

function Base.getindex(writer::DictWriter, key)
    return Base.getindex(writer.dict, key)
end

function Base.keys(writer::DictWriter)
    return keys(writer.dict)
end

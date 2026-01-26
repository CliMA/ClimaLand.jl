# common errors

@pub struct NotImplementedError <: Exception end
Base.showerror(io::IO, e::NotImplementedError) = print(io, "Not implemented yet")

"""
Raise a `NotImplementedError` exception. This is useful to mark
an interface method as not implemented yet.
"""
@pub not_implemented_error() = throw(NotImplementedError())

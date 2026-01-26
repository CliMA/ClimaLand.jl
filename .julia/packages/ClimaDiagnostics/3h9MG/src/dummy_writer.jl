"""The `DummyWriter` does nothing. It can be used to temporarily disable output of a
scheduled diagnostic. Mostly useful for testing and debugging ClimaDiagnostics.
"""
struct DummyWriter <: AbstractWriter end

# All the methods are the default methods for AbstractWriter

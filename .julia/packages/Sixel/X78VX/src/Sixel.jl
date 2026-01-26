module Sixel

export sixel_encode, sixel_decode

using ImageCore
using OffsetArrays
using IndirectArrays # sixel sequence is actually an indexed image format

import REPL: Terminals

include("interface.jl")
using .SixelInterface
const SEC = AbstractSixelEncoder
const SDC = AbstractSixelDecoder

include("terminaltools.jl")
using .TerminalTools

include("backend/libsixel/LibSixel.jl")
using .LibSixel


# Eventually we will rewrite everything in pure Julia :)
default_encoder(::AbstractArray) = LibSixel.LibSixelEncoder()
default_decoder() = LibSixel.LibSixelDecoder()


# The high-level API to deal with different Julia input types.
#
# We assume that encoder backends perfectly understands `io::IO` and `img::Array`. Any Julia-specific
# fancy array types (e.g., transpose, view) or lazy generators, if supported, should be handled here
# instead of in the backends.
include("encoder.jl")
include("decoder.jl")

# various frontend and input type supports
include("frontend/fileio.jl")

# Ref: https://invisible-island.net/xterm/ctlseqs/ctlseqs.html
"""
    is_sixel_supported(tty=stdout)::Bool

Check if given terminal `tty` supports sixel format.

!!! warning
    (Experiment) The return value is not fully tested on all terminals and all platforms.
"""
function is_sixel_supported(tty::Base.TTY=stdout)
    try
        device_attributes = TerminalTools.query_terminal("\033[0c", tty)
        "4" in split(chop(device_attributes), ';')
    catch
        return false
    end
end
is_sixel_supported(io::IO) = false

end # module

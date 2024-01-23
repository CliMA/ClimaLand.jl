# TimeVaryingInputs.jl
#
# This module contains structs and methods to process external data and evaluate it on the
# model. This module only concerns with evaluations in time, not in space.

# There are three possible sources of data:
# 1. Analytic functions that prescribe how a variable has to be set at a given time
# 2. 0D, single-site data, which is assumed small enough to be saved to memory
# 3. 2D, global data, which cannot be saved to memory in its entirety.
#
# The TimeVaryingInputs module introduces a shared interface for the three cases so that
# uses and developers do not have to worry about the details of what type of data will be
# provided. Behind the scenes, we introduce a new type AbstractTimeVaryingInput, that has
# three concrete implementation, corresponding to the three use cases described above.
# Constructors will automatically identify which of the three implementations to use based
# on the input data, and the existence of three concrete structs should be considered an
# implementation detail.
#
# The three TimeVaryingInputs are:
# - AnalyticTimeVaryingInput,
# - InterpolatingTimeVaryingInput0D,
# - InterpolatingTimeVaryingInput3D.
#
# Along side these TimeVaryingInputs, we also define InterpolationMethods that implement
# specific interpolation strategies (e.g., linear interpolation).
#
# In all cases, the TimeVaryingInputs work with simulation time (ie, seconds from the
# beginning of the reference date). It is up to the various TimeVaryingInputs to convert this
# information to an actual date (if needed).

module TimeVaryingInputs

import ..searchsortednearest

import Adapt
import CUDA
import ClimaComms
import ClimaCore: DeviceSideDevice, DeviceSideContext

export AbstractTimeVaryingInput,
    AbstractInterpolationMethod, TimeVaryingInput, evaluate!

"""
    AbstractTimeVaryingInput

Note
=====

`TimeVaryingInput`s should be considered implementation details. The exposed public interface
should only be considered
- `TimeVaryingInput(input; method, context)` for construction,
- `evaluate!(dest, input, time)` for evaluation
"""
abstract type AbstractTimeVaryingInput end

"""
    AbstractInterpolationMethod

Defines how to perform interpolation.

Not all the TimeVaryingInputs support all the interpolation methods (e.g., no interpolation
methods are supported when the given function is analytic).
"""
abstract type AbstractInterpolationMethod end

"""
    TimeVaryingInput(func)
    TimeVaryingInput(times, vals; method, context)

Construct on object that knows how to evaluate the given function/data on the model times.

When passing a function
=======================

When a function `func` is passed, the function has to be GPU-compatible (e.g., no splines).

When passing single-site data
=============================

When a `times` and `vals` are passed, `times` have to be sorted and the two arrays have to
have the same length.

"""
function TimeVaryingInput end

"""
    evaluate!(dest, input, time)

Evaluate the `input` at the given `time`, writing the output in-place to `dest`.

Depending on the details of `input`, this function might do I/O and communication.
"""
function evaluate! end

include("analytic_time_varying_input.jl")
include("interpolating_time_varying_input0d.jl")

end

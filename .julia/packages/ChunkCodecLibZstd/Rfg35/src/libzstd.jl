# Constants and c wrapper functions ported to Julia from zstd.h and zstd_errors.h
# https://github.com/facebook/zstd/blob/v1.5.6/lib/zstd.h
# https://github.com/facebook/zstd/blob/v1.5.6/lib/zstd_errors.h

"""
    ZSTD_minCLevel()::Int32

Return the minimum negative compression level allowed.
"""
function ZSTD_minCLevel()::Int32
    ccall((:ZSTD_minCLevel, libzstd), Cint, ())
end

"""
    ZSTD_maxCLevel()::Int32

Return the maximum compression level available.
"""
function ZSTD_maxCLevel()::Int32
    ccall((:ZSTD_maxCLevel, libzstd), Cint, ())
end

"""
    ZSTD_defaultCLevel()::Int32

Return the default compression level.
"""
function ZSTD_defaultCLevel()::Int32
    ccall((:ZSTD_defaultCLevel, libzstd), Cint, ())
end

"""
    ZSTD_versionNumber()::VersionNumber

Return the zstd runtime library version.
"""
function ZSTD_versionNumber()::VersionNumber
    v = ccall((:ZSTD_versionNumber, libzstd), Cuint, ())
    major = v ÷ 100_00
    minor = (v ÷ 100) % 100
    patch = v % 100
    VersionNumber(major, minor, patch)
end

"""
    ZSTD_isError(ret::Csize_t)::Bool

Return true if return `ret` is an error.
"""
function ZSTD_isError(ret::Csize_t)::Bool
    !iszero(ccall((:ZSTD_isError, libzstd), Cuint, (Csize_t,), ret))
end

"""
    struct ZSTD_bounds

Has the following fields:
- `error::Csize_t`
- `lowerBound::Cint`
- `upperBound::Cint`
"""
struct ZSTD_bounds
    error::Csize_t
    lowerBound::Cint
    upperBound::Cint
end

"""
    ZSTD_cParam_getBounds(cParam::Int32)::ZSTD_bounds

All parameters must belong to an interval with lower and upper bounds,
otherwise they will either trigger an error or be automatically clamped.
Return : a structure, [`ZSTD_bounds`](@ref), which contains
- an error status field, which must be tested using [`ZSTD_isError`](@ref)
- lower and upper bounds, both inclusive
"""
function ZSTD_cParam_getBounds(cParam::Int32)
    ccall((:ZSTD_cParam_getBounds, libzstd), ZSTD_bounds, (Cint,), cParam)
end

"""
    ZSTD_dParam_getBounds(dParam::Int32)::ZSTD_bounds

All parameters must belong to an interval with lower and upper bounds,
otherwise they will either trigger an error or be automatically clamped.
Return : a structure, [`ZSTD_bounds`](@ref), which contains
- an error status field, which must be tested using [`ZSTD_isError`](@ref)
- lower and upper bounds, both inclusive
"""
function ZSTD_dParam_getBounds(dParam::Int32)
    ccall((:ZSTD_dParam_getBounds, libzstd), ZSTD_bounds, (Cint,), dParam)
end


@assert Culonglong == UInt64
const ZSTD_CONTENTSIZE_UNKNOWN = -1%UInt64
const ZSTD_CONTENTSIZE_ERROR   = -2%UInt64

# convert a result into an error code, which can be compared to error enum list
function ZSTD_getErrorCode(ret::Csize_t)::Cint
    ccall((:ZSTD_getErrorCode, libzstd), Cint, (Csize_t,), ret)
end
# provides readable string from a function result
function ZSTD_getErrorName(ret::Csize_t)::String
    unsafe_string(ccall((:ZSTD_getErrorName, libzstd), Ptr{Cchar}, (Csize_t,), ret))
end

@enum ZSTD_ErrorCode::Cint begin
    ZSTD_error_no_error = 0
    ZSTD_error_GENERIC = 1
    ZSTD_error_prefix_unknown = 10
    ZSTD_error_version_unsupported = 12
    ZSTD_error_frameParameter_unsupported = 14
    ZSTD_error_frameParameter_windowTooLarge = 16
    ZSTD_error_corruption_detected = 20
    ZSTD_error_checksum_wrong = 22
    ZSTD_error_literals_headerWrong = 24
    ZSTD_error_dictionary_corrupted = 30
    ZSTD_error_dictionary_wrong = 32
    ZSTD_error_dictionaryCreation_failed = 34
    ZSTD_error_parameter_unsupported = 40
    ZSTD_error_parameter_combination_unsupported = 41
    ZSTD_error_parameter_outOfBound = 42
    ZSTD_error_tableLog_tooLarge = 44
    ZSTD_error_maxSymbolValue_tooLarge = 46
    ZSTD_error_maxSymbolValue_tooSmall = 48
    ZSTD_error_stabilityCondition_notRespected = 50
    ZSTD_error_stage_wrong = 60
    ZSTD_error_init_missing = 62
    ZSTD_error_memory_allocation = 64
    ZSTD_error_workSpace_tooSmall = 66
    ZSTD_error_dstSize_tooSmall = 70
    ZSTD_error_srcSize_wrong = 72
    ZSTD_error_dstBuffer_null = 74
    ZSTD_error_noForwardProgress_destFull = 80
    ZSTD_error_noForwardProgress_inputEmpty = 82
    # following error codes are __NOT STABLE__, they can be removed or changed in future versions
    # ZSTD_error_frameIndex_tooLarge = 100
    # ZSTD_error_seekableIO = 102
    # ZSTD_error_dstBuffer_wrong = 104
    # ZSTD_error_srcBuffer_wrong = 105
    # ZSTD_error_sequenceProducer_failed = 106
    # ZSTD_error_externalSequences_invalid = 107
    # ZSTD_error_maxCode = 120
end



# Just used to mark the type of pointers
mutable struct ZSTD_CCtx end
mutable struct ZSTD_DCtx end

# set parameter and throw error if it fails
# context must be initialized
function _set_parameter(cctx::Ptr{ZSTD_CCtx}, param::Cint, value::Cint)
    ret = ccall(
        (:ZSTD_CCtx_setParameter, libzstd), Csize_t,
        (Ptr{ZSTD_CCtx}, Cint, Cint),
        cctx, param, value,
    )
    if ZSTD_isError(ret)
        error("$(ZSTD_getErrorName(ret)): setting parameter $(param) to $(value)")
    end
    nothing
end
function _set_parameter(dctx::Ptr{ZSTD_DCtx}, param::Cint, value::Cint)
    ret = ccall(
        (:ZSTD_DCtx_setParameter, libzstd), Csize_t,
        (Ptr{ZSTD_CCtx}, Cint, Cint),
        dctx, param, value,
    )
    if ZSTD_isError(ret)
        error("$(ZSTD_getErrorName(ret)): setting parameter $(param) to $(value)")
    end
    nothing
end


# The following is the original license info from zstd.h and LICENSE file

#=
/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */
=#

#= contents of LICENSE
BSD License

For Zstandard software

Copyright (c) Meta Platforms, Inc. and affiliates. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither the name Facebook, nor Meta, nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#
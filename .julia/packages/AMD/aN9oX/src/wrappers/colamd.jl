function colamd_recommended(nnz, n_row, n_col)
  @ccall libcolamd.colamd_recommended(nnz::Cint, n_row::Cint, n_col::Cint)::Csize_t
end

function colamd_l_recommended(nnz, n_row, n_col)
  @ccall libcolamd.colamd_l_recommended(nnz::SS_Int, n_row::SS_Int, n_col::SS_Int)::Csize_t
end

function colamd_set_defaults(knobs)
  @ccall libcolamd.colamd_set_defaults(knobs::Ptr{Cdouble})::Cvoid
end

function colamd_l_set_defaults(knobs)
  @ccall libcolamd.colamd_l_set_defaults(knobs::Ptr{Cdouble})::Cvoid
end

function colamd(n_row, n_col, Alen, A, p, knobs, stats)
  @ccall libcolamd.colamd(
    n_row::Cint,
    n_col::Cint,
    Alen::Cint,
    A::Ptr{Cint},
    p::Ptr{Cint},
    knobs::Ptr{Cdouble},
    stats::Ptr{Cint},
  )::Cint
end

function colamd_l(n_row, n_col, Alen, A, p, knobs, stats)
  @ccall libcolamd.colamd_l(
    n_row::SS_Int,
    n_col::SS_Int,
    Alen::SS_Int,
    A::Ptr{SS_Int},
    p::Ptr{SS_Int},
    knobs::Ptr{Cdouble},
    stats::Ptr{SS_Int},
  )::SS_Int
end

function symamd(n, A, p, perm, knobs, stats, allocate, release)
  @ccall libcolamd.symamd(
    n::Cint,
    A::Ptr{Cint},
    p::Ptr{Cint},
    perm::Ptr{Cint},
    knobs::Ptr{Cdouble},
    stats::Ptr{Cint},
    allocate::Ptr{Cvoid},
    release::Ptr{Cvoid},
  )::Cint
end

function symamd_l(n, A, p, perm, knobs, stats, allocate, release)
  @ccall libcolamd.symamd_l(
    n::SS_Int,
    A::Ptr{SS_Int},
    p::Ptr{SS_Int},
    perm::Ptr{SS_Int},
    knobs::Ptr{Cdouble},
    stats::Ptr{SS_Int},
    allocate::Ptr{Cvoid},
    release::Ptr{Cvoid},
  )::SS_Int
end

function colamd_report(stats)
  @ccall libcolamd.colamd_report(stats::Ptr{Cint})::Cvoid
end

function colamd_l_report(stats)
  @ccall libcolamd.colamd_l_report(stats::Ptr{SS_Int})::Cvoid
end

function symamd_report(stats)
  @ccall libcolamd.symamd_report(stats::Ptr{Cint})::Cvoid
end

function symamd_l_report(stats)
  @ccall libcolamd.symamd_l_report(stats::Ptr{SS_Int})::Cvoid
end

const COLAMD_KNOBS = 20
const COLAMD_STATS = 20

const COLAMD_DENSE_ROW = 1
const COLAMD_DENSE_COL = 2
const COLAMD_AGGRESSIVE = 3
const COLAMD_DEFRAG_COUNT = 3
const COLAMD_STATUS = 4
const COLAMD_INFO1 = 5
const COLAMD_INFO2 = 6
const COLAMD_INFO3 = 7

const COLAMD_OK = 0
const COLAMD_OK_BUT_JUMBLED = 1

const COLAMD_ERROR_A_not_present = -1
const COLAMD_ERROR_p_not_present = -2
const COLAMD_ERROR_nrow_negative = -3
const COLAMD_ERROR_ncol_negative = -4
const COLAMD_ERROR_nnz_negative = -5
const COLAMD_ERROR_p0_nonzero = -6
const COLAMD_ERROR_A_too_small = -7
const COLAMD_ERROR_col_length_negative = -8
const COLAMD_ERROR_row_index_out_of_bounds = -9
const COLAMD_ERROR_out_of_memory = -10
const COLAMD_ERROR_internal_error = -999

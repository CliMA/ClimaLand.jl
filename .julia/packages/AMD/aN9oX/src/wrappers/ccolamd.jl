function ccolamd_recommended(nnz, n_row, n_col)
  @ccall libccolamd.ccolamd_recommended(nnz::Cint, n_row::Cint, n_col::Cint)::Csize_t
end

function ccolamd_l_recommended(nnz, n_row, n_col)
  @ccall libccolamd.ccolamd_l_recommended(nnz::SS_Int, n_row::SS_Int, n_col::SS_Int)::Csize_t
end

function ccolamd_set_defaults(knobs)
  @ccall libccolamd.ccolamd_set_defaults(knobs::Ptr{Cdouble})::Cvoid
end

function ccolamd_l_set_defaults(knobs)
  @ccall libccolamd.ccolamd_l_set_defaults(knobs::Ptr{Cdouble})::Cvoid
end

function ccolamd(n_row, n_col, Alen, A, p, knobs, stats, cmember)
  @ccall libccolamd.ccolamd(
    n_row::Cint,
    n_col::Cint,
    Alen::Cint,
    A::Ptr{Cint},
    p::Ptr{Cint},
    knobs::Ptr{Cdouble},
    stats::Ptr{Cint},
    cmember::Ptr{Cint},
  )::Cint
end

function ccolamd_l(n_row, n_col, Alen, A, p, knobs, stats, cmember)
  @ccall libccolamd.ccolamd_l(
    n_row::SS_Int,
    n_col::SS_Int,
    Alen::SS_Int,
    A::Ptr{SS_Int},
    p::Ptr{SS_Int},
    knobs::Ptr{Cdouble},
    stats::Ptr{SS_Int},
    cmember::Ptr{SS_Int},
  )::SS_Int
end

function csymamd(n, A, p, perm, knobs, stats, allocate, release, cmember, stype)
  @ccall libccolamd.csymamd(
    n::Cint,
    A::Ptr{Cint},
    p::Ptr{Cint},
    perm::Ptr{Cint},
    knobs::Ptr{Cdouble},
    stats::Ptr{Cint},
    allocate::Ptr{Cvoid},
    release::Ptr{Cvoid},
    cmember::Ptr{Cint},
    stype::Cint,
  )::Cint
end

function csymamd_l(n, A, p, perm, knobs, stats, allocate, release, cmember, stype)
  @ccall libccolamd.csymamd_l(
    n::SS_Int,
    A::Ptr{SS_Int},
    p::Ptr{SS_Int},
    perm::Ptr{SS_Int},
    knobs::Ptr{Cdouble},
    stats::Ptr{SS_Int},
    allocate::Ptr{Cvoid},
    release::Ptr{Cvoid},
    cmember::Ptr{SS_Int},
    stype::SS_Int,
  )::SS_Int
end

function ccolamd_report(stats)
  @ccall libccolamd.ccolamd_report(stats::Ptr{Cint})::Cvoid
end

function ccolamd_l_report(stats)
  @ccall libccolamd.ccolamd_l_report(stats::Ptr{SS_Int})::Cvoid
end

function csymamd_report(stats)
  @ccall libccolamd.csymamd_report(stats::Ptr{Cint})::Cvoid
end

function csymamd_l_report(stats)
  @ccall libccolamd.csymamd_l_report(stats::Ptr{SS_Int})::Cvoid
end

function ccolamd2(
  n_row,
  n_col,
  Alen,
  A,
  p,
  knobs,
  stats,
  Front_npivcol,
  Front_nrows,
  Front_ncols,
  Front_parent,
  Front_cols,
  p_nfr,
  InFront,
  cmember,
)
  @ccall libccolamd.ccolamd2(
    n_row::Cint,
    n_col::Cint,
    Alen::Cint,
    A::Ptr{Cint},
    p::Ptr{Cint},
    knobs::Ptr{Cdouble},
    stats::Ptr{Cint},
    Front_npivcol::Ptr{Cint},
    Front_nrows::Ptr{Cint},
    Front_ncols::Ptr{Cint},
    Front_parent::Ptr{Cint},
    Front_cols::Ptr{Cint},
    p_nfr::Ptr{Cint},
    InFront::Ptr{Cint},
    cmember::Ptr{Cint},
  )::Cint
end

function ccolamd2_l(
  n_row,
  n_col,
  Alen,
  A,
  p,
  knobs,
  stats,
  Front_npivcol,
  Front_nrows,
  Front_ncols,
  Front_parent,
  Front_cols,
  p_nfr,
  InFront,
  cmember,
)
  @ccall libccolamd.ccolamd2_l(
    n_row::SS_Int,
    n_col::SS_Int,
    Alen::SS_Int,
    A::Ptr{SS_Int},
    p::Ptr{SS_Int},
    knobs::Ptr{Cdouble},
    stats::Ptr{SS_Int},
    Front_npivcol::Ptr{SS_Int},
    Front_nrows::Ptr{SS_Int},
    Front_ncols::Ptr{SS_Int},
    Front_parent::Ptr{SS_Int},
    Front_cols::Ptr{SS_Int},
    p_nfr::Ptr{SS_Int},
    InFront::Ptr{SS_Int},
    cmember::Ptr{SS_Int},
  )::SS_Int
end

function ccolamd_apply_order(Front, Order, Temp, nn, nfr)
  @ccall libccolamd.ccolamd_apply_order(
    Front::Ptr{Cint},
    Order::Ptr{Cint},
    Temp::Ptr{Cint},
    nn::Cint,
    nfr::Cint,
  )::Cvoid
end

function ccolamd_l_apply_order(Front, Order, Temp, nn, nfr)
  @ccall libccolamd.ccolamd_l_apply_order(
    Front::Ptr{SS_Int},
    Order::Ptr{SS_Int},
    Temp::Ptr{SS_Int},
    nn::SS_Int,
    nfr::SS_Int,
  )::Cvoid
end

function ccolamd_fsize(nn, MaxFsize, Fnrows, Fncols, Parent, Npiv)
  @ccall libccolamd.ccolamd_fsize(
    nn::Cint,
    MaxFsize::Ptr{Cint},
    Fnrows::Ptr{Cint},
    Fncols::Ptr{Cint},
    Parent::Ptr{Cint},
    Npiv::Ptr{Cint},
  )::Cvoid
end

function ccolamd_l_fsize(nn, MaxFsize, Fnrows, Fncols, Parent, Npiv)
  @ccall libccolamd.ccolamd_l_fsize(
    nn::SS_Int,
    MaxFsize::Ptr{SS_Int},
    Fnrows::Ptr{SS_Int},
    Fncols::Ptr{SS_Int},
    Parent::Ptr{SS_Int},
    Npiv::Ptr{SS_Int},
  )::Cvoid
end

function ccolamd_postorder(
  nn,
  Parent,
  Npiv,
  Fsize,
  Order,
  Child,
  Sibling,
  Stack,
  Front_cols,
  cmember,
)
  @ccall libccolamd.ccolamd_postorder(
    nn::Cint,
    Parent::Ptr{Cint},
    Npiv::Ptr{Cint},
    Fsize::Ptr{Cint},
    Order::Ptr{Cint},
    Child::Ptr{Cint},
    Sibling::Ptr{Cint},
    Stack::Ptr{Cint},
    Front_cols::Ptr{Cint},
    cmember::Ptr{Cint},
  )::Cvoid
end

function ccolamd_l_postorder(
  nn,
  Parent,
  Npiv,
  Fsize,
  Order,
  Child,
  Sibling,
  Stack,
  Front_cols,
  cmember,
)
  @ccall libccolamd.ccolamd_l_postorder(
    nn::SS_Int,
    Parent::Ptr{SS_Int},
    Npiv::Ptr{SS_Int},
    Fsize::Ptr{SS_Int},
    Order::Ptr{SS_Int},
    Child::Ptr{SS_Int},
    Sibling::Ptr{SS_Int},
    Stack::Ptr{SS_Int},
    Front_cols::Ptr{SS_Int},
    cmember::Ptr{SS_Int},
  )::Cvoid
end

function ccolamd_post_tree(root, k, Child, Sibling, Order, Stack)
  @ccall libccolamd.ccolamd_post_tree(
    root::Cint,
    k::Cint,
    Child::Ptr{Cint},
    Sibling::Ptr{Cint},
    Order::Ptr{Cint},
    Stack::Ptr{Cint},
  )::Cint
end

function ccolamd_l_post_tree(root, k, Child, Sibling, Order, Stack)
  @ccall libccolamd.ccolamd_l_post_tree(
    root::SS_Int,
    k::SS_Int,
    Child::Ptr{SS_Int},
    Sibling::Ptr{SS_Int},
    Order::Ptr{SS_Int},
    Stack::Ptr{SS_Int},
  )::SS_Int
end

const CCOLAMD_KNOBS = 20
const CCOLAMD_STATS = 20

const CCOLAMD_DENSE_ROW = 1
const CCOLAMD_DENSE_COL = 2
const CCOLAMD_AGGRESSIVE = 3
const CCOLAMD_LU = 4
const CCOLAMD_DEFRAG_COUNT = 3
const CCOLAMD_STATUS = 4
const CCOLAMD_INFO1 = 5
const CCOLAMD_INFO2 = 6
const CCOLAMD_INFO3 = 7
const CCOLAMD_EMPTY_ROW = 8
const CCOLAMD_EMPTY_COL = 9
const CCOLAMD_NEWLY_EMPTY_ROW = 10
const CCOLAMD_NEWLY_EMPTY_COL = 11

const CCOLAMD_OK = 0
const CCOLAMD_OK_BUT_JUMBLED = 1

const CCOLAMD_ERROR_A_not_present = -1
const CCOLAMD_ERROR_p_not_present = -2
const CCOLAMD_ERROR_nrow_negative = -3
const CCOLAMD_ERROR_ncol_negative = -4
const CCOLAMD_ERROR_nnz_negative = -5
const CCOLAMD_ERROR_p0_nonzero = -6
const CCOLAMD_ERROR_A_too_small = -7
const CCOLAMD_ERROR_col_length_negative = -8
const CCOLAMD_ERROR_row_index_out_of_bounds = -9
const CCOLAMD_ERROR_out_of_memory = -10
const CCOLAMD_ERROR_invalid_cmember = -11
const CCOLAMD_ERROR_internal_error = -999

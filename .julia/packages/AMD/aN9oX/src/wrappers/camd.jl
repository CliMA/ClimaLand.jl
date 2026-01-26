function camd_order(n, Ap, Ai, P, Control, Info, C)
  @ccall libcamd.camd_order(
    n::Cint,
    Ap::Ptr{Cint},
    Ai::Ptr{Cint},
    P::Ptr{Cint},
    Control::Ptr{Cdouble},
    Info::Ptr{Cdouble},
    C::Ptr{Cint},
  )::Cint
end

function camd_l_order(n, Ap, Ai, P, Control, Info, C)
  @ccall libcamd.camd_l_order(
    n::SS_Int,
    Ap::Ptr{SS_Int},
    Ai::Ptr{SS_Int},
    P::Ptr{SS_Int},
    Control::Ptr{Cdouble},
    Info::Ptr{Cdouble},
    C::Ptr{SS_Int},
  )::SS_Int
end

function camd_2(
  n,
  Pe,
  Iw,
  Len,
  iwlen,
  pfree,
  Nv,
  Next,
  Last,
  Head,
  Elen,
  Degree,
  W,
  Control,
  Info,
  C,
  BucketSet,
)
  @ccall libcamd.camd_2(
    n::Cint,
    Pe::Ptr{Cint},
    Iw::Ptr{Cint},
    Len::Ptr{Cint},
    iwlen::Cint,
    pfree::Cint,
    Nv::Ptr{Cint},
    Next::Ptr{Cint},
    Last::Ptr{Cint},
    Head::Ptr{Cint},
    Elen::Ptr{Cint},
    Degree::Ptr{Cint},
    W::Ptr{Cint},
    Control::Ptr{Cdouble},
    Info::Ptr{Cdouble},
    C::Ptr{Cint},
    BucketSet::Ptr{Cint},
  )::Cvoid
end

function camd_l2(
  n,
  Pe,
  Iw,
  Len,
  iwlen,
  pfree,
  Nv,
  Next,
  Last,
  Head,
  Elen,
  Degree,
  W,
  Control,
  Info,
  C,
  BucketSet,
)
  @ccall libcamd.camd_l2(
    n::SS_Int,
    Pe::Ptr{SS_Int},
    Iw::Ptr{SS_Int},
    Len::Ptr{SS_Int},
    iwlen::SS_Int,
    pfree::SS_Int,
    Nv::Ptr{SS_Int},
    Next::Ptr{SS_Int},
    Last::Ptr{SS_Int},
    Head::Ptr{SS_Int},
    Elen::Ptr{SS_Int},
    Degree::Ptr{SS_Int},
    W::Ptr{SS_Int},
    Control::Ptr{Cdouble},
    Info::Ptr{Cdouble},
    C::Ptr{SS_Int},
    BucketSet::Ptr{SS_Int},
  )::Cvoid
end

function camd_valid(n_row, n_col, Ap, Ai)
  @ccall libcamd.camd_valid(n_row::Cint, n_col::Cint, Ap::Ptr{Cint}, Ai::Ptr{Cint})::Cint
end

function camd_l_valid(n_row, n_col, Ap, Ai)
  @ccall libcamd.camd_l_valid(
    n_row::SS_Int,
    n_col::SS_Int,
    Ap::Ptr{SS_Int},
    Ai::Ptr{SS_Int},
  )::SS_Int
end

function camd_cvalid(n, C)
  @ccall libcamd.camd_cvalid(n::Cint, C::Ptr{Cint})::Cint
end

function camd_l_cvalid(n, C)
  @ccall libcamd.camd_l_cvalid(n::SS_Int, C::Ptr{SS_Int})::SS_Int
end

function camd_defaults(Control)
  @ccall libcamd.camd_defaults(Control::Ptr{Cdouble})::Cvoid
end

function camd_l_defaults(Control)
  @ccall libcamd.camd_l_defaults(Control::Ptr{Cdouble})::Cvoid
end

function camd_control(Control)
  @ccall libcamd.camd_control(Control::Ptr{Cdouble})::Cvoid
end

function camd_l_control(Control)
  @ccall libcamd.camd_l_control(Control::Ptr{Cdouble})::Cvoid
end

function camd_info(Info)
  @ccall libcamd.camd_info(Info::Ptr{Cdouble})::Cvoid
end

function camd_l_info(Info)
  @ccall libcamd.camd_l_info(Info::Ptr{Cdouble})::Cvoid
end

const CAMD_CONTROL = 5
const CAMD_INFO = 20

const CAMD_DENSE = 1
const CAMD_AGGRESSIVE = 2
const CAMD_DEFAULT_DENSE = 10.0
const CAMD_DEFAULT_AGGRESSIVE = 2

const CAMD_STATUS = 1
const CAMD_N = 2
const CAMD_NZ = 3
const CAMD_SYMMETRY = 4
const CAMD_NZDIAG = 5
const CAMD_NZ_A_PLUS_AT = 6
const CAMD_NDENSE = 7
const CAMD_MEMORY = 8
const CAMD_NCMPA = 9
const CAMD_LNZ = 10
const CAMD_NDIV = 11
const CAMD_NMULTSUBS_LDL = 12
const CAMD_NMULTSUBS_LU = 13
const CAMD_DMAX = 14

const CAMD_OK = 0
const CAMD_OUT_OF_MEMORY = -1
const CAMD_INVALID = -2
const CAMD_OK_BUT_JUMBLED = 1

function amd_order(n, Ap, Ai, P, Control, Info)
  @ccall libamd.amd_order(
    n::Cint,
    Ap::Ptr{Cint},
    Ai::Ptr{Cint},
    P::Ptr{Cint},
    Control::Ptr{Cdouble},
    Info::Ptr{Cdouble},
  )::Cint
end

function amd_l_order(n, Ap, Ai, P, Control, Info)
  @ccall libamd.amd_l_order(
    n::SS_Int,
    Ap::Ptr{SS_Int},
    Ai::Ptr{SS_Int},
    P::Ptr{SS_Int},
    Control::Ptr{Cdouble},
    Info::Ptr{Cdouble},
  )::SS_Int
end

function amd_2(n, Pe, Iw, Len, iwlen, pfree, Nv, Next, Last, Head, Elen, Degree, W, Control, Info)
  @ccall libamd.amd_2(
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
  )::Cvoid
end

function amd_l2(n, Pe, Iw, Len, iwlen, pfree, Nv, Next, Last, Head, Elen, Degree, W, Control, Info)
  @ccall libamd.amd_l2(
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
  )::Cvoid
end

function amd_valid(n_row, n_col, Ap, Ai)
  @ccall libamd.amd_valid(n_row::Cint, n_col::Cint, Ap::Ptr{Cint}, Ai::Ptr{Cint})::Cint
end

function amd_l_valid(n_row, n_col, Ap, Ai)
  @ccall libamd.amd_l_valid(n_row::SS_Int, n_col::SS_Int, Ap::Ptr{SS_Int}, Ai::Ptr{SS_Int})::SS_Int
end

function amd_defaults(Control)
  @ccall libamd.amd_defaults(Control::Ptr{Cdouble})::Cvoid
end

function amd_l_defaults(Control)
  @ccall libamd.amd_l_defaults(Control::Ptr{Cdouble})::Cvoid
end

function amd_control(Control)
  @ccall libamd.amd_control(Control::Ptr{Cdouble})::Cvoid
end

function amd_l_control(Control)
  @ccall libamd.amd_l_control(Control::Ptr{Cdouble})::Cvoid
end

function amd_info(Info)
  @ccall libamd.amd_info(Info::Ptr{Cdouble})::Cvoid
end

function amd_l_info(Info)
  @ccall libamd.amd_l_info(Info::Ptr{Cdouble})::Cvoid
end

const AMD_CONTROL = 5
const AMD_INFO = 20

const AMD_DENSE = 1
const AMD_AGGRESSIVE = 2
const AMD_DEFAULT_DENSE = 10.0
const AMD_DEFAULT_AGGRESSIVE = 2

const AMD_STATUS = 1
const AMD_N = 2
const AMD_NZ = 3
const AMD_SYMMETRY = 4
const AMD_NZDIAG = 5
const AMD_NZ_A_PLUS_AT = 6
const AMD_NDENSE = 7
const AMD_MEMORY = 8
const AMD_NCMPA = 9
const AMD_LNZ = 10
const AMD_NDIV = 11
const AMD_NMULTSUBS_LDL = 12
const AMD_NMULTSUBS_LU = 13
const AMD_DMAX = 14

const AMD_OK = 0
const AMD_OUT_OF_MEMORY = -1
const AMD_INVALID = -2
const AMD_OK_BUT_JUMBLED = 1

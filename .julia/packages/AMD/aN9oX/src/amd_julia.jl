export Amd, amd_valid, amd

"""
    amd(A)
    amd(A, meta)

Given a sparse matrix `A` and an `Amd` structure `meta`, `p = amd(A, meta)` computes the approximate minimum degree ordering of `A + Aᵀ`.
The ordering is represented as a permutation vector `p`.
Factorizations of `A[p,p]` tend to be sparser than those of `A`.

The matrix `A` must be square and the sparsity pattern of `A + Aᵀ` is implicit.
Thus it is convenient to represent symmetric or hermitian matrices using one triangle only.
The diagonal of `A` may be present but will be ignored.

The ordering may be influenced by changing `meta.control[AMD_DENSE]` and `meta.control[AMD_AGGRESSIVE]`.

Statistics on the ordering appear in `meta.info`.
"""
function amd end

const amd_statuses = Dict(
  AMD_OK => "ok",
  AMD_OUT_OF_MEMORY => "out of memory",
  AMD_INVALID => "input invalid",
  AMD_OK_BUT_JUMBLED => "ok but jumbled",
)

"""
Base type to hold control and information related to a call to AMD.
`control` is a `Vector{Float64}` with components:

* d = control[AMD_DENSE]: rows with more than max(d√n, 16) entries are considered "dense" and appear last in the permutation.
If d < 0 no row will be treated as dense.

* control[AMD_AGGRESSIVE]: triggers aggressive absorption if nonzero.

`info` is a `Vector{Float64}` that contains statistics on the ordering.
"""
mutable struct Amd
  control::Vector{Cdouble}
  info::Vector{Cdouble}

  function Amd()
    control = zeros(Cdouble, AMD_CONTROL)
    info = zeros(Cdouble, AMD_INFO)
    amd_defaults(control)
    return new(control, info)
  end
end

function show(io::IO, meta::Amd)
  s = "Control:\n"
  s *= "  dense row parameter: $(meta.control[AMD_DENSE])\n"
  s *= "  aggressive absorption: $(meta.control[AMD_AGGRESSIVE])\n"
  s *= "Info:\n"
  s *= "  status: $(amd_statuses[meta.info[AMD_STATUS]])\n"
  s *= "  matrix size: $(meta.info[AMD_N])\n"
  s *= "  number of nonzeros: $(meta.info[AMD_NZ])\n"
  s *= "  pattern symmetry: $(meta.info[AMD_SYMMETRY])\n"
  s *= "  number of nonzeros on diagonal: $(meta.info[AMD_NZDIAG])\n"
  s *= "  number of nonzeros in A + Aᵀ: $(meta.info[AMD_NZ_A_PLUS_AT])\n"
  s *= "  number of dense columns: $(meta.info[AMD_NDENSE])\n"
  s *= "  memory used: $(meta.info[AMD_MEMORY])\n"
  s *= "  number of garbage collections: $(meta.info[AMD_NCMPA])\n"
  s *= "  approx number of nonzers in factor: $(meta.info[AMD_LNZ])\n"
  s *= "  number of float divides: $(meta.info[AMD_NDIV])\n"
  s *= "  number of float * or - for LDL: $(meta.info[AMD_NMULTSUBS_LDL])\n"
  s *= "  number of float * or - for LU: $(meta.info[AMD_NMULTSUBS_LU])\n"
  s *= "  max nonzeros in any column of factor: $(meta.info[AMD_DMAX])\n"
  print(io, s)
end

print(io::IO, meta::Amd) = show(io, meta)

for (validfn, typ) in ((:amd_valid, :Cint), (:amd_l_valid, :SS_Int))
  Base.Sys.WORD_SIZE == 32 && validfn == :amd_l_valid && continue
  @eval begin
    function amd_valid(A::SparseMatrixCSC{F, $typ}) where {F}
      nrow, ncol = size(A)
      colptr = A.colptr .- $typ(1)  # 0-based indexing
      rowval = A.rowval .- $typ(1)
      valid = $validfn(nrow, ncol, colptr, rowval)
      return valid == AMD_OK || valid == AMD_OK_BUT_JUMBLED
    end

    amd_valid(A::Symmetric{F, SparseMatrixCSC{F, $typ}}) where {F} = amd_valid(A.data)
    amd_valid(A::Hermitian{F, SparseMatrixCSC{F, $typ}}) where {F} = amd_valid(A.data)
  end
end

for (orderfn, typ) in ((:amd_order, :Cint), (:amd_l_order, :SS_Int))
  Base.Sys.WORD_SIZE == 32 && orderfn == :amd_l_order && continue
  @eval begin
    function amd(A::SparseMatrixCSC{F, $typ}, meta::Amd) where {F}
      nrow, ncol = size(A)
      nrow == ncol || error("AMD: input matrix must be square")
      colptr = A.colptr .- $typ(1)  # 0-based indexing
      rowval = A.rowval .- $typ(1)

      p = zeros($typ, nrow)
      valid = $orderfn(nrow, colptr, rowval, p, meta.control, meta.info)
      (valid == AMD_OK || valid == AMD_OK_BUT_JUMBLED) ||
        throw("amd_order returns: $(amd_statuses[valid])")
      p .+= 1
      return p
    end

    amd(A::Symmetric{F, SparseMatrixCSC{F, $typ}}, meta::Amd) where {F} = amd(A.data, meta)
    amd(A::Hermitian{F, SparseMatrixCSC{F, $typ}}, meta::Amd) where {F} = amd(A.data, meta)

    function amd(A::SparseMatrixCSC{F, $typ}) where {F}
      meta = Amd()
      amd(A, meta)
    end

    amd(A::Symmetric{F, SparseMatrixCSC{F, $typ}}) where {F} = amd(A.data)
    amd(A::Hermitian{F, SparseMatrixCSC{F, $typ}}) where {F} = amd(A.data)
  end
end

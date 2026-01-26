export Colamd, colamd, symamd

"""
    colamd(A)
    colamd(A, meta)

colamd computes a permutation vector `p` such that the Cholesky factorization of
  `A[:,p]' * A[:,p]` has less fill-in and requires fewer floating point operations than `AᵀA`.
"""
function colamd end

"""
    symamd(A)
    symamd(A, meta)

Given a symmetric or hermitian matrix `A`, symamd computes a permutation vector `p`
such that the Cholesky factorization of `A[p,p]` has less fill-in and
requires fewer floating point operations than that of A.
"""
function symamd end

const colamd_statuses = Dict(
  COLAMD_OK => "ok",
  COLAMD_OK_BUT_JUMBLED => "ok, but columns of input matrix were jumbled (unsorted columns or duplicate entries)",
  COLAMD_ERROR_A_not_present => "rowval is a null pointer",
  COLAMD_ERROR_p_not_present => "colptr is a null pointer",
  COLAMD_ERROR_nrow_negative => "nrow is negative",
  COLAMD_ERROR_ncol_negative => "nccol is negative",
  COLAMD_ERROR_nnz_negative => "number of nonzeros in matrix is negative",
  COLAMD_ERROR_p0_nonzero => "p[0] is nonzero",
  COLAMD_ERROR_A_too_small => "workspace is too small",
  COLAMD_ERROR_col_length_negative => "a column has a negative number of entries",
  COLAMD_ERROR_row_index_out_of_bounds => "a row index is out of bounds",
  COLAMD_ERROR_out_of_memory => "out of memory",
  COLAMD_ERROR_internal_error => "internal error",
)

mutable struct Colamd{T <: Union{Cint, SS_Int}}
  knobs::Vector{Cdouble}
  stats::Vector{T}

  function Colamd{T}() where {T <: Union{Cint, SS_Int}}
    knobs = zeros(Cdouble, COLAMD_KNOBS)
    stats = zeros(T, COLAMD_STATS)
    colamd_set_defaults(knobs)
    return new(knobs, stats)
  end
end

function show(io::IO, meta::Colamd)
  s = "dense row parameter: $(meta.knobs[COLAMD_DENSE_ROW])\n"
  s *= "dense col parameter: $(meta.knobs[COLAMD_DENSE_COL])\n"
  s *= "aggressive absorption: $(meta.knobs[COLAMD_AGGRESSIVE])\n"
  s *= "memory defragmentation: $(meta.stats[COLAMD_DEFRAG_COUNT])\n"
  s *= "status: $(colamd_statuses[meta.stats[COLAMD_STATUS]])\n"
  print(io, s)
end

print(io::IO, meta::Colamd) = show(io, meta)

for (fn, typ) in ((:colamd, :Cint), (:colamd_l, :SS_Int))
  Base.Sys.WORD_SIZE == 32 && fn == :colamd_l && continue
  @eval begin
    function colamd(A::SparseMatrixCSC{F, $typ}, meta::Colamd{$typ}) where {F}
      nrow, ncol = size(A)
      nnz = A.colptr[end] - 1
      p = A.colptr .- $typ(1)  # 0-based indexing
      len = colamd_recommended($typ(nnz), $typ(nrow), $typ(ncol))
      workspace = zeros($typ, len)
      workspace[1:length(A.rowval)] .= A.rowval .- $typ(1)
      valid = $fn(nrow, ncol, len, workspace, p, meta.knobs, meta.stats)
      Bool(valid) || throw("colamd status: $(colamd_statuses[meta.stats[COLAMD_STATUS]])")
      pop!(p)  # remove the number of nnz
      p .+= $typ(1)  # 1-based indexing
      return p
    end

    colamd(A::Symmetric{F, SparseMatrixCSC{F, $typ}}, meta::Colamd{$typ}) where {F} =
      colamd(A.data, meta)
    colamd(A::Hermitian{F, SparseMatrixCSC{F, $typ}}, meta::Colamd{$typ}) where {F} =
      colamd(A.data, meta)

    function colamd(A::SparseMatrixCSC{F, $typ}) where {F}
      meta = Colamd{$typ}()
      colamd(A, meta)
    end

    colamd(A::Symmetric{F, SparseMatrixCSC{F, $typ}}) where {F} = colamd(A.data)
    colamd(A::Hermitian{F, SparseMatrixCSC{F, $typ}}) where {F} = colamd(A.data)
  end
end

for (fn, typ) in ((:symamd, :Cint), (:symamd_l, :SS_Int))
  Base.Sys.WORD_SIZE == 32 && fn == :symamd_l && continue
  @eval begin
    function symamd(A::SparseMatrixCSC{F, $typ}, meta::Colamd{$typ}) where {F}
      nrow, ncol = size(A)
      colptr = A.colptr .- $typ(1)  # 0-based indexing
      rowval = A.rowval .- $typ(1)
      p = zeros($typ, nrow + 1) # p is used as a workspace during the ordering, which is why it must be of length n+1, not just n
      cfun_calloc = @cfunction(Base.Libc.calloc, Ptr{Cvoid}, ($typ, $typ))
      cfun_free = @cfunction(Base.Libc.free, Cvoid, (Ptr{Cvoid},))
      valid = $fn(nrow, rowval, colptr, p, meta.knobs, meta.stats, cfun_calloc, cfun_free)
      Bool(valid) || throw("symamd status: $(colamd_statuses[meta.stats[COLAMD_STATUS]])")
      pop!(p)
      p .+= $typ(1)  # 1-based indexing
      return p
    end

    symamd(A::Symmetric{F, SparseMatrixCSC{F, $typ}}, meta::Colamd{$typ}) where {F} =
      symamd(A.data, meta)
    symamd(A::Hermitian{F, SparseMatrixCSC{F, $typ}}, meta::Colamd{$typ}) where {F} =
      symamd(A.data, meta)

    function symamd(A::SparseMatrixCSC{F, $typ}) where {F}
      meta = Colamd{$typ}()
      symamd(A, meta)
    end

    symamd(A::Symmetric{F, SparseMatrixCSC{F, $typ}}) where {F} = symamd(A.data)
    symamd(A::Hermitian{F, SparseMatrixCSC{F, $typ}}) where {F} = symamd(A.data)
  end
end

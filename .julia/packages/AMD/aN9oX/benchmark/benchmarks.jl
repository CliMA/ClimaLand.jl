using BenchmarkTools
using DelimitedFiles, LinearAlgebra, Printf, SparseArrays
using Pkg.Artifacts

using MatrixMarket
using Metis, SymRCM, AMD
using LDLFactorizations, SolverBenchmark

const SUITE = BenchmarkGroup()
SUITE["no_ordering"] = BenchmarkGroup()
SUITE["amd"] = BenchmarkGroup()
SUITE["symamd"] = BenchmarkGroup()
SUITE["metis"] = BenchmarkGroup()
SUITE["symrcm"] = BenchmarkGroup()

# obtain path to SQD collection
const artifact_toml = joinpath(@__DIR__, "Artifacts.toml")
ensure_artifact_installed("sqdcollection", artifact_toml)
const sqd_hash = artifact_hash("sqdcollection", artifact_toml)
@assert artifact_exists(sqd_hash)
const sqd_path = joinpath(artifact_path(sqd_hash), "sqd-collection-0.1")

subdirs = readdir(sqd_path)
const formulations = ["2x2", "3x3"]

names_sqd = String[]
ratio_amd = Float64[]
ratio_symamd = Float64[]
ratio_metis = Float64[]
ratio_symrcm = Float64[]
ratio_classic = Float64[]
nnz_sqd = Int[]

for subdir ∈ subdirs
  subdir == ".git" && continue
  isdir(joinpath(sqd_path, subdir)) || continue  # ignore regular files
  for formulation ∈ formulations
    iterpath = joinpath(sqd_path, subdir, formulation, "iter_0")
    isdir(iterpath) || continue

    A = MatrixMarket.mmread(joinpath(iterpath, "K_0.mtx"))
    b = readdlm(joinpath(iterpath, "rhs_0.rhs"), Float64)[:]
    nnz_A = nnz(tril(A, -1))
    n = size(A, 1)
    if n ≤ 25000
      name = "$(subdir)_$(formulation)"
      push!(names_sqd, name)

      p_amd = amd(A)
      p_symamd = symamd(A)
      p_metis, _ = Metis.permutation(A)
      p_symrcm = symrcm(A)
      p_classic = collect(1:n)
      println(name)

      SUITE["no_ordering"][name] = @benchmarkable nnz(ldl($A, $p_classic).L) / $nnz_A
      SUITE["amd"][name] = @benchmarkable nnz(ldl($A, $p_amd).L) / $nnz_A
      SUITE["symamd"][name] = @benchmarkable nnz(ldl($A, $p_symamd).L) / $nnz_A
      SUITE["metis"][name] = @benchmarkable nnz(ldl($A, $p_metis).L) / $nnz_A
      SUITE["symrcm"][name] = @benchmarkable nnz(ldl($A, $p_symrcm).L) / $nnz_A
    end
  end
end

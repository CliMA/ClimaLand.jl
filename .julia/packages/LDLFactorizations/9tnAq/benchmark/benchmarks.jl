# benchmark file to compare two commits of LDLFactorizations

using BenchmarkTools
using MatrixMarket

using DelimitedFiles
using LinearAlgebra
using Pkg.Artifacts
using Printf
using SparseArrays

using LDLFactorizations

# obtain path to SQD collection
const artifact_toml = joinpath(@__DIR__, "Artifacts.toml")
ensure_artifact_installed("sqdcollection", artifact_toml)
const sqd_hash = artifact_hash("sqdcollection", artifact_toml)
@assert artifact_exists(sqd_hash)
const sqd_path = joinpath(artifact_path(sqd_hash), "sqd-collection-0.1")
subdirs = readdir(sqd_path)
const formulations = ("2x2", "3x3")
const iters = (0, 5, 10)

const SUITE = BenchmarkGroup()
SUITE["analyze"] = BenchmarkGroup()
SUITE["factorize"] = BenchmarkGroup()
SUITE["solve1"] = BenchmarkGroup()
SUITE["solve5"] = BenchmarkGroup()

for subdir ∈ subdirs
  subdir == ".git" && continue
  isdir(joinpath(sqd_path, subdir)) || continue  # ignore regular files
  for formulation ∈ formulations
    for iter ∈ iters
      iterpath = joinpath(sqd_path, subdir, formulation, "iter_$(iter)")
      isdir(iterpath) || continue
      A = MatrixMarket.mmread(joinpath(iterpath, "K_$(iter).mtx"))
      b = readdlm(joinpath(iterpath, "rhs_$(iter).rhs"))[:, 1]
      B = [b b b b b]
      name = "$(subdir)_$(formulation)_$(iter)"
      SUITE["analyze"][name] = @benchmarkable ldl_analyze($A)
      LDL = ldl_analyze(A)
      SUITE["factorize"][name] = @benchmarkable ldl_factorize!($A, $LDL)
      ldl_factorize!(A, LDL)
      x = similar(b)
      SUITE["solve1"][name] = @benchmarkable ldiv!($x, $LDL, $b)
      X = similar(B)
      SUITE["solve5"][name] = @benchmarkable ldiv!($X, $LDL, $B)
    end
  end
end

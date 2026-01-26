# benchmark LDLFactorizations vs QDLDL

using BenchmarkTools
using MatrixMarket

using LinearAlgebra
using Pkg.Artifacts
using Printf
using SparseArrays

using LDLFactorizations
using QDLDL

# obtain path to SQD collection
const artifact_toml = joinpath(@__DIR__, "Artifacts.toml")
const sqd_hash = artifact_hash("sqdcollection", artifact_toml)
@assert artifact_exists(sqd_hash)
const sqd_path = joinpath(artifact_path(sqd_hash), "sqd-collection-0.1")
subdirs = readdir(sqd_path)
const formulations = ("2x2", "3x3")
const iters = (0, 5, 10)

const SUITE = BenchmarkGroup()

for subdir ∈ subdirs
  subdir == ".git" && continue
  isdir(joinpath(sqd_path, subdir)) || continue  # ignore regular files
  for formulation ∈ formulations
    for iter ∈ iters
      iterpath = joinpath(sqd_path, subdir, formulation, "iter_$(iter)")
      isdir(iterpath) || continue
      A = MatrixMarket.mmread(joinpath(iterpath, "K_$(iter).mtx"))
      name = "$(subdir)_$(formulation)_$(iter)"
      SUITE[name] = BenchmarkGroup()
      SUITE[name]["LDL"] = @benchmarkable ldl($A)
      SUITE[name]["QDLDL"] = @benchmarkable qdldl($A)
    end
  end
end

using ADTypes
using Test

@test ADTypes.Auto(:ChainRules, 1) isa AutoChainRules{Int64}
@test ADTypes.Auto(:Diffractor) isa AutoDiffractor
@test ADTypes.Auto(:Enzyme) isa AutoEnzyme
@test ADTypes.Auto(:FastDifferentiation) isa AutoFastDifferentiation
@test ADTypes.Auto(:FiniteDiff) isa AutoFiniteDiff
@test ADTypes.Auto(:FiniteDifferences, 1.0) isa AutoFiniteDifferences{Float64}
@test ADTypes.Auto(:ForwardDiff) isa AutoForwardDiff
@test ADTypes.Auto(:Mooncake) isa AutoMooncake
@test ADTypes.Auto(:PolyesterForwardDiff) isa AutoPolyesterForwardDiff
@test ADTypes.Auto(:ReverseDiff) isa AutoReverseDiff
@test ADTypes.Auto(:Symbolics) isa AutoSymbolics
@test ADTypes.Auto(:Tapir) isa AutoTapir
@test ADTypes.Auto(:Tracker) isa AutoTracker
@test ADTypes.Auto(:Zygote) isa AutoZygote
@test ADTypes.Auto(nothing) isa NoAutoDiff

@test_throws MethodError ADTypes.Auto(:ThisPackageDoesNotExist)
@test_throws UndefKeywordError ADTypes.Auto(:ChainRules)
@test_throws UndefKeywordError ADTypes.Auto(:FiniteDifferences)

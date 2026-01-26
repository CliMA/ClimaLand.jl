using ADTypes
using ADTypes: AbstractADType,
               mode,
               ForwardMode,
               ForwardOrReverseMode,
               ReverseMode,
               SymbolicMode
using ADTypes: dense_ad,
               NoSparsityDetector,
               KnownJacobianSparsityDetector,
               KnownHessianSparsityDetector,
               sparsity_detector,
               jacobian_sparsity,
               hessian_sparsity,
               NoColoringAlgorithm,
               coloring_algorithm,
               column_coloring,
               row_coloring,
               symmetric_coloring
using Aqua: Aqua
using ChainRulesCore: ChainRulesCore, RuleConfig,
                      HasForwardsMode, HasReverseMode,
                      NoForwardsMode, NoReverseMode
using EnzymeCore: EnzymeCore
using JET: JET
using Test

## Backend-specific

struct CustomTag end

struct ForwardRuleConfig <: RuleConfig{Union{HasForwardsMode, NoReverseMode}} end
struct ReverseRuleConfig <: RuleConfig{Union{NoForwardsMode, HasReverseMode}} end
struct ForwardOrReverseRuleConfig <: RuleConfig{Union{HasForwardsMode, HasReverseMode}} end

struct FakeSparsityDetector <: ADTypes.AbstractSparsityDetector end
struct FakeColoringAlgorithm <: ADTypes.AbstractColoringAlgorithm end

function every_ad()
    return [
        AutoChainRules(; ruleconfig = :rc),
        AutoDiffractor(),
        AutoEnzyme(),
        AutoFastDifferentiation(),
        AutoFiniteDiff(),
        AutoFiniteDifferences(; fdm = :fdm),
        AutoForwardDiff(),
        AutoGTPSA(),
        AutoPolyesterForwardDiff(),
        AutoReverseDiff(),
        AutoSymbolics(),
        AutoTapir(),
        AutoTracker(),
        AutoZygote()
    ]
end

function every_ad_with_options()
    return [
        AutoChainRules(; ruleconfig = :rc),
        AutoDiffractor(),
        AutoEnzyme(),
        AutoEnzyme(mode = :forward),
        AutoFastDifferentiation(),
        AutoFiniteDiff(),
        AutoFiniteDiff(fdtype = :fd, fdjtype = :fdj, fdhtype = :fdh,
            relstep = 1, absstep = 2, dir = false),
        AutoFiniteDifferences(; fdm = :fdm),
        AutoForwardDiff(),
        AutoForwardDiff(chunksize = 3, tag = :tag),
        AutoGTPSA(),
        AutoGTPSA(descriptor = Val(:descriptor)),
        AutoMooncake(; config = :config),
        AutoMooncakeForward(; config = :config),
        AutoPolyesterForwardDiff(),
        AutoPolyesterForwardDiff(chunksize = 3, tag = :tag),
        AutoReverseDiff(),
        AutoReverseDiff(compile = true),
        AutoSymbolics(),
        AutoTapir(),
        AutoTapir(safe_mode = false),
        AutoTracker(),
        AutoZygote()
    ]
end

## Tests

@testset verbose=true "ADTypes.jl" begin
    if VERSION >= v"1.10"
        @testset "Aqua.jl" begin
            Aqua.test_all(ADTypes; deps_compat = (check_extras = false,))
        end
        @testset "JET.jl" begin
            JET.test_package(ADTypes, target_defined_modules = true)
        end
    end
    @testset "Dense" begin
        include("dense.jl")
    end
    @testset "Sparse" begin
        include("sparse.jl")
    end
    @testset "Symbols" begin
        include("symbols.jl")
    end
    @testset "Legacy" begin
        include("legacy.jl")
    end
    @testset "Miscellaneous" begin
        include("misc.jl")
    end
    if VERSION >= v"1.11.0-DEV.469"
        @testset "Public" begin
            include("public.jl")
        end
    end
end

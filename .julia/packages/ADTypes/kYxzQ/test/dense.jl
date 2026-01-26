@testset "AutoChainRules" begin
    ad = AutoChainRules(; ruleconfig = ForwardOrReverseRuleConfig())
    @test ad isa AbstractADType
    @test ad isa AutoChainRules{ForwardOrReverseRuleConfig}
    @test mode(ad) isa ForwardOrReverseMode
    @test ad.ruleconfig == ForwardOrReverseRuleConfig()

    ad = AutoChainRules(; ruleconfig = ForwardRuleConfig())
    @test ad isa AbstractADType
    @test ad isa AutoChainRules{ForwardRuleConfig}
    @test mode(ad) isa ForwardMode
    @test ad.ruleconfig == ForwardRuleConfig()

    ad = AutoChainRules(; ruleconfig = ReverseRuleConfig())
    @test ad isa AbstractADType
    @test ad isa AutoChainRules{ReverseRuleConfig}
    @test mode(ad) isa ReverseMode
    @test ad.ruleconfig == ReverseRuleConfig()
end

@testset "AutoDiffractor" begin
    ad = AutoDiffractor()
    @test ad isa AbstractADType
    @test ad isa AutoDiffractor
    @test mode(ad) isa ForwardOrReverseMode
end

@testset "AutoEnzyme" begin
    ad = AutoEnzyme()
    @test ad isa AbstractADType
    @test ad isa AutoEnzyme{Nothing, Nothing}
    @test mode(ad) isa ForwardOrReverseMode
    @test ad.mode === nothing

    ad = AutoEnzyme(; mode = EnzymeCore.Forward)
    @test ad isa AbstractADType
    @test ad isa AutoEnzyme{typeof(EnzymeCore.Forward), Nothing}
    @test mode(ad) isa ForwardMode
    @test ad.mode == EnzymeCore.Forward

    ad = AutoEnzyme(; function_annotation = EnzymeCore.Const)
    @test ad isa AbstractADType
    @test ad isa AutoEnzyme{Nothing, EnzymeCore.Const}
    @test mode(ad) isa ForwardOrReverseMode
    @test ad.mode === nothing

    ad = AutoEnzyme(;
        mode = EnzymeCore.Reverse, function_annotation = EnzymeCore.Duplicated)
    @test ad isa AbstractADType
    @test ad isa AutoEnzyme{typeof(EnzymeCore.Reverse), EnzymeCore.Duplicated}
    @test mode(ad) isa ReverseMode
    @test ad.mode == EnzymeCore.Reverse
end

@testset "AutoFastDifferentiation" begin
    ad = AutoFastDifferentiation()
    @test ad isa AbstractADType
    @test ad isa AutoFastDifferentiation
    @test mode(ad) isa SymbolicMode
end

@testset "AutoFiniteDiff" begin
    ad = AutoFiniteDiff()
    @test ad isa AbstractADType
    @test ad isa AutoFiniteDiff
    @test mode(ad) isa ForwardMode
    @test ad.fdtype === Val(:forward)
    @test ad.fdjtype === Val(:forward)
    @test ad.fdhtype === Val(:hcentral)
    @test ad.relstep === nothing
    @test ad.absstep === nothing
    @test ad.dir

    ad = AutoFiniteDiff(; fdtype = Val(:central), fdjtype = Val(:forward),
        relstep = 1e-3, absstep = 1e-4, dir = false)
    @test ad isa AbstractADType
    @test ad isa AutoFiniteDiff
    @test mode(ad) isa ForwardMode
    @test ad.fdtype === Val(:central)
    @test ad.fdjtype === Val(:forward)
    @test ad.fdhtype === Val(:hcentral)
    @test ad.relstep == 1e-3
    @test ad.absstep == 1e-4
    @test !ad.dir
end

@testset "AutoFiniteDifferences" begin
    ad = AutoFiniteDifferences(; fdm = nothing)
    @test ad isa AbstractADType
    @test ad isa AutoFiniteDifferences{Nothing}
    @test mode(ad) isa ForwardMode
    @test ad.fdm === nothing

    ad = AutoFiniteDifferences(; fdm = Val(:forward_fdm))
    @test ad isa AbstractADType
    @test ad isa AutoFiniteDifferences{Val{:forward_fdm}}
    @test mode(ad) isa ForwardMode
    @test ad.fdm == Val(:forward_fdm)
end

@testset "AutoForwardDiff" begin
    ad = AutoForwardDiff()
    @test ad isa AbstractADType
    @test ad isa AutoForwardDiff{nothing, Nothing}
    @test mode(ad) isa ForwardMode
    @test ad.tag === nothing

    ad = AutoForwardDiff(; chunksize = 10, tag = CustomTag())
    @test ad isa AbstractADType
    @test ad isa AutoForwardDiff{10, CustomTag}
    @test mode(ad) isa ForwardMode
    @test ad.tag == CustomTag()
end

@testset "AutoGTPSA" begin
    ad = AutoGTPSA(; descriptor = nothing)
    @test ad isa AbstractADType
    @test ad isa AutoGTPSA{Nothing}
    @test mode(ad) isa ForwardMode
    @test ad.descriptor === nothing

    ad = AutoGTPSA(; descriptor = Val(:descriptor))
    @test ad isa AbstractADType
    @test ad isa AutoGTPSA{Val{:descriptor}}
    @test mode(ad) isa ForwardMode
    @test ad.descriptor == Val(:descriptor)
end

@testset "AutoMooncake" begin
    ad = AutoMooncake(; config = :config)
    @test ad isa AbstractADType
    @test ad isa AutoMooncake
    @test mode(ad) isa ReverseMode
    @test ad.config == :config
    ad = AutoMooncake()
    @test ad.config === nothing
end

@testset "AutoMooncakeForward" begin
    ad = AutoMooncakeForward(; config = :config)
    @test ad isa AbstractADType
    @test ad isa AutoMooncakeForward
    @test mode(ad) isa ForwardMode
    @test ad.config === :config
    ad = AutoMooncakeForward()
    @test ad.config === nothing
end

@testset "AutoPolyesterForwardDiff" begin
    ad = AutoPolyesterForwardDiff()
    @test ad isa AbstractADType
    @test ad isa AutoPolyesterForwardDiff{nothing, Nothing}
    @test mode(ad) isa ForwardMode
    @test ad.tag === nothing

    ad = AutoPolyesterForwardDiff(; chunksize = 10, tag = CustomTag())
    @test ad isa AbstractADType
    @test ad isa AutoPolyesterForwardDiff{10, CustomTag}
    @test mode(ad) isa ForwardMode
    @test ad.tag == CustomTag()
end

@testset "AutoReverseDiff" begin
    ad = @inferred AutoReverseDiff()
    @test ad isa AbstractADType
    @test ad isa AutoReverseDiff{false}
    @test mode(ad) isa ReverseMode
    @test !ad.compile
    @test_deprecated ad.compile

    ad = AutoReverseDiff(; compile = true)
    @test ad isa AbstractADType
    @test ad isa AutoReverseDiff{true}
    @test mode(ad) isa ReverseMode
    @test ad.compile
    @test_deprecated ad.compile

    ad = @inferred AutoReverseDiff(; compile = Val(true))
    @test ad isa AbstractADType
    @test ad isa AutoReverseDiff{true}
    @test mode(ad) isa ReverseMode
    @test ad.compile
    @test_deprecated ad.compile
end

@testset "AutoSymbolics" begin
    ad = AutoSymbolics()
    @test ad isa AbstractADType
    @test ad isa AutoSymbolics
    @test mode(ad) isa SymbolicMode
end

@testset "AutoTapir" begin
    ad = AutoTapir()
    @test ad isa AbstractADType
    @test ad isa AutoTapir
    @test mode(ad) isa ReverseMode
    @test ad.safe_mode

    ad = AutoTapir(; safe_mode = false)
    @test !ad.safe_mode
end

@testset "AutoTaylorDiff" begin
    ad = AutoTaylorDiff(; order = 2)
    @test ad isa AbstractADType
    @test ad isa AutoTaylorDiff{2}
    @test mode(ad) isa ForwardMode

    ad = AutoTaylorDiff()
    @test ad isa AbstractADType
    @test ad isa AutoTaylorDiff{1}
    @test mode(ad) isa ForwardMode
end

@testset "AutoTracker" begin
    ad = AutoTracker()
    @test ad isa AbstractADType
    @test ad isa AutoTracker
    @test mode(ad) isa ReverseMode
end

@testset "AutoZygote" begin
    ad = AutoZygote()
    @test ad isa AbstractADType
    @test ad isa AutoZygote
    @test mode(ad) isa ReverseMode
end

@testset "NoAutoDiff" begin
    ad = NoAutoDiff()
    @test ad isa AbstractADType
    @test ad isa NoAutoDiff
    @test_throws NoAutoDiffSelectedError mode(ad)
end

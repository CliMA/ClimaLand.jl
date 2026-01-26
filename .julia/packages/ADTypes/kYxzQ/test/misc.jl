@testset "Broadcasting" begin
    for ad in every_ad()
        @test identity.(ad) == ad
    end
end

@testset "Printing" begin
    for ad in every_ad_with_options()
        @test startswith(string(ad), "Auto")
        @test contains(string(ad), "(")
        @test endswith(string(ad), ")")
    end

    sparse_backend1 = AutoSparse(AutoForwardDiff())
    sparse_backend2 = AutoSparse(
        AutoForwardDiff();
        sparsity_detector = FakeSparsityDetector(),
        coloring_algorithm = FakeColoringAlgorithm()
    )
    @test contains(string(sparse_backend1), string(AutoForwardDiff()))
    @test length(string(sparse_backend1)) < length(string(sparse_backend2))
end

#=
The following tests are only for visual assessment of the printing behavior.
They do not correspond to proper use of ADTypes constructors.
Please refer to the docstrings for that.
=#
for backend in [
    # dense
    ADTypes.AutoChainRules(; ruleconfig = :rc),
    ADTypes.AutoDiffractor(),
    ADTypes.AutoEnzyme(),
    ADTypes.AutoEnzyme(mode = :forward),
    ADTypes.AutoEnzyme(function_annotation = Val{:forward}),
    ADTypes.AutoEnzyme(mode = :reverse, function_annotation = Val{:duplicated}),
    ADTypes.AutoFastDifferentiation(),
    ADTypes.AutoFiniteDiff(),
    ADTypes.AutoFiniteDiff(fdtype = :fd, fdjtype = :fdj, fdhtype = :fdh),
    ADTypes.AutoFiniteDifferences(; fdm = :fdm),
    ADTypes.AutoForwardDiff(),
    ADTypes.AutoForwardDiff(chunksize = 3, tag = :tag),
    ADTypes.AutoGTPSA(),
    ADTypes.AutoGTPSA(; descriptor = Val(:descriptor)),
    ADTypes.AutoMooncake(),
    ADTypes.AutoMooncake(; config = :config),
    ADTypes.AutoMooncakeForward(),
    ADTypes.AutoMooncakeForward(; config = :config),
    ADTypes.AutoPolyesterForwardDiff(),
    ADTypes.AutoPolyesterForwardDiff(chunksize = 3, tag = :tag),
    ADTypes.AutoReverseDiff(),
    ADTypes.AutoReverseDiff(compile = true),
    ADTypes.AutoSymbolics(),
    ADTypes.AutoTapir(),
    ADTypes.AutoTapir(safe_mode = false),
    ADTypes.AutoTracker(),
    ADTypes.AutoZygote(),
    # sparse
    ADTypes.AutoSparse(ADTypes.AutoForwardDiff()),
    ADTypes.AutoSparse(
        ADTypes.AutoForwardDiff();
        sparsity_detector = FakeSparsityDetector(),
        coloring_algorithm = FakeColoringAlgorithm()
    )
]
    println(backend)
end

using Setfield

@testset "Setfield compatibility" begin
    ad = AutoEnzyme()
    @test ad.mode === nothing
    @set! ad.mode = EnzymeCore.Reverse
    @test ad.mode isa EnzymeCore.ReverseMode

    struct CustomTestTag end

    ad = AutoForwardDiff()
    @test ad.tag === nothing
    @set! ad.tag = CustomTestTag()
    @test ad.tag isa CustomTestTag

    ad = AutoForwardDiff(; chunksize = 10)
    @test ad.tag === nothing
    @set! ad.tag = CustomTestTag()
    @test ad.tag isa CustomTestTag
    @test ad isa AutoForwardDiff{10}

    ad = AutoPolyesterForwardDiff()
    @test ad.tag === nothing
    @set! ad.tag = CustomTestTag()
    @test ad.tag isa CustomTestTag

    ad = AutoPolyesterForwardDiff(; chunksize = 10)
    @test ad.tag === nothing
    @set! ad.tag = CustomTestTag()
    @test ad.tag isa CustomTestTag
    @test ad isa AutoPolyesterForwardDiff{10}
end

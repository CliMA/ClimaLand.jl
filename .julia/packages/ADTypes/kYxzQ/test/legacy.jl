@testset "AutoModelingToolkit" begin
    ad_sparse1 = @test_deprecated AutoModelingToolkit(;
        obj_sparse = true, cons_sparse = false)
    ad_sparse2 = @test_deprecated AutoModelingToolkit(true, false)

    ad_dense1 = @test_deprecated AutoModelingToolkit(;
        obj_sparse = false, cons_sparse = false)
    ad_dense2 = @test_deprecated AutoModelingToolkit(false, false)
    ad_dense3 = @test_deprecated AutoModelingToolkit()

    @test all(
        isa.((ad_sparse1, ad_sparse2, ad_dense1, ad_dense2, ad_dense3), AbstractADType))
    @test all(isa.((ad_sparse1, ad_sparse2), AutoSparse{<:AutoSymbolics}))
    @test all(isa.((ad_dense1, ad_dense2, ad_dense3), AutoSymbolics))
end

@testset "AutoSparseFastDifferentiation" begin
    ad = @test_deprecated AutoSparseFastDifferentiation()
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoFastDifferentiation
end

@testset "AutoSparseFiniteDiff" begin
    ad = @test_deprecated AutoSparseFiniteDiff()
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoFiniteDiff
    @test dense_ad(ad).fdtype === Val(:forward)
    @test dense_ad(ad).fdjtype === Val(:forward)
    @test dense_ad(ad).fdhtype === Val(:hcentral)
end

@testset "AutoSparseForwardDiff" begin
    ad = @test_deprecated AutoSparseForwardDiff()
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoForwardDiff{nothing, Nothing}

    ad = @test_deprecated AutoSparseForwardDiff(; chunksize = 10, tag = CustomTag())
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoForwardDiff{10, CustomTag}
end

@testset "AutoSparsePolyesterForwardDiff" begin
    ad = @test_deprecated AutoSparsePolyesterForwardDiff(;
        chunksize = 10, tag = CustomTag())
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoPolyesterForwardDiff{10, CustomTag}
end

@testset "AutoSparseReverseDiff" begin
    ad = @test_deprecated AutoSparseReverseDiff(; compile = true)
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoReverseDiff
    @test dense_ad(ad).compile

    ad = @test_deprecated AutoSparseReverseDiff(true)
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoReverseDiff
    @test dense_ad(ad).compile
end

@testset "AutoSparseZygote" begin
    ad = @test_deprecated AutoSparseZygote()
    @test ad isa AbstractADType
    @test dense_ad(ad) isa AutoZygote
end

@testset "AutoReverseDiff without kwarg" begin
    ad = @test_deprecated AutoReverseDiff(true)
    @test ad.compile
end

@testset "AutoTapir" begin
    @test_deprecated AutoTapir()
    @test_deprecated AutoTapir(; safe_mode = false)
end

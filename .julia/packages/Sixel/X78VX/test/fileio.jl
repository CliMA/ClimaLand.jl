@testset "fileio" begin
    tmp_file = tempname() * ".sixel"
    for img in (
        repeat(Gray.(0:0.1:0.9), inner=(10, 50)),
        repeat(RGB.(0:0.1:0.9), inner=(10, 50))
    )
        save(tmp_file, img)
        img_readback = load(tmp_file)
        @test eltype(img_readback) == RGB{N0f8}
        @test img_readback isa IndirectArray

        @test assess_psnr(img, eltype(img).(img_readback)) > 30
    end
end

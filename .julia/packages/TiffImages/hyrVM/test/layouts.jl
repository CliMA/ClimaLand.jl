@testset "Interpretations" begin
    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.BITSPERSAMPLE] = [16]
    ifd[TiffImages.SAMPLEFORMAT] = [3]
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.SAMPLESPERPIXEL] = 1

    @test TiffImages.interpretation(ifd) == Gray{Float16}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.SAMPLESPERPIXEL] = 2
    ifd[TiffImages.BITSPERSAMPLE] = [32, 32]
    ifd[TiffImages.SAMPLEFORMAT] = [3, 3]

    @test TiffImages.interpretation(ifd) == TiffImages.WidePixel{Gray{Float32}, Tuple{Float32}}

    # handle case where SamplesPerPixel is missing, issue #56
    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.BITSPERSAMPLE] = [8]

    @test TiffImages.interpretation(ifd) == Gray{N0f8}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.SAMPLESPERPIXEL] = 2
    ifd[TiffImages.BITSPERSAMPLE] = [8, 8]
    ifd[TiffImages.SAMPLEFORMAT] = [2, 2]
    ifd[TiffImages.EXTRASAMPLES] = TiffImages.EXTRASAMPLE_ASSOCALPHA

    @test TiffImages.interpretation(ifd) == GrayA{Q0f7}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.SAMPLESPERPIXEL] = 2
    ifd[TiffImages.BITSPERSAMPLE] = [8, 8]
    ifd[TiffImages.SAMPLEFORMAT] = [2, 2]
    ifd[TiffImages.EXTRASAMPLES] = TiffImages.EXTRASAMPLE_UNASSALPHA

    @test TiffImages.interpretation(ifd) == TiffImages.WidePixel{Gray{Q0f7}, Tuple{Q0f7}}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.SAMPLESPERPIXEL] = 2
    ifd[TiffImages.BITSPERSAMPLE] = [8, 8]
    ifd[TiffImages.SAMPLEFORMAT] = [1, 1]
    ifd[TiffImages.EXTRASAMPLES] = TiffImages.EXTRASAMPLE_UNSPECIFIED

    @test TiffImages.interpretation(ifd) == TiffImages.WidePixel{Gray{N0f8}, Tuple{N0f8}}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.SAMPLESPERPIXEL] = 3
    ifd[TiffImages.BITSPERSAMPLE] = [12, 13, 14]
    ifd[TiffImages.SAMPLEFORMAT] = [1, 1, 1]
    ifd[TiffImages.EXTRASAMPLES] = [0, 0]

    @test TiffImages.interpretation(ifd) == TiffImages.WidePixel{Gray{N4f12}, Tuple{N3f13, N2f14}}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.SAMPLESPERPIXEL] = 3
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_RGB
    ifd[TiffImages.SAMPLEFORMAT] = [1, 1, 1]
    ifd[TiffImages.BITSPERSAMPLE] = [12, 12, 12]

    @test TiffImages.interpretation(ifd) == RGB{N4f12}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_RGB
    ifd[TiffImages.SAMPLESPERPIXEL] = 4
    ifd[TiffImages.BITSPERSAMPLE] = [8, 8, 8, 8]
    ifd[TiffImages.SAMPLEFORMAT] = [1, 1, 1, 1]
    ifd[TiffImages.EXTRASAMPLES] = TiffImages.EXTRASAMPLE_UNSPECIFIED

    @test TiffImages.interpretation(ifd) == TiffImages.WidePixel{RGB{N0f8}, Tuple{N0f8}}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_RGB
    ifd[TiffImages.SAMPLESPERPIXEL] = 4
    ifd[TiffImages.BITSPERSAMPLE] = [8, 8, 8, 8]
    ifd[TiffImages.SAMPLEFORMAT] = [1, 1, 1, 1]
    ifd[TiffImages.EXTRASAMPLES] = TiffImages.EXTRASAMPLE_ASSOCALPHA

    @test TiffImages.interpretation(ifd) == RGBA{N0f8}

    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_RGB
    ifd[TiffImages.SAMPLESPERPIXEL] = 6
    ifd[TiffImages.BITSPERSAMPLE] = [16, 16, 16, 16, 16, 16]
    ifd[TiffImages.SAMPLEFORMAT] = [1, 1, 1, 1, 1, 1]
    ifd[TiffImages.EXTRASAMPLES] = [1, 0, 2]

    @test TiffImages.interpretation(ifd) == TiffImages.WidePixel{RGBA{N0f16}, Tuple{N0f16, N0f16}}
end

@testset "Sample Types" begin
    @test TiffImages.bitspersample(RGB{Float32}) == [32, 32, 32]
    @test TiffImages.bitspersample(RGB{N1f7}) == [7, 7, 7]
    @test TiffImages.bitspersample(RGB{Q4f11}) == [12, 12, 12]
    @test TiffImages.bitspersample(TiffImages.WidePixel{RGB{N0f16}, Tuple{N0f16, N0f16}}) == [16, 16, 16, 16 ,16]

    @test TiffImages.sampleformat(RGB{Float32}) == [3, 3, 3]
    @test TiffImages.sampleformat(RGB{N1f7}) == [1, 1, 1]
    @test TiffImages.sampleformat(RGB{Q4f11}) == [2, 2, 2]
    @test TiffImages.sampleformat(TiffImages.WidePixel{RGB{N0f16}, Tuple{N0f16, N0f16}}) == [1, 1, 1, 1, 1]

    @test TiffImages.extrasamples(RGBA{Float64}) == 1
    @test TiffImages.extrasamples(RGB{Float64}) == nothing
    @test TiffImages.extrasamples(TiffImages.WidePixel{RGBA{Float64}, Tuple{Float64, Float64}}) == [1, 0, 0]
    @test TiffImages.extrasamples(TiffImages.WidePixel{RGB{Float64}, Tuple{Float64, Float64}}) == [0, 0]
end

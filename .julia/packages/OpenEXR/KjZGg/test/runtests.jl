using Test, OpenEXR, FileIO, Colors, FixedPointNumbers

# helpers for testing channelwise bit-depth truncation
function channelwise_max_diff(a::RGBA, b::RGBA)
    max(abs(a.r - b.r), abs(a.g - b.g), abs(a.b - b.b), abs(a.alpha - b.alpha))
end

function channelwise_max_diff(a::RGB, b::RGB)
    max(abs(a.r - b.r), abs(a.g - b.g), abs(a.b - b.b))
end

function channelwise_max_diff(a::GrayA, b::GrayA)
    max(abs(a.val - b.val), abs(a.alpha - b.alpha))
end

function channelwise_max_diff(a::Gray, b::Gray)
    abs(a.val - b.val)
end

@testset "RoundTrip" begin

    @testset "Identical" begin
        for typ in (RGBA{Float16}, RGB{Float16}, GrayA{Float16}, Gray{Float16})
            img = rand(typ, 256, 512)
            fn = File{DataFormat{:EXR}}(tempname())
            OpenEXR.save(fn, img)
            try
                loaded_img = OpenEXR.load(fn)
                @test typeof(loaded_img) === Array{typ,2}
                @test img == loaded_img
            finally
                rm(fn.filename)
            end
        end
    end

    @testset "Lossless with type conversion" begin
        for (save_type, load_type) in (
            (RGBA{N0f8}, RGBA{Float16}),
            (RGB{N0f8}, RGB{Float16}),
            (GrayA{N0f8}, GrayA{Float16}),
            (Gray{N0f8}, Gray{Float16}),
        )
            img = rand(save_type, 256, 512)
            fn = File{DataFormat{:EXR}}(tempname())
            OpenEXR.save(fn, img)
            try
                loaded_img = OpenEXR.load(fn)
                @test typeof(loaded_img) === Array{load_type,2}
                converted_img = (c -> convert(save_type, c)).(loaded_img)
                @test converted_img == img
            finally
                rm(fn.filename)
            end
        end
    end

    @testset "Lossy with type conversion" begin

        @testset "Bit depth truncation" begin
            for (save_type, load_type) in (
                (RGBA{Float32}, RGBA{Float16}),
                (RGB{Float32}, RGB{Float16}),
                (GrayA{Float32}, GrayA{Float16}),
                (Gray{Float32}, Gray{Float16}),
            )
                img = rand(save_type, 256, 512)
                fn = File{DataFormat{:EXR}}(tempname())
                OpenEXR.save(fn, img)
                try
                    loaded_img = OpenEXR.load(fn)

                    # return type is as expected
                    @test typeof(loaded_img) === Array{load_type,2}

                    # differences are bounded channelwise
                    diffs = map(channelwise_max_diff, img, loaded_img)
                    @test maximum(diffs) <= 2.5e-4
                finally
                    rm(fn.filename)
                end
            end
        end

        @testset "Colorspace conversion" begin
            for save_type in
                (HSV{Float16}, HSL{Float16}, Lab{Float16}, LCHab{Float16}, YIQ{Float16})
                img = rand(save_type, 256, 512)
                fn = File{DataFormat{:EXR}}(tempname())
                OpenEXR.save(fn, img)
                try
                    loaded_img = OpenEXR.load(fn)

                    # return type is as expected
                    @test typeof(loaded_img) === Array{RGB{Float16},2}

                    # differences are bounded channelwise
                    converted_img = (c -> convert(RGB{Float16}, c)).(img)
                    diffs = map(channelwise_max_diff, converted_img, loaded_img)
                    @test maximum(diffs) <= 2.5e-4
                finally
                    rm(fn.filename)
                end
            end
        end
    end
end

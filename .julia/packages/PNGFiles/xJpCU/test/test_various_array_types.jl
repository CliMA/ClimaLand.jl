import AxisArrays: AxisArray
using PaddedViews
using OffsetArrays

@testset "Non-Standard Arrays" begin
    function _load_and_save(io, img)
        PNGFiles.save(io, img)
        seek(io, 0)
        PNGFiles.load(io)
    end

    @testset "OffsetArrays" begin
        for CT in (Gray{N0f8}, RGB{N0f8}, Gray{Float32}, RGB{Float32})
            img = OffsetArray(rand(CT, 4, 4), -1, -1)

            ref = _load_and_save(IOBuffer(), parent(img))
            actual = _load_and_save(IOBuffer(), img)
            @test ref == actual
            open(joinpath(PNG_TEST_PATH, "test_offsetarray.png"), "w+") do fio
                ref = _load_and_save(fio, parent(img))
                actual = _load_and_save(fio, img)
                @test ref == actual
            end
        end
    end

    @testset "AxisArrays" begin
        for CT in (Gray{N0f8}, RGB{N0f8}, Gray{Float32}, RGB{Float32})
            img = AxisArray(rand(CT, 4, 4))

            ref = _load_and_save(IOBuffer(), convert(Array, img))
            actual = _load_and_save(IOBuffer(), img)
            @test ref == actual
            open(joinpath(PNG_TEST_PATH, "test_axisarray.png"), "w+") do fio
                ref = _load_and_save(fio, convert(Array, img))
                actual = _load_and_save(fio, img)
                @test ref == actual
            end
        end
    end

    @testset "PaddedViews" begin
        for CT in (Gray{N0f8}, RGB{N0f8}, Gray{Float32}, RGB{Float32})
            img = PaddedView(zero(CT), rand(CT, 4, 4), (-1:5, -1:5))
            ref = _load_and_save(IOBuffer(), collect(img))
            actual = _load_and_save(IOBuffer(), img)
            @test ref == actual
            open(joinpath(PNG_TEST_PATH, "test_paddedviews.png"), "w+") do fio
                ref = _load_and_save(fio, collect(img))
                actual = _load_and_save(fio, img)
                @test ref == actual
            end
        end
    end

    @testset "indexed image with special palatte type" begin
        index = rand(1:5, 3, 4)
        palatte = colorview(RGB{N0f8}, rand(UInt8, 3, 5))
        img = IndirectArray(index, palatte)
        ref = _load_and_save(IOBuffer(), collect(img))
        actual = _load_and_save(IOBuffer(), img)
        @test ref == actual
        @test actual isa IndirectArray{RGB{N0f8}, 2}
    end
end

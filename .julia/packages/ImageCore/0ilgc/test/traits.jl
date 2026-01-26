using ImageCore, Colors, FixedPointNumbers, ColorVectorSpace, MappedArrays, OffsetArrays
using Test
using ImageCore: Pixel, NumberLike, GenericImage, GenericGrayImage, default_names

@testset "Image traits" begin
    for (B, swap) in ((rand(UInt16(1):UInt16(20), 3, 5), false),
                      (rand(Gray{Float32}, 3, 5), false),
                      (rand(RGB{Float16}, 3, 5), false),
                      (bitrand(3, 5), false),
                      (rand(UInt32, 3, 5), false),
                      (view(rand(3, 2, 5), :, 1, :), false),
                      (OffsetArray(rand(3, 5), -1:1, -2:2), false),
                      (PermutedDimsArray(rand(5, 3), (2, 1)), true),
                      (mappedarray(identity, PermutedDimsArray(rand(5, 3), (2, 1))), true),
                      (colorview(RGB, zeros(3, 5), zeroarray, zeros(3, 5)), false))
        @test pixelspacing(B) == (1,1)
        if !isa(B, SubArray)
            @test spacedirections(B) == (swap ? ((0,1),(1,0)) : ((1,0),(0,1)))
        else
            @test spacedirections(B) == ((1,0,0), (0,0,1))
        end
        @test sdims(B) == 2
        @test coords_spatial(B) == (swap ? (2,1) : (1,2))
        @test nimages(B) == 1
        @test size_spatial(B) == (3,5)
        if isa(B, OffsetArray)
            @test indices_spatial(B) == (-1:1, -2:2)
        else
            @test indices_spatial(B) == (Base.OneTo(3), Base.OneTo(5))
        end
        assert_timedim_last(B)
    end
end

# delibrately written in a redundant way
@testset "*Like traits" begin
    @testset "Pixel" begin
        @test NumberLike <: Pixel
        @test Number <: Pixel
        @test Gray <: Pixel
        @test RGB <: Pixel

        @test isa(oneunit(Gray), Pixel)
        @test isa(RGB(1.0, 0.0, 0.0), Pixel)
    end

    @testset "NumberLike" begin
        @test Number <: NumberLike
        @test Real <: NumberLike
        @test AbstractFloat <: NumberLike
        @test FixedPoint <: NumberLike
        @test Integer <: NumberLike
        @test Bool <: NumberLike

        @test Gray <: NumberLike
        @test Gray{<:AbstractFloat} <: NumberLike
        @test Gray{<:Bool} <: NumberLike

        @test isa(oneunit(Gray), NumberLike)
    end

    @testset "GenericImage" begin
        @test GenericGrayImage <: GenericImage
        for sz in Any[(3, 3), (3, 3, 3)]
            @test isa(rand(Bool, sz), GenericImage)
            @test isa(rand(N0f8, sz), GenericImage)
            @test isa(rand(Float32, sz), GenericImage)

            @test isa(rand(Gray, sz), GenericImage)
            @test isa(rand(Gray{Bool}, sz), GenericImage)
            @test isa(rand(Gray{N0f8}, sz), GenericImage)
            @test isa(rand(Gray{Float32}, sz), GenericImage)

            @test isa(rand(GrayA, sz), GenericImage)
            @test isa(rand(GrayA{N0f8}, sz), GenericImage)
            @test isa(rand(GrayA{Float32}, sz), GenericImage)

            @test isa(rand(RGB, sz), GenericImage)
            @test isa(rand(RGB{N0f8}, sz), GenericImage)
            @test isa(rand(RGB{Float32}, sz), GenericImage)

            @test isa(rand(RGBA, sz), GenericImage)
            @test isa(rand(RGBA{N0f8}, sz), GenericImage)
            @test isa(rand(RGBA{Float32}, sz), GenericImage)

            @test isa(rand(Lab, sz), GenericImage)
            @test isa(rand(Lab{Float32}, sz), GenericImage)
        end
    end

    @testset "GrayImage" begin
        for sz in [(3, 3), (3, 3, 3)]
            @test isa(rand(Bool, sz), GenericGrayImage)
            @test isa(rand(N0f8, sz), GenericGrayImage)
            @test isa(rand(Float32, sz), GenericGrayImage)

            @test isa(rand(Gray, sz), GenericGrayImage)
            @test isa(rand(Gray{Bool}, sz), GenericGrayImage)
            @test isa(rand(Gray{N0f8}, sz), GenericGrayImage)
            @test isa(rand(Gray{Float32}, sz), GenericGrayImage)
        end
    end

    @testset "dispatch" begin
        begin
            whatis(::GenericImage) = "GenericImage"
            whatis(::GenericGrayImage) = "GenericGrayImage"
            whatis(::GenericImage{<:AbstractRGB}) = "GenericRGBImage"

            whatis(::GenericImage{<:Pixel, 2}) = "Generic2dImage"
            whatis(::GenericGrayImage{<:NumberLike, 2}) = "Gray2dImage"
            whatis(::GenericImage{<:AbstractRGB, 2}) = "RGB2dImage"

            @test whatis(rand(Lab, 2, 2, 2)) == "GenericImage"

            @test whatis(rand(Lab, 2, 2)) == "Generic2dImage"

            @test whatis(rand(2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(N0f8, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Bool, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Float32, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Int64, 2, 2, 2)) == "GenericGrayImage"

            @test whatis(rand(Gray, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Gray{N0f8}, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Gray{Bool}, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Gray{Float32}, 2, 2, 2)) == "GenericGrayImage"

            @test whatis(rand(2, 2)) == "Gray2dImage"
            @test whatis(rand(N0f8, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Bool, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Float32, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Int64, 2, 2)) == "Gray2dImage"

            @test whatis(rand(Gray, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Gray{N0f8}, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Gray{Bool}, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Gray{Float32}, 2, 2)) == "Gray2dImage"

            @test whatis(rand(RGB, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(RGB{N0f8}, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(RGB{Float32}, 2, 2, 2)) == "GenericRGBImage"

            @test whatis(rand(BGR, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(BGR{N0f8}, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(BGR{Float32}, 2, 2, 2)) == "GenericRGBImage"

            @test whatis(rand(RGB, 2, 2)) == "RGB2dImage"
            @test whatis(rand(RGB{N0f8}, 2, 2)) == "RGB2dImage"
            @test whatis(rand(RGB{Float32}, 2, 2)) == "RGB2dImage"

            @test whatis(rand(BGR, 2, 2)) == "RGB2dImage"
            @test whatis(rand(BGR{N0f8}, 2, 2)) == "RGB2dImage"
            @test whatis(rand(BGR{Float32}, 2, 2)) == "RGB2dImage"
        end
    end
end

struct RowVector{T,P} <: AbstractVector{T}
    v::Vector{T}
    p::P
end

ImageCore.HasDimNames(::Type{<:RowVector}) = HasDimNames{true}()

ImageCore.HasProperties(::Type{<:RowVector}) = HasProperties{true}()

Base.names(::RowVector) = (:row,)
Base.axes(rv::RowVector) = axes(rv.v)


@testset "Trait Interface" begin
    img = reshape(1:24, 2,3,4)
    @test @inferred(namedaxes(img)) == NamedTuple{(:dim_1, :dim_2, :dim_3)}(axes(img))
    @test @inferred(HasDimNames(img)) == HasDimNames{false}()
    @test @inferred(HasProperties(img)) == HasProperties{false}()

    rv = RowVector([1:10...], Dict{String,Any}())
    @test @inferred(HasDimNames(rv)) == HasDimNames{true}()
    @test @inferred(HasProperties(rv)) == HasProperties{true}()
    @test @inferred(namedaxes(rv)) == NamedTuple{(:row,)}((Base.OneTo(10),))

    # default names
    @test @inferred(default_names(Val(3))) == (:dim_1, :dim_2, :dim_3)
end

nothing

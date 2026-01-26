@info "Beginning deprecation tests"

@testset "Deprecations" begin
    @testset "copy/similar" begin
        img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3])
        img2 = copy(img)
        @test img2.data == img.data
        img2[2,2] = -1
        @test img2[2,2] < 0
        @test img[2,2] >= 0
        img2["prop2"][2] = -2
        @test img2["prop2"] == [1,-2,3]
        @test img["prop2"] == [1,2,3]
        img2 = similar(img)
        @test img2["prop1"] == 1
        @test img2["prop2"] == [1,2,3]
        @test img2.data != img.data
        img2["prop3"] = 7
        @test !haskey(img, "prop3")
        img2 = similar(img, RGB{Float16}, (Base.OneTo(5),))
        @test img2["prop1"] == 1
        @test img2["prop2"] == [1,2,3]
        @test eltype(img2) == RGB{Float16}
        @test size(img2) == (5,)
        @test img2.data != img.data
        img2["prop3"] = 7
        @test !haskey(img, "prop3")

        A = AxisArray(rand(3,5), :y, :x)
        B = ImageMeta(A, info="blah")
        C = similar(B)
        @test isa(C, ImageMeta)
        @test isa(C.data, AxisArray)
        @test eltype(C) == Float64
        C = similar(B, RGB{Float16})
        @test isa(C, ImageMeta)
        @test isa(C.data, AxisArray)
        @test eltype(C) == RGB{Float16}
    end

    @testset "copy/shareproperties/viewim" begin
        img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3])
        @test !isempty(properties(img))
        v = view(img, 1:2, 1:2)
        c = img[1:2, 1:2]
        @test v["prop1"] == 1
        @test c["prop1"] == 1
        img2 = copyproperties(img, reshape(1:15, 5, 3))
        @test size(img2) == (5,3)
        img2["prop1"] = -1
        @test img["prop1"] == 1
        img2 = shareproperties(img, reshape(1:15, 5, 3))
        @test size(img2) == (5,3)
        img2["prop1"] = -1
        @test img["prop1"] == -1
        @test v["prop1"] == -1
        @test c["prop1"] == 1
        imgb = ImageMeta(rand(RGB{N0f8}, 2, 2), propa = "hello", propb = [1,2])
        copy!(img, imgb, "propa", "propb")
        @test img["propa"] == "hello"
        @test img["propb"] == [1,2]
        img["propb"][2] = 10
        @test img["propb"] == [1,10]
        @test imgb["propb"] == [1,2]
        delete!(img, "propb")
        @test  haskey(img, "propa")
        @test !haskey(img, "propb")

        img = ImageMeta(AxisArray(rand(3,5,8),
                                  Axis{:x}(1:3),
                                  Axis{:y}(1:5),
                                  Axis{:time}(0.1:0.1:0.8));
                        prop1 = 1, prop2 = [1,2,3])
        v = view(img, Axis{:time}(0.25..0.5))
        c = img[Axis{:time}(0.25..0.5)]
        @test v["prop1"] == 1
        @test c["prop1"] == 1
    end

    @testset "traits" begin
        img = ImageMeta(AxisArray(rand(3,5,8),
                                  Axis{:x}(1:3),
                                  Axis{:y}(1:5),
                                  Axis{:time}(0.1:0.1:0.8)))
        @test @inferred(timedim(img)) == 3
        @test @inferred(pixelspacing(img)) === (1,1)
        @test spacedirections(img) === ((1,0),(0,1))
        img["spacedirections"] = "broken"
        @test spacedirections(img) == "broken"
        @test sdims(img) == 2
        @test @inferred(coords_spatial(img)) == (1,2)
        @test spatialorder(img) == (:x, :y)
        @test nimages(img) == 8
        @test @inferred(size_spatial(img)) == (3,5)
        @test @inferred(indices_spatial(img)) == (Base.OneTo(3), Base.OneTo(5))
        assert_timedim_last(img)
        @test timeaxis(img) == Axis{:time}(0.1:0.1:0.8)
    end

    @testset "spatialprops" begin
        img = ImageMeta(rand(3,5),
                        spatialproperties=Set(["vector","matrix","tuple"]),
                        vector=[1,2],
                        matrix=[1 3; 2 4],
                        tuple=(1,2))
        for imgp in (img', permutedims(img, (2,1)), PermutedDimsArray(img, (2,1)))
            @test imgp.data == img.data'
            @test imgp["vector"] == [2,1]
            @test imgp["matrix"] == [4 2; 3 1]
            @test imgp["tuple"]  == (2,1)
            @test img["vector"] == [1,2]
            @test img["matrix"] == [1 3; 2 4]
            @test img["tuple"]  == (1,2)
        end
    end

    @testset "AbstractDicts" begin
        img = ImageMeta(AxisArray(rand(3,5,8),
                                  Axis{:x}(1:3),
                                  Axis{:y}(1:5),
                                  Axis{:time}(0.1:0.1:0.8)),
                        IdDict{String,Any}("pixelspacing" => (1,1)))
        @test @inferred(pixelspacing(img)) === (1,1)

        rand_img = rand(10,10)
        dict = Dict("attr1" => 3, "attr2" => 4)
        @test ImageMeta(rand_img, dict)[1, 1] == rand_img[1, 1]
    end
end

nothing

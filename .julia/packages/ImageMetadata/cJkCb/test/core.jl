using ImageCore, SimpleTraits, ImageAxes, ImageMetadata, OffsetArrays, ImageBase
using Test
import Dates: now
using Unitful: m
import AxisArrays
using AxisArrays: AxisArray, axisnames, (..)

@testset "indexing" begin
    # 1d images
    for A = (rand(3),
             view(rand(6), 1:2:5),
             view(rand(6), [1,2,4]),
             rand(RGB{N0f8}, 3),
             rand(Gray{N0f8}, 3))
        img = ImageMeta(A, prop1=1, prop2=[1,2,3])
        @test eltype(img) == eltype(A)
        @test IndexStyle(img) == IndexStyle(A)
        @test ndims(img) == 1
        @test size(img) == (3,)
        @test arraydata(img) === A
        for i = 1:3
            @test @inferred(img[i]) === A[i]
        end
        for I in eachindex(img)
            @test @inferred(img[I]) === A[I]
        end
        k = 0
        for a in img
            @test a == A[k+=1]
        end
        img[2] = zero(eltype(img))
        @test A[2] == zero(eltype(A))
        img[3] = oneunit(eltype(img))
        @test A[3] == oneunit(eltype(A))
        @test_throws BoundsError img[0]
        @test_throws BoundsError img[4]
        @test img.prop1 == 1
        @test img.prop2 == [1,2,3]
        img.prop1 = -1
        @test img.prop1 == -1
    end
    for A in (rand(3,5),
              view(rand(4,6), 1:3, 1:5),
              view(rand(4,5), [1,2,4], :),
              reshape(1:15, 3, 5),
              rand(RGB{N0f8}, 3, 5),
              rand(Gray{Float32}, 3, 5))
        img = ImageMeta(A, prop1=1, prop2=[1,2,3]) # TODO: add @inferred (see julia #17719)
        @test eltype(img) == eltype(A)
        @test IndexStyle(img) == IndexStyle(A)
        @test ndims(img) == 2
        @test size(img) == (3,5)
        @test arraydata(img) === A
        for j = 1:5, i = 1:3
            @test @inferred(img[i,j]) === A[i,j]
        end
        for k = 1:15
            @test @inferred(img[k]) === A[k]
        end
        for I in eachindex(img)
            @test @inferred(img[I]) === A[I]
        end
        k = 0
        for a in img
            @test a == A[k+=1]
        end
        if !isa(A, typeof(reshape(1:15, 3, 5)))
            img[2,3] = zero(eltype(img))
            @test A[2,3] == zero(eltype(A))
            img[4] = oneunit(eltype(img))
            @test A[4] == oneunit(eltype(A))
        end
        @test_throws BoundsError img[0,0]
        @test_throws BoundsError img[4,1]
        @test_throws BoundsError img[1,6]
        @test img.prop1 == 1
        @test img.prop2 == [1,2,3]
        img.prop1 = -1
        @test img.prop1 == -1
        # vector-indexing
        @test isa(@inferred(img[:,:]), ImageMeta) && img[:,:] == img
        @test isa(@inferred(img[:]), ImageMeta) && img[:] == A[:]
        @test isa(@inferred(img[1:2,1:2]), ImageMeta) && img[1:2,1:2] == A[1:2,1:2]
        @test isa(@inferred(img[1:2,:]), ImageMeta) && img[1:2,:] == A[1:2,:]
        @test isa(@inferred(img[:,1:2]), ImageMeta) && img[:,1:2] == A[:,1:2]
        @test isa(@inferred(view(img,:,:)), ImageMeta) && view(img,:,:) == img
        @test isa(@inferred(view(img,:)), ImageMeta) && view(img,:) == A[:]
        @test isa(@inferred(view(img,1:2,1:2)), ImageMeta) && view(img,1:2,1:2) == A[1:2,1:2]
        @test isa(@inferred(view(img,1:2,:)), ImageMeta) && view(img,1:2,:) == A[1:2,:]
        @test isa(@inferred(view(img,:,1:2)), ImageMeta) && view(img,:,1:2) == A[:,1:2]
    end
    # Test bounds-checking removal by @inbounds
    if Base.JLOptions().check_bounds != 1 && Base.JLOptions().can_inline == 1
        set5!(x) = @inbounds x[5] = 1.234
        get5(x) = @inbounds x[5]
        aa = zeros(3)
        sizehint!(aa, 10)  # make sure we don't cause a segfault
        a = ImageMeta(aa)
        @test_throws BoundsError a[5]
        set5!(a)
        val = get5(a)
        @test val == 1.234
        a = zeros(3,5)
    end

    # Issue #10
    A = AxisArray(rand(3,5), :y, :x)
    B = ImageMeta(A, info="blah")
    Broi = B[2:3, 2:3]
    @test isa(Broi, ImageMeta)
    @test axisnames(Broi) == (:y, :x)
    A1, B1 = A[2:7], B[2:7]
    @test isa(B1, ImageMeta) && A1 == B1
    Broi = view(B, 2:3, 2:3)
    @test isa(Broi, ImageMeta)
    @test axisnames(Broi) == (:y, :x)

    # Non-1 indexing
    Ao = OffsetArray(rand(2,3), 0:1, -1:1)
    A = ImageMeta(Ao, info="blah")
    @test axes(A) === axes(Ao)
end

@testset "convert" begin
    A = rand(3,5)
    M = convert(ImageMeta, A)
    @test isa(M, ImageMeta)
    @test convert(ImageMeta, M) === M
    @test convert(ImageMeta{Float64}, M) === M
    @test eltype(convert(ImageMeta{Gray{Float64}}, M)) == Gray{Float64}
    metagray = Gray.(M)
    @test eltype(metagray) == Gray{Float64} && isa(metagray, ImageMeta)
    @test convert(ImageMeta{Gray{Float64}}, M) == M
    @test convert(ImageMeta{Gray{Float64}}, A) == A
end

@testset "reinterpret" begin
    # It's possible that reinterpretc shouldn't be defined for ImageMeta, but...
    A = rand(Float32, 4, 5)
    M = ImageMeta(A, meta=true)
    Mr = reinterpretc(Gray, M)
    @test eltype(Mr) == Gray{Float32}
    @test Mr.meta = true
    # Ensure that it gets defined for the formerly un-reinterpretable
    A = zeros(Float32, 4, 5)
    M = ImageMeta(view(A, 1:2:3, 1:4), meta=true)
    Mr = reinterpretc(Gray, M)
    Mr = reinterpretc(Gray{Float32}, M)
    @test eltype(Mr) == Gray{Float32}
    @test Mr.meta = true
    Mr[:] .= Gray{Float32}(0.5)
    @test all(A[1:2:3,1:4] .== 0.5)
end

@testset "copy/similar" begin
    img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3])
    img2 = copy(img)
    @test arraydata(img2) == arraydata(img)
    img2[2,2] = -1
    @test img2[2,2] < 0
    @test img[2,2] >= 0
    img2.prop2[2] = -2
    @test img2.prop2 == [1,-2,3]
    @test img.prop2 == [1,2,3]
    img2 = similar(img)
    @test img2.prop1 == 1
    @test img2.prop2 == [1,2,3]
    @test arraydata(img2) != arraydata(img)
    img2.prop3 = 7
    @test !hasproperty(img, :prop3)
    img2 = similar(img, RGB{Float16}, (Base.OneTo(5),))
    @test img2.prop1 == 1
    @test img2.prop2 == [1,2,3]
    @test eltype(img2) == RGB{Float16}
    @test size(img2) == (5,)
    @test arraydata(img2) != arraydata(img)
    img2.prop3 = 7
    @test !hasproperty(img, :prop3)

    A = AxisArray(rand(3,5), :y, :x)
    B = ImageMeta(A, info="blah")
    C = similar(B)
    @test isa(C, ImageMeta)
    @test isa(arraydata(C), AxisArray)
    @test eltype(C) == Float64
    C = similar(B, RGB{Float16})
    @test isa(C, ImageMeta)
    @test isa(arraydata(C), AxisArray)
    @test eltype(C) == RGB{Float16}
end

@testset "copy/shareproperties/viewim" begin
    img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3])
    @test !isempty(properties(img))
    v = view(img, 1:2, 1:2)
    c = img[1:2, 1:2]
    @test v.prop1 == 1
    @test c.prop1 == 1
    img2 = copyproperties(img, reshape(1:15, 5, 3))
    @test size(img2) == (5,3)
    img2.prop1 = -1
    @test img.prop1 == 1
    img2 = shareproperties(img, reshape(1:15, 5, 3))
    @test size(img2) == (5,3)
    img2.prop1 = -1
    @test img.prop1 == -1
    @test v.prop1 == -1
    @test c.prop1 == 1
    imgb = ImageMeta(rand(RGB{N0f8}, 2, 2), propa = "hello", propb = [1,2])
    copy!(img, imgb, :propa, :propb)
    @test img.propa == "hello"
    @test img.propb == [1,2]
    img.propb[2] = 10
    @test img.propb == [1,10]
    @test imgb.propb == [1,2]
    delete!(img, :propb)
    @test  hasproperty(img, :propa)
    @test !hasproperty(img, :propb)

    img = ImageMeta(AxisArray(rand(3,5,8),
                              Axis{:x}(1:3),
                              Axis{:y}(1:5),
                              Axis{:time}(0.1:0.1:0.8));
                    prop1 = 1, prop2 = [1,2,3])
    v = view(img, Axis{:time}(0.25..0.5))
    c = img[Axis{:time}(0.25..0.5)]
    @test v.prop1 == 1
    @test c.prop1 == 1
end

@testset "views" begin
    for A1 in (rand(Gray{N0f8}, 4, 5), rand(RGB{Float32}, 4, 5))
        t1 = now()
        M1 = ImageMeta(A1, date=t1)
        vM1 = channelview(M1)
        @test isa(vM1, ImageMeta)
        @test vM1.date == t1
        @test arraydata(vM1) == channelview(A1)
        @test isa(rawview(vM1), ImageMeta)
        if ndims(rawview(vM1)) == 3
            @test rawview(vM1)[1,2,1] === rawview(channelview(A1))[1,2,1]
        else
            @test rawview(vM1)[1,2] === rawview(channelview(A1))[1,2]
        end
    end
    for (A1,C) in ((rand(UInt8, 4, 5), Gray), (rand(UInt8, 3, 4, 5), RGB))
        M1 = ImageMeta(A1)
        vM1 = normedview(M1)
        @test isa(vM1, ImageMeta)
        @test arraydata(vM1) == normedview(A1)
        cvM = colorview(C, vM1)
        @test isa(cvM, ImageMeta)
        @test cvM[1,2] === colorview(C, normedview(A1))[1,2]
    end
    A = AxisArray(rand(RGB{N0f8}, 3, 5), :y, :x)
    t = now()
    M = ImageMeta(A, date=t)
    vM = channelview(M)
    @test isa(vM, ImageMeta)
    @test colordim(vM) == 1
    pvM = permutedims(vM, (2,3,1))
    @test colordim(pvM) == 3
    @test pvM.date == t
    sleep(0.1)
    pvM.date = now()
    @test M.date == t && pvM.date != t
    vM = PermutedDimsArray(M, (2,1))
    @test vM.date == t
    @test axisnames(vM) == (:x, :y)
end

@testset "meta-axes" begin
    A = AxisArray(rand(3,5), :y, :x)
    M = ImageMeta(A)
    @test AxisArrays.HasAxes(M) == AxisArrays.HasAxes{true}()
    @test AxisArrays.axes(M) == (Axis{:y}(1:3), Axis{:x}(1:5))
    @test AxisArrays.axes(M, 2) == Axis{:x}(1:5)
    @test AxisArrays.axes(M, Axis{:y}) == Axis{:y}(1:3)
    @test AxisArrays.axes(M, Axis{:x}()) == Axis{:x}(1:5)
    @test axisdim(M, Axis{:y}) == 1
    @test size(M, Axis{:y}) == 3
    @test size(M, Axis{:x}()) == 5
    @test axes(M, Axis{:y}) === Base.OneTo(3)
    @test axes(M, Axis{:x}()) == Base.OneTo(5)
    @test axisdim(M, Axis{:x}) == 2
    @test_throws ErrorException axisdim(M, Axis{:z})
    @test axisnames(M) == (:y, :x)
    @test axisvalues(M) == (1:3, 1:5)
    @test M[Axis{:y}(2:3)] == A[Axis{:y}(2:3)]
    @test view(M, Axis{:y}(2:3)) == A[Axis{:y}(2:3)]
    M[Axis{:y}(2:3), Axis{:x}(1)] .= -5
    @test all(A[2:3,1] .== -5)
end

@testset "show" begin
    for supp in (Set([:prop3]), :prop3)
        img = ImageMeta(rand(3,5); prop1 = 1, prop2 = [1,2,3], suppress = supp, prop3 = "hide")
        str = string(img)
        @test occursin("ImageMeta with", str)
        @test occursin("prop1", str)
        @test occursin("prop2", str)
        @test occursin("$([1,2,3])", str)
        @test occursin("prop3", str)
        @test !occursin("hide", str)
        @test occursin("<suppressed>", str)
        @test !occursin("suppress:", str)
        io = IOBuffer()
        show(io, MIME("text/plain"), img)
        str2 = String(take!(io))
        @test str == str2
    end
end

@testset "traits" begin
    img = ImageMeta(AxisArray(rand(3,5,8),
                              Axis{:x}(1:3),
                              Axis{:y}(1:5),
                              Axis{:time}(0.1:0.1:0.8)))
    @test @inferred(timedim(img)) == 3
    @test @inferred(pixelspacing(img)) === (1,1)
    @test spacedirections(img) === ((1,0),(0,1))
    img.spacedirections = "broken"
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
                    spatialproperties=Set([:vector,:matrix,:tuple]),
                    vector=[1,2],
                    matrix=[1 3; 2 4],
                    tuple=(1,2))
    for imgp in (img', permutedims(img, (2,1)), PermutedDimsArray(img, (2,1)))
        @test arraydata(imgp) == arraydata(img)'
        @test imgp.vector == [2,1]
        @test imgp.matrix == [4 2; 3 1]
        @test imgp.tuple  == (2,1)
        @test img.vector == [1,2]
        @test img.matrix == [1 3; 2 4]
        @test img.tuple  == (1,2)
    end
end

@testset "dimchange" begin
    M = ImageMeta([1,2,3,4])
    Mp = M'
    @test ndims(Mp) == 2 && isa(Mp, ImageMeta)

    M = ImageMeta([1,2,3,4],
                  spatialproperties=[:vector],
                  vector=[1])
    @test_throws ErrorException M'
end

@testset "Inference" begin
    # Heavily-nested types tax Julia's inference capabilities
    a = rand(N2f14, 3, 5, 2)
    a1 = a[1,1,1]
    aa = AxisArray(a, Axis{:y}(0:2), Axis{:x}(0:4), Axis{:z}(0:1))
    cv = colorview(RGB, aa, zeroarray, aa)
    @test @inferred(cv[1,1,1]) == RGB(a1, zero(a1), a1)
    au = AxisArray(a, Axis{:y}(0m:1m:2m), Axis{:x}(0m:1m:4m), Axis{:z}(0m:1m:1m))
    cv = colorview(RGB, au, zeroarray, au)
    @test @inferred(cv[1,1,1]) == RGB(a1, zero(a1), a1)
    am = ImageMeta(au)
    cv = colorview(RGB, am, zeroarray, am)
    @inferred cv[1,1,1]
end

@testset "AxisArray_CartesianIndex" begin
    M = reshape([1,2,3,4], 2,2)
    M = AxisArray(M)
    img = ImageMeta(M)
    @test 1 == img[1, CartesianIndex(1)] #Int and cartesianindex
    @test maximum(img,dims=2) == reshape([3,4],2,1)
    @test mean(img,dims=1) == reshape([1.5,3.5], 1,2)
end

@testset "AbstractDicts" begin
    img = ImageMeta(AxisArray(rand(3,5,8),
                              Axis{:x}(1:3),
                              Axis{:y}(1:5),
                              Axis{:time}(0.1:0.1:0.8)),
                    IdDict{Symbol,Any}(:pixelspacing => (1,1)))
    @test @inferred(pixelspacing(img)) === (1,1)

    rand_img = rand(10,10)
    dict = Dict(:attr1 => 3, :attr2 => 4)
    @test ImageMeta(rand_img, dict)[1, 1] == rand_img[1, 1]
end

@testset "OffsetArrays" begin
    if isdefined(OffsetArrays, :centered)
        # OffsetArrays v1.9+ only
        @testset "centered" begin
            check_range(r, f, l) = (@test first(r) == f; @test last(r) == l)
            check_range_axes(r, f, l) = check_range(axes(r)[1], f, l)

            check_range(Base.axes(OffsetArrays.centered(1:3))[1], -1, 1)
            a = AxisArray(rand(3, 3), Axis{:y}(0.1:0.1:0.3), Axis{:x}(1:3))
            am = ImageMeta(a; prop1="simple")
            ca = OffsetArrays.centered(am)
            axs = Base.axes(ca)
            check_range(axs[1], -1, 1)
            check_range(axs[2], -1, 1)
            @test Base.axes(arraydata(ca)) == axs
        end
    end
end

@testset "restrict" begin
    imgcol = rand(RGB{Float32}, 5, 6)

    for dims in [(), (1, ), (1, 2)]
        imgmeta = ImageMeta(imgcol, myprop=1)

        out = restrict(imgmeta, dims)
        @test out isa ImageMeta
        @test out.myprop == 1
        @test arraydata(out) == restrict(arraydata(imgmeta), dims)
    end

    imgcolax = AxisArray(imgcol, :y, :x)
    for ax in [Axis{:x}, Axis{:y}]
        imgmetaax = ImageMeta(imgcolax, myprop=1)
        out = restrict(imgmetaax, ax)
        @test out isa ImageMeta
        @test out.myprop == 1
        @test arraydata(out) == restrict(arraydata(imgmetaax), ax)
    end
end

nothing

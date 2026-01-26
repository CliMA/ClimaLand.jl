@testset "restrict" begin
    imgcol = rand(RGB{Float64}, 5, 6)
    imgcolax = AxisArray(imgcol, :y, :x)
    imgr = restrict(imgcolax, (1,2))
    @test pixelspacing(imgr) == (2,2)
    @test pixelspacing(imgcolax) == (1,1)  # issue https://github.com/JuliaImages/Images.jl/issues/347

    imgr = @inferred(restrict(imgcolax, Axis{:y}))
    @test pixelspacing(imgr) == (2, 1)
    imgr = @inferred(restrict(imgcolax, Axis{:x}))
    @test pixelspacing(imgr) == (1, 2)

    @testset "OffsetArray" begin
        A = ones(4, 5)
        Ao = OffsetArray(A, 0:3, -2:2)
        Aor = restrict(Ao)
        @test axes(Aor) == axes(OffsetArray(restrict(A), 0:2, -1:1))

        # restricting AxisArrays with offset axes
        AA = AxisArray(Ao, Axis{:y}(0:3), Axis{:x}(-2:2))
        @test axes(AA) === axes(Ao)
        AAr = restrict(AA)
        axs = axisvalues(AAr)
        @test axes(axs[1])[1] == axes(Aor)[1]
        @test axes(axs[2])[1] == axes(Aor)[2]
        AAA = AxisArray(Aor, axs)  # just test that it's constructable (another way of enforcing agreement)
    end
end

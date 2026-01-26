using StrideArraysCore, ThreadingUtilities, Aqua
# using InteractiveUtils
using Test

@inline relu(x) = (y = zero(x); ifelse(x < y, y, x))

function closeopensum(x)
  s = zero(eltype(x))
  @inbounds @simd for i ∈ StrideArraysCore.CloseOpen(length(x))
    s += x[i+1]
  end
  s
end
function closeopensumfastmath(x)
  s = zero(eltype(x))
  @inbounds @fastmath for i ∈ StrideArraysCore.CloseOpen(length(x))
    s += x[i+1]
  end
  s
end
function cartesianindexsum(A)
  s = zero(eltype(A))
  @inbounds @simd for I ∈ CartesianIndices(A)
    s += A[I]
  end
  s
end
allocated_cartesianindexsum(x) = @allocated cartesianindexsum(x)
@noinline function conv!(f::F, _C, _A, _K, _b) where {F}
  C = StrideArraysCore.LayoutPointers.zero_offsets(_C)
  A = StrideArraysCore.LayoutPointers.zero_offsets(_A)
  K = StrideArraysCore.LayoutPointers.zero_offsets(_K)
  b = StrideArraysCore.LayoutPointers.zero_offsets(_b)

  for j₁ ∈ axes(C, 1), j₂ ∈ axes(C, 2), o ∈ axes(K, 4)
    s = zero(eltype(C))
    for k₁ ∈ axes(K, 1), k₂ ∈ axes(K, 2), i ∈ axes(K, 3)
      s += A[j₁+k₁, j₂+k₂, i] * K[k₁, k₂, i, o]
    end
    C[j₁, j₂, o] = f(s + b[o])
  end

  return C
end

function conv_group!(f::F, C, A, K, b, range) where {F}
  conv!(
    f,
    view(C, :, :, range),
    view(A, :, :, range),
    view(K, :, :, :, range),
    view(b, range),
  )
  return C
end
randn_stridearray(size...) = StrideArray(randn(Float32, size...), static.(size))

@noinline function conv_grouped!(f::F, C, A, K, b) where {F}
  conv_group!(f, C, A, K, b, static(1 + 0 * 4):static(1 * 4))
  conv_group!(f, C, A, K, b, static(1 + 1 * 4):static(2 * 4))
  conv_group!(f, C, A, K, b, static(1 + 2 * 4):static(3 * 4))
  conv_group!(f, C, A, K, b, static(1 + 3 * 4):static(4 * 4))
  conv_group!(f, C, A, K, b, static(1 + 4 * 4):static(5 * 4))
  conv_group!(f, C, A, K, b, static(1 + 5 * 4):static(6 * 4))
  conv_group!(f, C, A, K, b, static(1 + 6 * 4):static(7 * 4))
  conv_group!(f, C, A, K, b, static(1 + 7 * 4):static(8 * 4))
  return C
end

@inline function unpad(A::AbstractArray{<:Any,3})

  return view(A, static(2):static(9), static(2):static(9), :)
end

@testset "StrideArraysCore.jl" begin
  # Currently StrideArraysCore commits piracy with zero_offsets(A::AbstractArray) and preserve_buffer(A::MemoryBuffer)
  Aqua.test_all(StrideArraysCore; piracies = false)

  @testset "StrideArrays Basic" begin

    Acomplex = StrideArray{Complex{Float64}}(undef, (StaticInt(4), StaticInt(5)))
    @test @inferred(StrideArraysCore.ArrayInterface.known_size(Acomplex)) === (4, 5)
    Acomplex .= rand.(Complex{Float64})
    @test StrideArray(Acomplex) ===
          reinterpret(reshape, Complex{Float64}, reinterpret(reshape, Float64, Acomplex))
    @test StrideArray(Acomplex) ===
          reinterpret(Complex{Float64}, reinterpret(Float64, Acomplex))
    @test StrideArray(Acomplex) === reinterpret(
      reshape,
      Complex{Float64},
      reinterpret(reshape, Complex{UInt64}, Acomplex),
    )
    @test StrideArraysCore.static_size(Acomplex) ===
          StrideArraysCore.static_size(StrideArray(Acomplex)) ===
          (StaticInt(4), StaticInt(5))
    @test StrideArraysCore.static_size(reinterpret(reshape, Float64, Acomplex)) ===
          (StaticInt(2), StaticInt(4), StaticInt(5))
    A = rand(100, 100)
    B = copy(A)
    C = StrideArraysCore.PtrArray(A)
    @test typeof(@inferred(StrideArraysCore.axes(Acomplex))) ==
          @inferred(StrideArraysCore.ArrayInterface.axes_types(typeof(Acomplex)))
    @test typeof(@inferred(StrideArraysCore.axes(C))) ==
          @inferred(StrideArraysCore.ArrayInterface.axes_types(typeof(C)))
    @test typeof(
      @inferred(
        StrideArraysCore.axes(StrideArraysCore.LayoutPointers.zero_offsets(Acomplex))
      )
    ) == @inferred(
      StrideArraysCore.ArrayInterface.axes_types(
        typeof(StrideArraysCore.LayoutPointers.zero_offsets(Acomplex)),
      )
    )

    @test similar(C) isa StrideArraysCore.StrideArray
    @test view(C, :) === vec(C)
    let D = similar(C, Float32)
      @test D isa StrideArraysCore.StrideArray
      @test eltype(D) === Float32
    end
    GC.@preserve A begin
      #TODO: eliminate need to use `@inbounds` for simd
      # @test closeopensum(C) == closeopensum(A)
      # @code_llvm closeopensum(C)

      # @code_llvm closeopensum(A)
      @test closeopensumfastmath(C) == closeopensumfastmath(A)
      @test foldl(+, C) == foldl(+, A)
      # @test sum(C) == sum(A)
      @test closeopensum(C) ≈ closeopensumfastmath(C) ≈ sum(C)
      @test A == B
      C .*= 3
      @test A == 3 .* B
      @test C == A
      @test C == 3 .* B

      D = copy(A)
      Cslice = view(C, 23:48, 17:89)
      @test Base.stride(Cslice, 1) == 1
      Base.stride(C, 1) ==
      StrideArraysCore.static_strides(Cslice, 1) ==
      StrideArraysCore.static_strides(Cslice, static(1)) ==
      StrideArraysCore.static_strides(C, 1) ==
      StrideArraysCore.static_strides(C, static(1))
      @test Base.stride(Cslice, 2) ==
            Base.stride(C, 2) ==
            StrideArraysCore.static_strides(Cslice, 2) ==
            StrideArraysCore.static_strides(Cslice, static(2)) ==
            StrideArraysCore.static_strides(C, 2) ==
            StrideArraysCore.static_strides(C, static(2))
      if VERSION >= v"1.9.0-DEV.569"
        @test Base.stride(C, 3) ==
              StrideArraysCore.static_strides(C, 3) ==
              StrideArraysCore.static_strides(C, static(3))

      end
      @test Base.stride(C, 3) == StrideArraysCore.static_strides(C, 3)
      @test Base.stride(Cslice, 3) == StrideArraysCore.static_strides(Cslice, 3)
      @test Base.stride(Cslice, 3) == StrideArraysCore.static_strides(Cslice, static(3))
      @test Base.strides(Cslice) ==
            Base.strides(C) ==
            StrideArraysCore.static_strides(Cslice) ==
            StrideArraysCore.static_strides(C)

      Cslice .= 2
      @test D != C
      D[23:48, 17:89] .= 2
      @test D == C
      @test C === view(C, :, :)
      @test @inferred(size(view(C, StaticInt(1):StaticInt(8), :), 1)) === 8
      @test @inferred(
        StrideArraysCore.static_size(view(C, StaticInt(1):StaticInt(8), :), StaticInt(1))
      ) === StaticInt(8)
      @test @inferred(
        (static ∘ size)(view(C, StaticInt(1):StaticInt(8), :), StaticInt(1))
      ) === StaticInt(8)
      @test @inferred(
        StrideArraysCore.static_size(view(C, StaticInt(1):StaticInt(8), :), 1)
      ) === 8
      @test @inferred(
        StrideArraysCore.static_size(view(C, StaticInt(1):StaticInt(8), :), StaticInt(1))
      ) === StaticInt(8)
      @test @inferred(
        StrideArraysCore.ArrayInterface.known_size(
          typeof(view(C, StaticInt(1):StaticInt(8), :)),
        )
      ) === (8, nothing)

      @test C isa PtrArray
      @test C' isa PtrArray
      @test @inferred(permutedims(C, Val((2, 1)))) isa PtrArray
      @test @inferred(adjoint(C)) == D'
      @test @inferred(axes(C)) == axes(D)
      @test @inferred(axes(C')) == axes(D')
      @test @inferred(eachindex(view(C, :, 2:6))) == 1:(5*size(C, 1))
      @test @inferred(eachindex(view(C', 2:6, :)')) == 1:(5*size(C, 1))
      @test @inferred(eachindex(view(C, 2:6, :))) == CartesianIndices((5, size(C, 2)))
      @test @inferred(
        axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), StaticInt(1))
      ) === StrideArraysCore.CloseOpen(StaticInt(5))
      @test @inferred(
        axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), StaticInt(1))
      ) === StrideArraysCore.CloseOpen(StaticInt(5))
      @test @inferred(
        length(axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), StaticInt(1)))
      ) === 5
      if VERSION >= v"1.6"
        @test @inferred(
          StrideArraysCore.static_length(
            axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), StaticInt(1)),
          )
        ) === StaticInt(5)
      else
        @test @inferred(
          StrideArraysCore.static_length(
            axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), StaticInt(1)),
          )
        ) === 5
      end
      @test @inferred(
        StrideArraysCore.static_length(
          axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), StaticInt(1)),
        )
      ) === StaticInt(5)
      @test @inferred(
        (static ∘ length)(
          axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), StaticInt(1)),
        )
      ) === StaticInt(5)
      @test @inferred(
        length(axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :), 1))
      ) === 5
      ax1, ax2 = axes(StrideArraysCore.zview(C, StaticInt(2):StaticInt(6), :))
      @test StrideArraysCore.ArrayInterface.known_length(ax1) == 5
      @test StrideArraysCore.ArrayInterface.known_length(ax2) ===
            StrideArraysCore.ArrayInterface.known_length(1:1)
      @test StrideArraysCore.ArrayInterface.known_first(ax1) == 0
      @test StrideArraysCore.ArrayInterface.known_first(ax2) == 0
      @test StrideArraysCore.ArrayInterface.known_step(ax1) == 1
      @test StrideArraysCore.ArrayInterface.known_step(ax2) == 1
      @test StrideArraysCore.ArrayInterface.known_last(ax1) == 4
      @test StrideArraysCore.ArrayInterface.known_last(ax2) ===
            StrideArraysCore.ArrayInterface.known_length(1:1)
    end
    W = rand(2, 3, 4)
    X = PtrArray(W)
    @test @inferred(StrideArraysCore.ArrayInterface.known_size(X)) ===
          (nothing, nothing, nothing)
    GC.@preserve W begin
      @test W == X
      @test permutedims(W, (1, 2, 3)) == @inferred(permutedims(X, Val((1, 2, 3))))
      @test permutedims(W, (1, 3, 2)) == @inferred(permutedims(X, Val((1, 3, 2))))
      @test permutedims(W, (2, 1, 3)) == @inferred(permutedims(X, Val((2, 1, 3))))
      @test permutedims(W, (2, 3, 1)) == @inferred(permutedims(X, Val((2, 3, 1))))
      @test permutedims(W, (3, 1, 2)) == @inferred(permutedims(X, Val((3, 1, 2))))
      @test permutedims(W, (3, 2, 1)) == @inferred(permutedims(X, Val((3, 2, 1))))
      if Base.JLOptions().check_bounds == 1
        @test_throws BoundsError X[length(X)+1]
        @test_throws BoundsError X[-4]
        @test_throws BoundsError X[2, 5, 3]
      end
      @test cartesianindexsum(W) ≈ cartesianindexsum(X)
      @test iszero(allocated_cartesianindexsum(X))
    end
    @test X === PtrArray(pointer(X), size(X))
    y = rand(77)
    GC.@preserve y begin
      z = PtrArray(y)
      @test y == z
      @test y' == z'
      @test PtrArray(y') === z'
      @test PtrArray(transpose(y)) === transpose(PtrArray(y))
      @test pointer(y) === pointer(z)
      if Base.JLOptions().check_bounds == 1
        @test_throws BoundsError z[-8]
        @test_throws BoundsError z[88]
      end
    end
  end
  @testset "ThreadingUtilities" begin
    xu = zeros(UInt, 100)
    x = rand(100)
    y = rand(100)
    z = rand(100)
    t = (x, y, z)
    pt, gt = StrideArraysCore.object_and_preserve(t)
    greet = ["Hello", "world", "hang", "in", "there"]
    pg, gg = StrideArraysCore.object_and_preserve(greet)
    ph, gh = StrideArraysCore.object_and_preserve(greet[1])
    pht, ght = StrideArraysCore.object_and_preserve((3, greet[1]))
    GC.@preserve xu gt gg gh ght begin
      ThreadingUtilities.store!(pointer(xu), pt, 0)
      @test @inferred(ThreadingUtilities.load(pointer(xu), typeof(pt), 0)) ==
            (6sizeof(UInt), t)
      offset = sizeof(UInt)
      for a ∈ t
        _p, g = StrideArraysCore.object_and_preserve(a)
        ThreadingUtilities.store!(pointer(xu), _p, offset)
        offset, p = @inferred(ThreadingUtilities.load(pointer(xu), typeof(_p), offset))
        @test g === a
        @test p isa PtrArray
        @test pointer(p) === pointer(a)
        @test p !== a
        @test p == a
      end
      new_offset = ThreadingUtilities.store!(pointer(xu), pg, offset)
      @test @inferred(ThreadingUtilities.load(pointer(xu), typeof(pg), offset)) ==
            (new_offset, greet)
      ThreadingUtilities.store!(pointer(xu), ph, offset)
      @test @inferred(ThreadingUtilities.load(pointer(xu), typeof(ph), offset)) ===
            (offset + sizeof(UInt), greet[1])
      ThreadingUtilities.store!(pointer(xu), pht, offset)
      @test @inferred(ThreadingUtilities.load(pointer(xu), typeof(pht), offset)) ===
            (offset + 2sizeof(UInt), (3, greet[1]))
    end

    greetsa = StrideArray(greet)
    @test greetsa == greet
    greetsa[1] = "howdy"
    greetsa[5] = "there!"
    @test greet[1] == "howdy"
    @test greet[5] == "there!"
    A = rand(ComplexF64, 4, 5)
    Asa = StrideArray(A)
    @test A == Asa
    Asa[3, 4] = 1234.5 + 678.9im
    @test A[3, 4] == 1234.5 + 678.9im
  end
  @testset "StrideArrays Initialization" begin
    A = StrideArray{Float64}(undef, (3, 5))
    @test StrideArraysCore.static_size(A) === (3, 5)
    @test StrideArraysCore.static_strides(A) === (StaticInt(1), 3)
    @test StrideArraysCore.static_length(A) === 15

    B = StrideArray(undef, (StaticInt(3), 5))
    @test StrideArraysCore.static_strides(B) === (StaticInt(1), StaticInt(3))
    @test StrideArraysCore.static_size(B) === (StaticInt(3), 5)
    @test StrideArraysCore.static_length(B) === 15

    D = StrideArray(undef, (StaticInt(3), StaticInt(5)))
    @test StrideArraysCore.static_strides(D) === (StaticInt(1), StaticInt(3))
    @test StrideArraysCore.static_size(D) === (StaticInt(3), StaticInt(5))
    @test StrideArraysCore.static_length(D) === StaticInt(15)
    @test D isa StrideArraysCore.StaticStrideArray
    for C ∈ Any[A, B, D]
      @test strides(C) === (1, 3)
      @test size(C) === (3, 5)
      @test StrideArraysCore.offsets(C) === (StaticInt(1), StaticInt(1))
      @test StrideArraysCore.offsets(StrideArraysCore.zero_offsets(C)) ===
            (StaticInt(0), StaticInt(0))
      C[2, 3] = 4
      StrideArraysCore.zero_offsets(C)[2, 3] = -10
      @test C[2, 3] === StrideArraysCore.zero_offsets(C)[1, 2] === 4.0
      @test C[3, 4] === StrideArraysCore.zero_offsets(C)[2, 3] === -10.0
    end
    @test all(iszero, StrideArray(zero, static(4), static(8)))
    @test all(iszero, StrideArray(zero, 1000, 2000))
    @test all(isone, StrideArray(one, static(4), static(8)))
    @test all(isone, StrideArray(one, 100, 200))
  end
  @testset "views" begin
    B0 = reshape(collect(1:12), 3, 4)
    B1 = StrideArray(B0)
    @test view(B0, :, 4:-1:1) == view(B1, :, 4:-1:1) == B1[:, 4:-1:1]
    @test view(B0, :, 1:2:4) == view(B1, :, 1:2:4) == B1[:, 1:2:4]
    @test view(B1, :, 4:-1:1) === B1[:, 4:-1:1]
    @test view(B1, 2, 1:2:4) === B1[2, 1:2:4]
    @test view(B0, 2, 4:-1:1) == view(B1, 2, 4:-1:1) == B1[2, 4:-1:1]
    @test view(B0, 2, 1:2:4) == view(B1, 2, 1:2:4) == B1[2, 1:2:4]
    @test view(B1, 2, 4:-1:1) === B1[2, 4:-1:1]
    @test view(B1, 2, 1:2:4) === B1[2, 1:2:4]
    A = StrideArray{Float64}(undef, (100, 100)) .= rand.()
    vA = view(A, 3:40, 2:50)
    @test vA === A[3:40, 2:50]
    vAslice = view(A, :, 2:50)
    vaslice = view(A, 2:50, 4)
    x = StrideArray{Float64}(undef, (100,)) .= rand.()
    vxslice = view(x, 2:50)
    @test vxslice === x[2:50]
    for (i, n) in enumerate(2:50)
      @test vaslice[i] == A[n, 4]
      @test vxslice[i] == x[n]
      for (j, m) in enumerate(3:40)
        @test A[m, n] == vA[j, i]
      end
      for j = 1:50
        @test A[j, n] == vAslice[j, i]
      end
    end
    a = StrideArray(zeros(77))
    b = a[10:20]
    c = @view(a[20:30])
    b[4] = 7
    @test a[13] == 7
    c[4] = 37
    @test a[23] == 37
    @test typeof(a) === typeof(b) === typeof(c)
  end
  @testset "BitPtrArray" begin
    b = collect(1:10) .> 5
    @test sprint((io, t) -> show(io, t), StrideArray(b)) == """
              Bool[0, 0, 0, 0, 0, 1, 1, 1, 1, 1]"""
    @test sprint((io, t) -> show(io, t), StrideArray(b)') == """
              Bool[0 0 0 0 0 1 1 1 1 1]"""
  end
  @testset "PtrArray0" begin
    x = collect(0:3)
    pzx = StrideArraysCore.PtrArray0(pointer(x), (4,))
    GC.@preserve x begin
      for i = 0:3
        @test pzx[i] == pzx[i, 1] == i
      end
    end
  end

  @testset "conv" begin


    A = randn_stridearray(10, 10, 32)
    C = randn_stridearray(10, 10, 32)
    K = randn_stridearray(3, 3, 4, 32)
    b = randn_stridearray(32)
    Ca = Array(C)
    @test Ca == C
    @test unpad(Ca) == unpad(C)
    conv_grouped!(relu, unpad(C), A, K, b)
    #Ca = Array(C);
    Aa = Array(A)
    Ka = Array(K)
    ba = Array(b)
    @test Aa == A
    @test Ka == K
    @test ba == b
    for i = 1:3
      @test axes(A, i) == axes(Aa, i)
      @test axes(K, i) == axes(Ka, i)
      @test axes(C, i) == axes(Ca, i)
    end
    @test axes(C, 4) == axes(Ca, 4)
    @test axes(b, 1) == axes(ba, 1)
    conv_grouped!(relu, unpad(Ca), Aa, Ka, ba)
    @test C ≈ Ca

  end
  @testset "square" begin
    A = rand(5, 5)
    GC.@preserve A begin
      B = StrideArraysCore.SquarePtrMatrix(pointer(A), 5)
      @test A == B
    end
  end

  @testset "PtrArray slicing" begin
    # Testing for the bug in https://github.com/JuliaSIMD/StrideArrays.jl/issues/88
    x = [1 3; 2 4]
    GC.@preserve x begin
      ptrarrays = (
        StrideArraysCore.PtrArray(pointer(x), (2, 2)),
        StrideArraysCore.PtrArray(pointer(x), (StaticInt(2), StaticInt(2))),
      )
      for y in ptrarrays
        # This triggered the bug, a 2D index into a Vector, which is used by display/show
        y1 = y[1, :]
        @test y1[2, 1] == x[1, 2]

        # Some general tests with slicing and bounds checking
        @test y[1, 1] == x[1, 1]
        @test y[2, 1] == x[2, 1]
        @test y[1, 2] == x[1, 2]
        @test y[2, 2] == x[2, 2]
        @test y1[1] == x[1, 1]
        @test y1[2] == x[1, 2]
        y2 = y[2, :]
        @test y2[1] == x[2, 1]
        @test y2[2] == x[2, 2]
        y3 = y[:, 1]
        @test y3[1] == x[1, 1]
        @test y3[2] == x[2, 1]
        @test_throws BoundsError y[1:5]
        @test_throws BoundsError y[3, :]
        @test_throws BoundsError y[:, 3]
        @test_throws BoundsError y[2:3, :]
        @test_throws BoundsError y1[3]
        @test_throws BoundsError y1[1:3]
        @test_throws BoundsError y2[3]
        @test_throws BoundsError y2[1:3]
      end
    end
  end
end

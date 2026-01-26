# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function matrix_vector_division_test(x)
    if size(x) == (3, 3)
        A = [
            2 1 0
            1 2 1
            0 1 2
        ]
        @test_rewrite A * x / 2
    end
end

function _xAx_test(x::AbstractVector, A::AbstractMatrix)
    for t in [transpose, adjoint]
        @test_rewrite(t(x) * A)
        # Complex expression
        @test_rewrite(t(x) * ones(Int, size(A)...))
        @test_rewrite(t(x) * A * x)
        # Complex expression
        @test_rewrite(t(x) * ones(Int, size(A)...) * x)
        @test_rewrite reshape(x, (1, length(x))) * A * x .- 1
        @test_rewrite t(x) * A * x .- 1
        @test_rewrite t(x) * A * x - 1
        @test_rewrite t(x) * x + t(x) * A * x
        @test_rewrite t(x) * x - t(x) * A * x
        @test MA.promote_operation(*, typeof(t(x)), typeof(A), typeof(x)) ==
              typeof(t(x) * A * x)
        @test MA.promote_operation(*, typeof(t(x)), typeof(x)) ==
              typeof(t(x) * x)
        @test_rewrite t(x) * x + 2 * t(x) * A * x
        @test_rewrite t(x) * x - 2 * t(x) * A * x
        @test_rewrite t(x) * A * x + 2 * t(x) * x
        @test_rewrite t(x) * A * x - 2 * t(x) * x
        @test MA.promote_operation(
            *,
            Int,
            typeof(t(x)),
            typeof(A),
            typeof(x),
        ) == typeof(2 * t(x) * A * x)
        @test MA.promote_operation(*, Int, typeof(t(x)), typeof(x)) ==
              typeof(2 * t(x) * x)
    end
end
function _xABx_test(x::AbstractVector, A::AbstractMatrix, B::AbstractMatrix)
    for t in [transpose, adjoint]
        @test_rewrite t(t(x) * A) + 2B * x
        @test_rewrite t(t(x) * A) + 2B * x .- 1
        @test_rewrite t(t(x) * A) + 2B * x .- [length(x):-1:1;]
        @test_rewrite t(t(x) * A) + 2B * x - [length(x):-1:1;]
    end
end

function _matrix_vector_test(x::AbstractVector, A::AbstractMatrix)
    _xAx_test(x, A)
    B = sparse(A)
    _xAx_test(x, B)

    @test MA.isequal_canonical(A * x, B * x)
    @test MA.isequal_canonical(A * x, MA.@rewrite(B * x))
    @test MA.isequal_canonical(MA.@rewrite(A * x), MA.@rewrite(B * x))
    @test MA.isequal_canonical(x' * A, x' * B)
    @test MA.isequal_canonical(x' * A, MA.@rewrite(x' * B))
    @test MA.isequal_canonical(MA.@rewrite(x' * A), MA.@rewrite(x' * B))
    @test MA.isequal_canonical(x'A * x, x' * B * x)
    @test MA.isequal_canonical(x' * A * x, MA.@rewrite(x' * B * x))
    @test MA.isequal_canonical(MA.@rewrite(x' * A * x), MA.@rewrite(x' * B * x))

    _xABx_test(x, A, A)
    _xABx_test(x, A, B)
    _xABx_test(x, B, A)
    return _xABx_test(x, B, B)
end

function matrix_vector_test(x)
    add_test(x, x)
    if size(x) != (3,)
        return
    end
    A = [
        2 1 0
        1 2 1
        0 1 2
    ]

    @test_rewrite x .+ A * x
    @test_rewrite A * x .+ A * x
    @test_rewrite x .- A * x
    @test_rewrite A * x .- A * x
    @test_rewrite A .+ (A + A)^2
    @test_rewrite A .- (A + A)^2
    @test_rewrite A * x .+ (A + A)^2 * x
    @test_rewrite A * x .- (A + A)^2 * x
    @test_rewrite A * x .- (A + A)^2 * x

    @test MA.isequal_canonical(-x, [-x[1], -x[2], -x[3]])
    xAx =
        2x[1] * x[1] + 2x[1] * x[2] + 2x[2] * x[2] + 2x[2] * x[3] + 2x[3] * x[3]
    @test MA.isequal_canonical(x' * A * x, xAx)
    y = A * x
    @test MA.isequal_canonical(x' * y, xAx)
    @test MA.isequal_canonical(
        y,
        [
            2x[1] + x[2]
            2x[2] + x[1] + x[3]
            x[2] + 2x[3]
        ],
    )
    @test MA.isequal_canonical(
        -y,
        [
            -2x[1] - x[2]
            -x[1] - 2x[2] - x[3]
            -x[2] - 2x[3]
        ],
    )
    @test MA.isequal_canonical(
        x' * A,
        [
            2x[1] + x[2]
            x[1] + 2x[2] + x[3]
            x[2] + 2x[3]
        ]',
    )

    @test MA.isequal_canonical(
        y .+ 1,
        [
            2x[1] + x[2] + 1
            x[1] + 2x[2] + x[3] + 1
            x[2] + 2x[3] + 1
        ],
    )
    @test MA.isequal_canonical(
        y .- 1,
        [
            2x[1] + x[2] - 1
            x[1] + 2x[2] + x[3] - 1
            x[2] + 2x[3] - 1
        ],
    )
    @test MA.isequal_canonical(
        y .+ 2ones(Int, 3),
        [
            2x[1] + x[2] + 2
            x[1] + 2x[2] + x[3] + 2
            x[2] + 2x[3] + 2
        ],
    )
    @test MA.isequal_canonical(
        y .- 2ones(Int, 3),
        [
            2x[1] + x[2] - 2
            x[1] + 2x[2] + x[3] - 2
            x[2] + 2x[3] - 2
        ],
    )
    @test MA.isequal_canonical(
        2ones(Int, 3) .+ y,
        [
            2x[1] + x[2] + 2
            x[1] + 2x[2] + x[3] + 2
            x[2] + 2x[3] + 2
        ],
    )
    @test MA.isequal_canonical(
        2ones(Int, 3) .- y,
        [
            -2x[1] - x[2] + 2
            -x[1] - 2x[2] - x[3] + 2
            -x[2] - 2x[3] + 2
        ],
    )
    @test MA.isequal_canonical(
        y .+ x,
        [
            3x[1] + x[2]
            x[1] + 3x[2] + x[3]
            x[2] + 3x[3]
        ],
    )
    @test MA.isequal_canonical(
        x .+ y,
        [
            3x[1] + x[2]
            x[1] + 3x[2] + x[3]
            x[2] + 3x[3]
        ],
    )
    @test MA.isequal_canonical(
        2y .+ 2x,
        [
            6x[1] + 2x[2]
            2x[1] + 6x[2] + 2x[3]
            2x[2] + 6x[3]
        ],
    )
    @test MA.isequal_canonical(
        y .- x,
        [
            x[1] + x[2]
            x[1] + x[2] + x[3]
            x[2] + x[3]
        ],
    )
    @test MA.isequal_canonical(
        x .- y,
        [
            -x[1] - x[2]
            -x[1] - x[2] - x[3]
            -x[2] - x[3]
        ],
    )
    @test MA.isequal_canonical(
        y .+ x[:],
        [
            3x[1] + x[2]
            x[1] + 3x[2] + x[3]
            x[2] + 3x[3]
        ],
    )
    @test MA.isequal_canonical(
        x[:] .+ y,
        [
            3x[1] + x[2]
            x[1] + 3x[2] + x[3]
            x[2] + 3x[3]
        ],
    )

    _matrix_vector_test(x, A)

    A = [
        1 2 3
        0 4 5
        6 0 7
    ]
    return _matrix_vector_test(x, A)
end

_constant(x) = reshape(collect(1:length(x)), size(x)...)

function non_array_test(x, x2)
    # This is needed to compare arrays that have nonstandard indexing
    function elements_equal(
        A::AbstractArray{T,N},
        B::AbstractArray{T,N},
    ) where {T,N}
        return all(MA.isequal_canonical(a, b) for (a, b) in zip(A, B))
    end

    @test elements_equal(+x, +x2)
    @test elements_equal(-x, -x2)
    @test elements_equal(x .+ first(x), x2 .+ first(x2))
    @test elements_equal(x .- first(x), x2 .- first(x2))
    @test elements_equal(first(x) .- x, first(x2) .- x2)
    @test elements_equal(first(x) .+ x, first(x2) .+ x2)
    @test elements_equal(2 .* x, 2 .* x2)
    @test elements_equal(first(x) .+ x2, first(x2) .+ x)
    @test MA.isequal_canonical(sum(x), sum(x2))
    if !MA._one_indexed(x2)
        @test_throws DimensionMismatch x + x2
    end
    if x isa AbstractVector
        # `diagm` only define with `Pair` in Julia v1.0 and v1.1
        @testset "diagm" begin
            if !MA._one_indexed(x2) && eltype(x2) <: MA.AbstractMutable
                @test_throws AssertionError diagm(x2)
            else
                @test MA.isequal_canonical(diagm(0 => x), diagm(0 => x2))
                @test MA.isequal_canonical(diagm(x), diagm(x2))
            end
        end
    end
end

function non_array_test(x)
    non_array_test(x, x)
    if x isa AbstractVector
        non_array_test(x, view(x, :))
    end
    if x isa AbstractVector || x isa AbstractMatrix
        non_array_test(x, sparse(x))
    end
end

function dot_test(x)
    @test_rewrite dot(first(x), first(x))
    @test_rewrite dot(2, first(x))
    @test_rewrite dot(first(x), 2)

    A = _constant(x)
    @test_rewrite dot(A, x)
    @test_rewrite dot(x, A)

    y = repeat(x, outer = (one.(size(x))..., size(x, 1)))
    @test_rewrite dot(x, _constant(x)) - dot(y, _constant(y))
end

function sum_test(matrix)
    @test_rewrite sum(matrix)
    if matrix isa AbstractMatrix
        @test_rewrite sum([
            2matrix[i, j] for i in axes(matrix, 1), j in axes(matrix, 2)
        ])
        @test_rewrite sum(
            2matrix[i, j] for i in axes(matrix, 1), j in axes(matrix, 2)
        )
    end
end

function sum_multiplication_test(matrix)
    if matrix isa AbstractMatrix
        @test_rewrite sum([
            2matrix[i, j]^2 for i in axes(matrix, 1), j in axes(matrix, 2)
        ])
        @test_rewrite sum(
            2matrix[i, j]^2 for i in axes(matrix, 1), j in axes(matrix, 2)
        )
    end
end

function _is_supported(op::Function, Ts::Type...)
    try
        op(zero.(Ts)...)
        true
    catch
        false
    end
end

function transpose_test(x)
    if x isa Vector
        y = reshape(x, 1, length(x))
        @test MA.isequal_canonical(x', y)
        @test MA.isequal_canonical(copy(transpose(x)), y)
    end
    # If the element type does not support multiplication, e.g.
    # JuMP or MOI quadratic functions, then we should skip these tests.
    if x isa AbstractMatrix && _is_supported(*, eltype(x), eltype(x))
        y = [x[i, j] for j in axes(x, 2), i in axes(x, 1)]
        @test MA.isequal_canonical(x', y)
        @test MA.isequal_canonical(copy(transpose(x)), y)
        if size(x, 1) == size(x, 2)
            yy = y * y
            @test MA.isequal_canonical(x' * x', yy)
            @test MA.isequal_canonical(transpose(x) * transpose(x), yy)
            z = _constant(x)
            zt = copy(z')
            @test MA.isequal_canonical(adjoint(x) * adjoint(z), y * zt)
            @test MA.isequal_canonical(transpose(x) * transpose(z), y * zt)
            @test MA.isequal_canonical(adjoint(z) * adjoint(x), zt * y)
            @test MA.isequal_canonical(transpose(z) * transpose(x), zt * y)
        end
    end
    if x isa AbstractVector || x isa AbstractMatrix
        @test_rewrite(x')
        @test_rewrite(transpose(x))
        @test (x')' == x
        @test transpose(transpose(x)) == x
    end
end

_matrix(x::Matrix) = x
_matrix(x::AbstractMatrix) = Matrix(x)
function _matrix_equal(x::AbstractMatrix, y::AbstractMatrix)
    return MA.isequal_canonical(_matrix(x), _matrix(y))
end

function _broadcast_test(x, A)
    B = sparse(A)
    y = SparseMatrixCSC(
        size(x)...,
        copy(B.colptr),
        copy(B.rowval),
        collect(vec(x)),
    )

    # `SparseMatrixCSC .+ Array` give `SparseMatrixCSC` so we cast it with `Matrix`
    # before comparing.
    @test _matrix_equal(A .+ x, B .+ x)
    @test _matrix_equal(A .+ x, A .+ y)
    @test _matrix_equal(A .+ y, B .+ y)
    @test _matrix_equal(x .+ A, x .+ B)
    @test _matrix_equal(x .+ A, y .+ A)
    @test _matrix_equal(y .+ A, y .+ B)

    @test _matrix_equal(A .- x, B .- x)
    @test _matrix_equal(A .- x, A .- y)
    @test _matrix_equal(A .- y, B .- y)
    @test _matrix_equal(x .- A, x .- B)
    @test _matrix_equal(x .- A, y .- A)
    @test _matrix_equal(y .- A, y .- B)
end
function broadcast_test(x)
    if !(x isa AbstractMatrix)
        return
    end
    A = reshape(1:length(x), size(x)...)
    if size(x) == (2, 2)
        @test MA.isequal_canonical(
            A .+ x,
            [
                1+x[1, 1] 3+x[1, 2]
                2+x[2, 1] 4+x[2, 2]
            ],
        )
        @test MA.isequal_canonical(
            x .+ A,
            [
                1+x[1, 1] 3+x[1, 2]
                2+x[2, 1] 4+x[2, 2]
            ],
        )
        @test MA.isequal_canonical(
            x .+ x,
            [2x[1, 1] 2x[1, 2]; 2x[2, 1] 2x[2, 2]],
        )

        @test MA.isequal_canonical(
            A .- x,
            [
                1-x[1, 1] 3-x[1, 2]
                2-x[2, 1] 4-x[2, 2]
            ],
        )
        @test MA.isequal_canonical(
            x .- x,
            [zero(typeof(x[1] - x[1])) for _1 in 1:2, _2 in 1:2],
        )
        @test MA.isequal_canonical(
            x .- A,
            [
                -1+x[1, 1] -3+x[1, 2]
                -2+x[2, 1] -4+x[2, 2]
            ],
        )
    end
    return _broadcast_test(x, A)
end

function _broadcast_multiplication_test(x, A)
    B = sparse(A)
    y = SparseMatrixCSC(
        size(x)...,
        copy(B.colptr),
        copy(B.rowval),
        collect(vec(x)),
    )

    @test _matrix_equal(A .* x, B .* x)
    @test _matrix_equal(A .* x, A .* y)
    @test _matrix_equal(A .* y, B .* y)
    @test _matrix_equal(x .* A, x .* B)
    @test _matrix_equal(x .* A, y .* A)
    @test _matrix_equal(y .* A, y .* B)
end

function broadcast_multiplication_test(x)
    if !(x isa AbstractMatrix)
        return
    end
    A = reshape(1:length(x), size(x)...)
    if size(x) == (2, 2)
        @test MA.isequal_canonical(
            A .* x,
            [
                1*x[1, 1] 3*x[1, 2]
                2*x[2, 1] 4*x[2, 2]
            ],
        )
        @test MA.isequal_canonical(
            x .* A,
            [
                1*x[1, 1] 3*x[1, 2]
                2*x[2, 1] 4*x[2, 2]
            ],
        )
        @test MA.isequal_canonical(
            x .* x,
            [x[1, 1]^2 x[1, 2]^2; x[2, 1]^2 x[2, 2]^2],
        )

        # TODO: Refactor to avoid calling the internal JuMP function
        # `_densify_with_jump_eltype`.
        #z = JuMP._densify_with_jump_eltype((2 .* y) ./ 3)
        #@test MA.isequal_canonical((2 .* x) ./ 3, z)
        #z = JuMP._densify_with_jump_eltype(2 * (y ./ 3))
        #@test MA.isequal_canonical(2 .* (x ./ 3), z)
        #z = JuMP._densify_with_jump_eltype((x[1,1],) .* B)
        #@test MA.isequal_canonical((x[1,1],) .* A, z)
    end
    return _broadcast_multiplication_test(x, A)
end

function _broadcast_division_test(x, A)
    B = sparse(A)
    y = SparseMatrixCSC(
        size(x)...,
        copy(B.colptr),
        copy(B.rowval),
        collect(vec(x)),
    )

    @test _matrix_equal(x ./ A, x ./ B)
    @test _matrix_equal(x ./ A, y ./ A)
end
function broadcast_division_test(x)
    if !(x isa AbstractMatrix)
        return
    end
    A = reshape(1:length(x), size(x)...)
    if size(x) == (2, 2)
        @test MA.isequal_canonical(
            x ./ A,
            [
                x[1, 1]/1 x[1, 2]/3
                x[2, 1]/2 x[2, 2]/4
            ],
        )
    end
    return _broadcast_division_test(x, A)
end

function symmetric_unary_test(x)
    if x isa AbstractMatrix && size(x, 1) == size(x, 2)
        unary_test(LinearAlgebra.Symmetric(x))
    end
end

function symmetric_add_test(x)
    if x isa AbstractMatrix && size(x, 1) == size(x, 2)
        y = LinearAlgebra.Symmetric(x)
        add_test(y, y)
        add_test(x, y)
    end
end

function matrix_uniform_scaling_test(x)
    if !(x isa AbstractMatrix && size(x, 1) == size(x, 2))
        return
    end
    add_test(x, I)
    @test_rewrite (x .+ 1) + I
    @test_rewrite (x .- 1) - I
    @test_rewrite I + (x .+ 1)
    @test_rewrite I - (x .- 1)
    @test_rewrite I * x
    @test_rewrite I * (x .+ 1)
    @test_rewrite (x .+ 1) * I
    @test_rewrite (x .+ 1) + I * I
    @test_rewrite (x .+ 1) + 2 * I
    @test_rewrite (x .+ 1) + I * 2
    @test_rewrite I + 2I * x
    @test_rewrite x + 2I * x
end

function symmetric_matrix_uniform_scaling_test(x)
    if x isa AbstractMatrix && size(x, 1) == size(x, 2)
        matrix_uniform_scaling_test(LinearAlgebra.Symmetric(x))
    end
end

function triangular_test(x)
    if !(x isa AbstractMatrix && size(x, 1) == size(x, 2))
        return
    end
    n = LinearAlgebra.checksquare(x)
    ut = LinearAlgebra.UpperTriangular(x)
    add_test(ut, ut)
    y = Matrix(ut)
    for i in 1:n
        for j in 1:(i-1)
            @test iszero(y[i, j])
            @test MA.iszero!!(y[i, j])
        end
    end
    lt = LinearAlgebra.LowerTriangular(x)
    add_test(lt, lt)
    z = Matrix(lt)
    for j in 1:n
        for i in 1:(j-1)
            @test iszero(z[i, j])
            @test MA.iszero!!(z[i, j])
        end
    end
    return
end

function diagonal_test(x)
    if !(x isa AbstractVector && MA._one_indexed(x))
        return
    end
    d = LinearAlgebra.Diagonal(x)
    add_test(d, d)
    y = Matrix(d)
    t = LinearAlgebra.Tridiagonal(x[2:end], x, x[2:end])
    add_test(t, t)
    z = Matrix(t)
    for i in eachindex(x)
        @test MA.isequal_canonical(y[i, i], convert(eltype(y), x[i]))
        @test MA.isequal_canonical(z[i, i], convert(eltype(z), x[i]))
    end
    n = length(x)
    for j in 1:n
        for i in 1:(j-1)
            @test iszero(y[i, j])
            @test MA.iszero!!(y[i, j])
            @test iszero(y[j, i])
            @test MA.iszero!!(y[j, i])
            if abs(i - j) > 1
                @test iszero(z[i, j])
                @test MA.iszero!!(z[i, j])
                @test iszero(z[j, i])
                @test MA.iszero!!(z[j, i])
            end
        end
    end
end

const array_tests = Dict(
    "matrix_vector_division" => matrix_vector_division_test,
    "non_array" => non_array_test,
    "matrix_vector" => matrix_vector_test,
    "dot" => dot_test,
    "sum" => sum_test,
    "sum_multiplication" => sum_multiplication_test,
    "transpose" => transpose_test,
    "broadcast" => broadcast_test,
    "broadcast_multiplication" => broadcast_multiplication_test,
    "broadcast_division" => broadcast_division_test,
    "unary" => unary_test,
    "symmetric_unary" => symmetric_unary_test,
    "symmetric_add" => symmetric_add_test,
    "matrix_uniform_scaling" => matrix_uniform_scaling_test,
    "symmetric_matrix_uniform_scaling" =>
        symmetric_matrix_uniform_scaling_test,
    "triangular" => triangular_test,
    "diagonal" => diagonal_test,
)

@test_suite array

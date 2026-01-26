# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function _mat_mat_test(A, B)
    @test _matrix_equal(copy(transpose(A)) * B, A' * B)
    @test _matrix_equal(A * copy(transpose(B)), A * B')
    @test _matrix_equal(copy(transpose(B)) * A, B' * A)
    @test _matrix_equal(B * copy(transpose(A)), B * A')
end

function sparse_linear_test(X11, X23, Xd)
    X = sparse([1, 2], [1, 3], [X11, X23], 3, 3)

    #@test MA.isequal_canonical([X11 0. 0.; 0. 0. X23; 0. 0. 0.], @inferred MA._densify_with_jump_eltype(X))

    Y = sparse([1, 2], [1, 3], [2X11, 4X23], 3, 3)
    Yd = [
        2X11 0 0
        0 0 4X23
        0 0 0
    ]

    function _test_broadcast(A, B)
        @test_rewrite(A .+ B)
        @test_rewrite(B .+ A)
        @test_rewrite(A .- B)
        @test_rewrite(B .- A)
    end

    _test_broadcast(Xd, X)
    _test_broadcast(Xd, Y)
    _test_broadcast(Yd, X)
    _test_broadcast(Yd, Y)
    _test_broadcast(X, Y)

    add_test(Xd, Yd)
    add_test(Xd, Y)
    add_test(Xd, Xd)

    v = [4, 5, 6]
    @test MA.isequal_canonical(X * v, [4X11, 6X23, 0])
    @test MA.isequal_canonical((v' * X)', [4X11, 0, 5X23])
    @test MA.isequal_canonical(
        transpose(copy(transpose(v)) * X),
        [4X11, 0, 5X23],
    )
    @test MA.isequal_canonical(X' * v, [4X11, 0, 5X23])
    @test MA.isequal_canonical(copy(transpose(X)) * v, [4X11, 0, 5X23])

    A = [
        2 1 0
        1 2 1
        0 1 2
    ]

    add_test(X, A)
    _mat_mat_test(X, A)
    add_test(Y, A)
    _mat_mat_test(Y, A)

    @test _matrix_equal(
        X * A,
        [
            2X11 X11 0
            0 X23 2X23
            0 0 0
        ],
    )
    @test _matrix_equal(
        A * X,
        [
            2X11 0 X23
            X11 0 2X23
            0 0 X23
        ],
    )
    @test _matrix_equal(
        A * X',
        [
            2X11 0 0
            X11 X23 0
            0 2X23 0
        ],
    )
    @test _matrix_equal(
        X' * A,
        [
            2X11 X11 0
            0 0 0
            X23 2X23 X23
        ],
    )
    @test _matrix_equal(
        A' * X,
        [
            2X11 0 X23
            X11 0 2X23
            0 0 X23
        ],
    )

    B = sparse(A)
    @test _matrix_equal(X * A, X * B)
    @test _matrix_equal(A * X, B * X)
    @test _matrix_equal(A * X', B * X')
    @test _matrix_equal(A' * X, B' * X)
end

function sparse_quadratic_test(X11, X23, Xd)
    Y = sparse([1, 2], [1, 3], [2X11, 4X23], 3, 3) # for testing GenericAffExpr
    Yd = [
        2X11 0 0
        0 0 4X23
        0 0 0
    ]

    Z = sparse([1, 2], [1, 3], [X11^2, 2X23^2], 3, 3)
    Zd = [
        X11^2 0 0
        0 0 2X23^2
        0 0 0
    ]

    _mat_mat_test(Xd, Yd)
    _mat_mat_test(Xd, Y)
    _mat_mat_test(Xd, Xd)

    add_test(Yd, Zd)
    add_test(Zd, Xd)

    A = [
        2 1 0
        1 2 1
        0 1 2
    ]

    add_test(Z, A)
    return _mat_mat_test(Z, A)
end

const sparse_tests = Dict(
    "sparse_linear" => sparse_linear_test,
    "sparse_quadratic" => sparse_quadratic_test,
)

@test_suite sparse

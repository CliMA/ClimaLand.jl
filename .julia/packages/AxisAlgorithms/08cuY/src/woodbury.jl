function _A_ldiv_B_md!(dest, W::Woodbury, src,  R1, R2)
    _A_ldiv_B_md!(dest, W.A, src, R1, R2)
    tmp1 = _A_mul_B_md(W.V, dest, R1, R2)
    tmp2 = _A_mul_B_md(W.Cp, tmp1, R1, R2)
    tmp3 = _A_mul_B_md(W.U, tmp2, R1, R2)
    # TODO?: would be nice to fuse the next two steps
    tmp4 = _A_ldiv_B_md(W.A, tmp3, R1, R2)
    sub!(dest, tmp4)
end

function sub!(A, B)
    for I in eachindex(A, B)
        A[I] -= B[I]
    end
    A
end

check_matrix(W::Woodbury) = check_matrix(W.A)

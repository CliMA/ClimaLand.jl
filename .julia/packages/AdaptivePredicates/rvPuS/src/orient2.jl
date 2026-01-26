function orient2fast(pa, pb, pc)
    @inbounds begin
        acx = pa[1] - pc[1]
        bcx = pb[1] - pc[1]
        acy = pa[2] - pc[2]
        bcy = pb[2] - pc[2]
        return acx * bcy - acy * bcx
    end
end

function orient2exact(pa, pb, pc)
    @inbounds begin
        axby1, axby0 = Two_Product(pa[1], pb[2])
        axcy1, axcy0 = Two_Product(pa[1], pc[2])
        aterms3, aterms2, aterms1, aterms0 = Two_Two_Diff(axby1, axby0, axcy1, axcy0)
        aterms = (aterms0, aterms1, aterms2, aterms3)

        bxcy1, bxcy0 = Two_Product(pb[1], pc[2])
        bxay1, bxay0 = Two_Product(pb[1], pa[2])
        bterms3, bterms2, bterms1, bterms0 = Two_Two_Diff(bxcy1, bxcy0, bxay1, bxay0)
        bterms = (bterms0, bterms1, bterms2, bterms3)

        cxay1, cxay0 = Two_Product(pc[1], pa[2])
        cxby1, cxby0 = Two_Product(pc[1], pb[2])
        cterms3, cterms2, cterms1, cterms0 = Two_Two_Diff(cxay1, cxay0, cxby1, cxby0)
        cterms = (cterms0, cterms1, cterms2, cterms3)

        T = eltype(pa)
        h8 = ntuple(_ -> zero(T), Val(8))
        h12 = ntuple(_ -> zero(T), Val(12))
        v, vlength = fast_expansion_sum_zeroelim(4, aterms, 4, bterms, h8)
        w, wlength = fast_expansion_sum_zeroelim(vlength, v, 4, cterms, h12)
        return w[wlength]
    end
end

function orient2slow(pa, pb, pc)
    @inbounds begin
        acx, acxtail = Two_Diff(pa[1], pc[1])
        acy, acytail = Two_Diff(pa[2], pc[2])
        bcx, bcxtail = Two_Diff(pb[1], pc[1])
        bcy, bcytail = Two_Diff(pb[2], pc[2])

        axby7, axby6, axby5, axby4, axby3, axby2, axby1, axby0 = Two_Two_Product(acx, acxtail, bcy, bcytail)
        axby = (axby0, axby1, axby2, axby3, axby4, axby5, axby6, axby7)
        negate = -acy
        negatetail = -acytail
        bxay7, bxay6, bxay5, bxay4, bxay3, bxay2, bxay1, bxay0 = Two_Two_Product(bcx, bcxtail, negate, negatetail)
        bxay = (bxay0, bxay1, bxay2, bxay3, bxay4, bxay5, bxay6, bxay7)

        T = eltype(pa)
        h16 = ntuple(_ -> zero(T), Val(16))
        deter, deterlen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, h16)
        return deter[deterlen]
    end
end

function orient2adapt(pa, pb, pc, detsum)
    @inbounds begin
        T = eltype(pa)
        acx = pa[1] - pc[1]
        bcx = pb[1] - pc[1]
        acy = pa[2] - pc[2]
        bcy = pb[2] - pc[2]

        detleft, detlefttail = Two_Product(acx, bcy)
        detright, detrighttail = Two_Product(acy, bcx)

        B3, B2, B1, B0 = Two_Two_Diff(detleft, detlefttail, detright, detrighttail)
        B = (B0, B1, B2, B3)

        det = estimate(4, B)
        errbound = ccwerrboundB(T) * detsum
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        acxtail = Two_Diff_Tail(pa[1], pc[1], acx)
        bcxtail = Two_Diff_Tail(pb[1], pc[1], bcx)
        acytail = Two_Diff_Tail(pa[2], pc[2], acy)
        bcytail = Two_Diff_Tail(pb[2], pc[2], bcy)

        if iszero(acxtail) && iszero(acytail) && iszero(bcxtail) && iszero(bcytail)
            return det
        end

        errbound = ccwerrboundC(T) * detsum + resulterrbound(T) * Absolute(det)
        det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail)
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        s1, s0 = Two_Product(acxtail, bcy)
        t1, t0 = Two_Product(acytail, bcx)
        u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
        u = (u0, u1, u2, u3)
        h8 = ntuple(_ -> zero(T), Val(8))
        C1, C1length = fast_expansion_sum_zeroelim(4, B, 4, u, h8)

        s1, s0 = Two_Product(acx, bcytail)
        t1, t0 = Two_Product(acy, bcxtail)
        u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
        u = (u0, u1, u2, u3)
        h12 = ntuple(_ -> zero(T), Val(12))
        C2, C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, h12)

        s1, s0 = Two_Product(acxtail, bcytail)
        t1, t0 = Two_Product(acytail, bcxtail)
        u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
        u = (u0, u1, u2, u3)
        h16 = ntuple(_ -> zero(T), Val(16))
        D, Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, h16)

        return D[Dlength]
    end
end

function orient2(pa, pb, pc)
    @inbounds begin
        T = eltype(pa)
        detleft = (pa[1] - pc[1]) * (pb[2] - pc[2])
        detright = (pa[2] - pc[2]) * (pb[1] - pc[1])
        det = detleft - detright

        if detleft > 0
            if detright ≤ 0
                return det
            else
                detsum = detleft + detright
            end
        elseif detleft < 0
            if detright ≥ 0
                return det
            else
                detsum = -detleft - detright
            end
        else
            return det
        end

        errbound = ccwerrboundA(T) * detsum
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        return orient2adapt(pa, pb, pc, detsum)
    end
end
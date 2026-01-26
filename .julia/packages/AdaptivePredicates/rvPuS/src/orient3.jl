function orient3fast(pa, pb, pc, pd)
    @inbounds begin
        adx = pa[1] - pd[1]
        bdx = pb[1] - pd[1]
        cdx = pc[1] - pd[1]
        ady = pa[2] - pd[2]
        bdy = pb[2] - pd[2]
        cdy = pc[2] - pd[2]
        adz = pa[3] - pd[3]
        bdz = pb[3] - pd[3]
        cdz = pc[3] - pd[3]

        return adx * (bdy * cdz - bdz * cdy) +
               bdx * (cdy * adz - cdz * ady) +
               cdx * (ady * bdz - adz * bdy)
    end
end

function orient3exact(pa, pb, pc, pd, cache=nothing)
    @inbounds begin
        T = eltype(pa)
        h8 = ntuple(_ -> zero(T), Val(8))
        h12 = ntuple(_ -> zero(T), Val(12))
        h24 = ntuple(_ -> zero(T), Val(24))
        h32 = ntuple(_ -> zero(T), Val(32))

        axby1, axby0 = Two_Product(pa[1], pb[2])
        bxay1, bxay0 = Two_Product(pb[1], pa[2])
        ab3, ab2, ab1, ab0 = Two_Two_Diff(axby1, axby0, bxay1, bxay0)
        ab = (ab0, ab1, ab2, ab3)

        bxcy1, bxcy0 = Two_Product(pb[1], pc[2])
        cxby1, cxby0 = Two_Product(pc[1], pb[2])
        bc3, bc2, bc1, bc0 = Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0)
        bc = (bc0, bc1, bc2, bc3)

        cxdy1, cxdy0 = Two_Product(pc[1], pd[2])
        dxcy1, dxcy0 = Two_Product(pd[1], pc[2])
        cd3, cd2, cd1, cd0 = Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0)
        cd = (cd0, cd1, cd2, cd3)

        dxay1, dxay0 = Two_Product(pd[1], pa[2])
        axdy1, axdy0 = Two_Product(pa[1], pd[2])
        da3, da2, da1, da0 = Two_Two_Diff(dxay1, dxay0, axdy1, axdy0)
        da = (da0, da1, da2, da3)

        axcy1, axcy0 = Two_Product(pa[1], pc[2])
        cxay1, cxay0 = Two_Product(pc[1], pa[2])
        ac3, ac2, ac1, ac0 = Two_Two_Diff(axcy1, axcy0, cxay1, cxay0)
        ac = (ac0, ac1, ac2, ac3)

        bxdy1, bxdy0 = Two_Product(pb[1], pd[2])
        dxby1, dxby0 = Two_Product(pd[1], pb[2])
        bd3, bd2, bd1, bd0 = Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0)
        bd = (bd0, bd1, bd2, bd3)

        temp8, templen = fast_expansion_sum_zeroelim(4, cd, 4, da, h8)
        cda, cdalen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, h12)
        temp8, templen = fast_expansion_sum_zeroelim(4, da, 4, ab, h8)
        dab, dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, h12)
        bd = (-bd[1], -bd[2], -bd[3], -bd[4])
        ac = (-ac[1], -ac[2], -ac[3], -ac[4])
        temp8, templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, h8)
        abc, abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, h12)
        temp8, templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, h8)
        bcd, bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, h12)

        adet, alen = scale_expansion_zeroelim(bcdlen, bcd, pa[3], h24)
        bdet, blen = scale_expansion_zeroelim(cdalen, cda, -pb[3], h24)
        cdet, clen = scale_expansion_zeroelim(dablen, dab, pc[3], h24)
        ddet, dlen = scale_expansion_zeroelim(abclen, abc, -pd[3], h24)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h32) # abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h48_1)
        cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, h32) # cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, h48_2)
        deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, h32) # deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, h96)

        if ablen < 32 && cdlen < 32 && deterlen < 32
            return deter[deterlen]
        else
            h48_1, h48_2, h96 = orient3exact_cache(T, cache)

            abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h48_1)
            cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, h48_2)
            deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, h96)

            return deter[deterlen]
        end
    end
end

function orient3slow(pa, pb, pc, pd, cache=nothing)
    @inbounds begin
        T = eltype(pa)
        h16 = ntuple(_ -> zero(T), Val(16))
        h32 = ntuple(_ -> zero(T), Val(32))

        adx, adxtail = Two_Diff(pa[1], pd[1])
        ady, adytail = Two_Diff(pa[2], pd[2])
        adz, adztail = Two_Diff(pa[3], pd[3])
        bdx, bdxtail = Two_Diff(pb[1], pd[1])
        bdy, bdytail = Two_Diff(pb[2], pd[2])
        bdz, bdztail = Two_Diff(pb[3], pd[3])
        cdx, cdxtail = Two_Diff(pc[1], pd[1])
        cdy, cdytail = Two_Diff(pc[2], pd[2])
        cdz, cdztail = Two_Diff(pc[3], pd[3])

        axby7, axby6, axby5, axby4, axby3, axby2, axby1, axby0 = Two_Two_Product(adx, adxtail, bdy, bdytail)
        axby = (axby0, axby1, axby2, axby3, axby4, axby5, axby6, axby7)
        negate = -ady
        negatetail = -adytail
        bxay7, bxay6, bxay5, bxay4, bxay3, bxay2, bxay1, bxay0 = Two_Two_Product(bdx, bdxtail, negate, negatetail)
        bxay = (bxay0, bxay1, bxay2, bxay3, bxay4, bxay5, bxay6, bxay7)
        bxcy7, bxcy6, bxcy5, bxcy4, bxcy3, bxcy2, bxcy1, bxcy0 = Two_Two_Product(bdx, bdxtail, cdy, cdytail)
        bxcy = (bxcy0, bxcy1, bxcy2, bxcy3, bxcy4, bxcy5, bxcy6, bxcy7)
        negate = -bdy
        negatetail = -bdytail
        cxby7, cxby6, cxby5, cxby4, cxby3, cxby2, cxby1, cxby0 = Two_Two_Product(cdx, cdxtail, negate, negatetail)
        cxby = (cxby0, cxby1, cxby2, cxby3, cxby4, cxby5, cxby6, cxby7)
        cxay7, cxay6, cxay5, cxay4, cxay3, cxay2, cxay1, cxay0 = Two_Two_Product(cdx, cdxtail, ady, adytail)
        cxay = (cxay0, cxay1, cxay2, cxay3, cxay4, cxay5, cxay6, cxay7)
        negate = -cdy
        negatetail = -cdytail
        axcy7, axcy6, axcy5, axcy4, axcy3, axcy2, axcy1, axcy0 = Two_Two_Product(adx, adxtail, negate, negatetail)
        axcy = (axcy0, axcy1, axcy2, axcy3, axcy4, axcy5, axcy6, axcy7)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, h16)
        temp32, temp32len = scale_expansion_zeroelim(temp16len, temp16, adz, h32)
        temp32t, temp32tlen = scale_expansion_zeroelim(temp16len, temp16, adztail, h32)
        adet, alen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h32) # adet, alen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h64_1)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, h16)
        temp32, temp32len = scale_expansion_zeroelim(temp16len, temp16, bdz, h32)
        temp32t, temp32tlen = scale_expansion_zeroelim(temp16len, temp16, bdztail, h32)
        bdet, blen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h32) # bdet, blen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h64_2)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, h16)
        temp32, temp32len = scale_expansion_zeroelim(temp16len, temp16, cdz, h32)
        temp32t, temp32tlen = scale_expansion_zeroelim(temp16len, temp16, cdztail, h32)
        cdet, clen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h32) # cdet, clen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h64_3)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h32) # abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h128)
        deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h32) # deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h192_1)

        if alen < 32 && blen < 32 && clen < 32 && ablen < 32 && deterlen < 32
            return deter[deterlen]
        else
            h64_1, h64_2, h64_3, h128, h192_1 = orient3slow_cache(T, cache)

            temp16, temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, h16)
            temp32, temp32len = scale_expansion_zeroelim(temp16len, temp16, adz, h32)
            temp32t, temp32tlen = scale_expansion_zeroelim(temp16len, temp16, adztail, h32)
            adet, alen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h64_1)

            temp16, temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, h16)
            temp32, temp32len = scale_expansion_zeroelim(temp16len, temp16, bdz, h32)
            temp32t, temp32tlen = scale_expansion_zeroelim(temp16len, temp16, bdztail, h32)
            bdet, blen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h64_2)

            temp16, temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, h16)
            temp32, temp32len = scale_expansion_zeroelim(temp16len, temp16, cdz, h32)
            temp32t, temp32tlen = scale_expansion_zeroelim(temp16len, temp16, cdztail, h32)
            cdet, clen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, h64_3)

            abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h128)
            deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h192_1)

            return deter[deterlen]
        end
    end
end

function orient3adapt(pa, pb, pc, pd, permanent, cache=nothing)
    @inbounds begin
        T = eltype(pa)
        h8 = ntuple(_ -> zero(T), Val(8))
        h12 = ntuple(_ -> zero(T), Val(12))
        h16 = ntuple(_ -> zero(T), Val(16))
        h32 = ntuple(_ -> zero(T), Val(32))

        adx = pa[1] - pd[1]
        bdx = pb[1] - pd[1]
        cdx = pc[1] - pd[1]
        ady = pa[2] - pd[2]
        bdy = pb[2] - pd[2]
        cdy = pc[2] - pd[2]
        adz = pa[3] - pd[3]
        bdz = pb[3] - pd[3]
        cdz = pc[3] - pd[3]

        bdxcdy1, bdxcdy0 = Two_Product(bdx, cdy)
        cdxbdy1, cdxbdy0 = Two_Product(cdx, bdy)
        bc3, bc2, bc1, bc0 = Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0)
        bc = (bc0, bc1, bc2, bc3)
        adet, alen = scale_expansion_zeroelim(4, bc, adz, h8)

        cdxady1, cdxady0 = Two_Product(cdx, ady)
        adxcdy1, adxcdy0 = Two_Product(adx, cdy)
        ca3, ca2, ca1, ca0 = Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0)
        ca = (ca0, ca1, ca2, ca3)
        bdet, blen = scale_expansion_zeroelim(4, ca, bdz, h8)

        adxbdy1, adxbdy0 = Two_Product(adx, bdy)
        bdxady1, bdxady0 = Two_Product(bdx, ady)
        ab3, ab2, ab1, ab0 = Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0)
        ab = (ab0, ab1, ab2, ab3)
        cdet, clen = scale_expansion_zeroelim(4, ab, cdz, h8)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h16)

        args = pa, pb, pc, pd, permanent, adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz, bc, adet, alen, ca, bdet, blen, ab, cdet, clen, abdet, ablen, cache

        @check_length _orient3dadapt args fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h32) # fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h192_1)

        det = estimate(finlength, fin1)
        errbound = o3derrboundB(T) * permanent
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        adxtail = Two_Diff_Tail(pa[1], pd[1], adx)
        bdxtail = Two_Diff_Tail(pb[1], pd[1], bdx)
        cdxtail = Two_Diff_Tail(pc[1], pd[1], cdx)
        adytail = Two_Diff_Tail(pa[2], pd[2], ady)
        bdytail = Two_Diff_Tail(pb[2], pd[2], bdy)
        cdytail = Two_Diff_Tail(pc[2], pd[2], cdy)
        adztail = Two_Diff_Tail(pa[3], pd[3], adz)
        bdztail = Two_Diff_Tail(pb[3], pd[3], bdz)
        cdztail = Two_Diff_Tail(pc[3], pd[3], cdz)

        if iszero(adxtail) && iszero(bdxtail) && iszero(cdxtail) &&
           iszero(adytail) && iszero(bdytail) && iszero(cdytail) &&
           iszero(adztail) && iszero(bdztail) && iszero(cdztail)
            return det
        end

        errbound = o3derrboundC(T) * permanent + resulterrbound(T) * Absolute(det)
        det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
                adztail * (bdx * cdy - bdy * cdx)) +
               (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
                bdztail * (cdx * ady - cdy * adx)) +
               (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
                cdztail * (adx * bdy - ady * bdx))
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        finnow = fin1
        finother = h32 # finother = h192_2

        if iszero(adxtail)
            if iszero(adytail)
                at_b = (zero(T), zero(T), zero(T), zero(T))
                at_blen = 1
                at_c = (zero(T), zero(T), zero(T), zero(T))
                at_clen = 1
            else
                negate = -adytail
                at_blarge, at_b0 = Two_Product(negate, bdx)
                at_b = (at_b0, at_blarge, zero(T), zero(T))
                at_blen = 2
                at_clarge, at_c0 = Two_Product(adytail, cdx)
                at_c = (at_c0, at_clarge, zero(T), zero(T))
                at_clen = 2
            end
        else
            if iszero(adytail)
                at_blarge, at_b0 = Two_Product(adxtail, bdy)
                at_b = (at_b0, at_blarge, zero(T), zero(T))
                at_blen = 2
                negate = -adxtail
                at_clarge, at_c0 = Two_Product(negate, cdy)
                at_c = (at_c0, at_clarge, zero(T), zero(T))
                at_clen = 2
            else
                adxt_bdy1, adxt_bdy0 = Two_Product(adxtail, bdy)
                adyt_bdx1, adyt_bdx0 = Two_Product(adytail, bdx)
                at_blarge, at_b2, at_b1, at_b0 = Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0)
                at_b = (at_b0, at_b1, at_b2, at_blarge)
                at_blen = 4
                adyt_cdx1, adyt_cdx0 = Two_Product(adytail, cdx)
                adxt_cdy1, adxt_cdy0 = Two_Product(adxtail, cdy)
                at_clarge, at_c2, at_c1, at_c0 = Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0)
                at_c = (at_c0, at_c1, at_c2, at_clarge)
                at_clen = 4
            end
        end
        if iszero(bdxtail)
            if iszero(bdytail)
                bt_c = (zero(T), zero(T), zero(T), zero(T))
                bt_clen = 1
                bt_a = (zero(T), zero(T), zero(T), zero(T))
                bt_alen = 1
            else
                negate = -bdytail
                bt_clarge, bt_c0 = Two_Product(negate, cdx)
                bt_c = (bt_c0, bt_clarge, zero(T), zero(T))
                bt_clen = 2
                bt_alarge, bt_a0 = Two_Product(bdytail, adx)
                bt_a = (bt_a0, bt_alarge, zero(T), zero(T))
                bt_alen = 2
            end
        else
            if iszero(bdytail)
                bt_clarge, bt_c0 = Two_Product(bdxtail, cdy)
                bt_c = (bt_c0, bt_clarge, zero(T), zero(T))
                bt_clen = 2
                negate = -bdxtail
                bt_alarge, bt_a0 = Two_Product(negate, ady)
                bt_a = (bt_a0, bt_alarge, zero(T), zero(T))
                bt_alen = 2
            else
                bdxt_cdy1, bdxt_cdy0 = Two_Product(bdxtail, cdy)
                bdyt_cdx1, bdyt_cdx0 = Two_Product(bdytail, cdx)
                bt_clarge, bt_c2, bt_c1, bt_c0 = Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0)
                bt_c = (bt_c0, bt_c1, bt_c2, bt_clarge)
                bt_clen = 4
                bdyt_adx1, bdyt_adx0 = Two_Product(bdytail, adx)
                bdxt_ady1, bdxt_ady0 = Two_Product(bdxtail, ady)
                bt_alarge, bt_a2, bt_a1, bt_a0 = Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0)
                bt_a = (bt_a0, bt_a1, bt_a2, bt_alarge)
                bt_alen = 4
            end
        end
        if iszero(cdxtail)
            if iszero(cdytail)
                ct_a = (zero(T), zero(T), zero(T), zero(T))
                ct_alen = 1
                ct_b = (zero(T), zero(T), zero(T), zero(T))
                ct_blen = 1
            else
                negate = -cdytail
                ct_alarge, ct_a0 = Two_Product(negate, adx)
                ct_a = (ct_a0, ct_alarge, zero(T), zero(T))
                ct_alen = 2
                ct_blarge, ct_b0 = Two_Product(cdytail, bdx)
                ct_b = (ct_b0, ct_blarge, zero(T), zero(T))
                ct_blen = 2
            end
        else
            if iszero(cdytail)
                ct_alarge, ct_a0 = Two_Product(cdxtail, ady)
                ct_a = (ct_a0, ct_alarge, zero(T), zero(T))
                ct_alen = 2
                negate = -cdxtail
                ct_blarge, ct_b0 = Two_Product(negate, bdy)
                ct_b = (ct_b0, ct_blarge, zero(T), zero(T))
                ct_blen = 2
            else
                cdxt_ady1, cdxt_ady0 = Two_Product(cdxtail, ady)
                cdyt_adx1, cdyt_adx0 = Two_Product(cdytail, adx)
                ct_alarge, ct_a2, ct_a1, ct_a0 = Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0)
                ct_a = (ct_a0, ct_a1, ct_a2, ct_alarge)
                ct_alen = 4
                cdyt_bdx1, cdyt_bdx0 = Two_Product(cdytail, bdx)
                cdxt_bdy1, cdxt_bdy0 = Two_Product(cdxtail, bdy)
                ct_blarge, ct_b2, ct_b1, ct_b0 = Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0)
                ct_b = (ct_b0, ct_b1, ct_b2, ct_blarge)
                ct_blen = 4
            end
        end
        bct, bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, h8)
        w, wlength = scale_expansion_zeroelim(bctlen, bct, adz, h16)
        @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
        finnow, finother = finother, finnow

        cat, catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, h8)
        w, wlength = scale_expansion_zeroelim(catlen, cat, bdz, w)
        @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
        finnow, finother = finother, finnow

        abt, abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, h8)
        w, wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w)
        @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
        finnow, finother = finother, finnow

        if !iszero(adztail)
            v, vlength = scale_expansion_zeroelim(4, bc, adztail, h12)
            @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdztail)
            v, vlength = scale_expansion_zeroelim(4, ca, bdztail, h12)
            @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdztail)
            v, vlength = scale_expansion_zeroelim(4, ab, cdztail, h12)
            @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother)
            finnow, finother = finother, finnow
        end

        if !iszero(adxtail)
            if !iszero(bdytail)
                adxt_bdyt1, adxt_bdyt0 = Two_Product(adxtail, bdytail)
                u3, u2, u1, u0 = Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz)
                u = (u0, u1, u2, u3)
                @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(cdztail)
                    u3, u2, u1, u0 = Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail)
                    u = (u0, u1, u2, u3)
                    @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
            if !iszero(cdytail)
                negate = -adxtail
                adxt_cdyt1, adxt_cdyt0 = Two_Product(negate, cdytail)
                u3, u2, u1, u0 = Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz)
                u = (u0, u1, u2, u3)
                @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(bdztail)
                    u3, u2, u1, u0 = Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail)
                    u = (u0, u1, u2, u3)
                    @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
        end
        if !iszero(bdxtail)
            if !iszero(cdytail)
                bdxt_cdyt1, bdxt_cdyt0 = Two_Product(bdxtail, cdytail)
                u3, u2, u1, u0 = Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz)
                u = (u0, u1, u2, u3)
                @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(adztail)
                    u3, u2, u1, u0 = Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail)
                    u = (u0, u1, u2, u3)
                    @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
            if !iszero(adytail)
                negate = -bdxtail
                bdxt_adyt1, bdxt_adyt0 = Two_Product(negate, adytail)
                u3, u2, u1, u0 = Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz)
                u = (u0, u1, u2, u3)
                @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(cdztail)
                    u3, u2, u1, u0 = Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail)
                    u = (u0, u1, u2, u3)
                    @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
        end
        if !iszero(cdxtail)
            if !iszero(adytail)
                cdxt_adyt1, cdxt_adyt0 = Two_Product(cdxtail, adytail)
                u3, u2, u1, u0 = Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz)
                u = (u0, u1, u2, u3)
                @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(bdztail)
                    u3, u2, u1, u0 = Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail)
                    u = (u0, u1, u2, u3)
                    @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
            if !iszero(bdytail)
                negate = -cdxtail
                cdxt_bdyt1, cdxt_bdyt0 = Two_Product(negate, bdytail)
                u3, u2, u1, u0 = Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz)
                u = (u0, u1, u2, u3)
                @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(adztail)
                    u3, u2, u1, u0 = Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail)
                    u = (u0, u1, u2, u3)
                    @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
        end

        if !iszero(adztail)
            w, wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w)
            @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdztail)
            w, wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w)
            @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdztail)
            w, wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w)
            @check_length _orient3dadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
            finnow, finother = finother, finnow
        end

        return finnow[finlength]
    end
end
function _orient3dadapt(pa, pb, pc, pd, permanent, adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz, bc, adet, alen, ca, bdet, blen, ab, cdet, clen, abdet, ablen, cache)
    @inbounds begin
        T = eltype(pa)
        h8 = ntuple(_ -> zero(T), Val(8))
        h12 = ntuple(_ -> zero(T), Val(12))
        h16 = ntuple(_ -> zero(T), Val(16))
        h192_1, h192_2 = orient3adapt_cache(T, cache)

        fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h192_1) 

        det = estimate(finlength, fin1)
        errbound = o3derrboundB(T) * permanent
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        adxtail = Two_Diff_Tail(pa[1], pd[1], adx)
        bdxtail = Two_Diff_Tail(pb[1], pd[1], bdx)
        cdxtail = Two_Diff_Tail(pc[1], pd[1], cdx)
        adytail = Two_Diff_Tail(pa[2], pd[2], ady)
        bdytail = Two_Diff_Tail(pb[2], pd[2], bdy)
        cdytail = Two_Diff_Tail(pc[2], pd[2], cdy)
        adztail = Two_Diff_Tail(pa[3], pd[3], adz)
        bdztail = Two_Diff_Tail(pb[3], pd[3], bdz)
        cdztail = Two_Diff_Tail(pc[3], pd[3], cdz)

        if iszero(adxtail) && iszero(bdxtail) && iszero(cdxtail) &&
           iszero(adytail) && iszero(bdytail) && iszero(cdytail) &&
           iszero(adztail) && iszero(bdztail) && iszero(cdztail)
            return det
        end

        errbound = o3derrboundC(T) * permanent + resulterrbound(T) * Absolute(det)
        det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
                adztail * (bdx * cdy - bdy * cdx)) +
               (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
                bdztail * (cdx * ady - cdy * adx)) +
               (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
                cdztail * (adx * bdy - ady * bdx))
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        finnow = fin1
        finother = h192_2

        if iszero(adxtail)
            if iszero(adytail)
                at_b = (zero(T), zero(T), zero(T), zero(T))
                at_blen = 1
                at_c = (zero(T), zero(T), zero(T), zero(T))
                at_clen = 1
            else
                negate = -adytail
                at_blarge, at_b0 = Two_Product(negate, bdx)
                at_b = (at_b0, at_blarge, zero(T), zero(T))
                at_blen = 2
                at_clarge, at_c0 = Two_Product(adytail, cdx)
                at_c = (at_c0, at_clarge, zero(T), zero(T))
                at_clen = 2
            end
        else
            if iszero(adytail)
                at_blarge, at_b0 = Two_Product(adxtail, bdy)
                at_b = (at_b0, at_blarge, zero(T), zero(T))
                at_blen = 2
                negate = -adxtail
                at_clarge, at_c0 = Two_Product(negate, cdy)
                at_c = (at_c0, at_clarge, zero(T), zero(T))
                at_clen = 2
            else
                adxt_bdy1, adxt_bdy0 = Two_Product(adxtail, bdy)
                adyt_bdx1, adyt_bdx0 = Two_Product(adytail, bdx)
                at_blarge, at_b2, at_b1, at_b0 = Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0)
                at_b = (at_b0, at_b1, at_b2, at_blarge)
                at_blen = 4
                adyt_cdx1, adyt_cdx0 = Two_Product(adytail, cdx)
                adxt_cdy1, adxt_cdy0 = Two_Product(adxtail, cdy)
                at_clarge, at_c2, at_c1, at_c0 = Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0)
                at_c = (at_c0, at_c1, at_c2, at_clarge)
                at_clen = 4
            end
        end
        if iszero(bdxtail)
            if iszero(bdytail)
                bt_c = (zero(T), zero(T), zero(T), zero(T))
                bt_clen = 1
                bt_a = (zero(T), zero(T), zero(T), zero(T))
                bt_alen = 1
            else
                negate = -bdytail
                bt_clarge, bt_c0 = Two_Product(negate, cdx)
                bt_c = (bt_c0, bt_clarge, zero(T), zero(T))
                bt_clen = 2
                bt_alarge, bt_a0 = Two_Product(bdytail, adx)
                bt_a = (bt_a0, bt_alarge, zero(T), zero(T))
                bt_alen = 2
            end
        else
            if iszero(bdytail)
                bt_clarge, bt_c0 = Two_Product(bdxtail, cdy)
                bt_c = (bt_c0, bt_clarge, zero(T), zero(T))
                bt_clen = 2
                negate = -bdxtail
                bt_alarge, bt_a0 = Two_Product(negate, ady)
                bt_a = (bt_a0, bt_alarge, zero(T), zero(T))
                bt_alen = 2
            else
                bdxt_cdy1, bdxt_cdy0 = Two_Product(bdxtail, cdy)
                bdyt_cdx1, bdyt_cdx0 = Two_Product(bdytail, cdx)
                bt_clarge, bt_c2, bt_c1, bt_c0 = Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0)
                bt_c = (bt_c0, bt_c1, bt_c2, bt_clarge)
                bt_clen = 4
                bdyt_adx1, bdyt_adx0 = Two_Product(bdytail, adx)
                bdxt_ady1, bdxt_ady0 = Two_Product(bdxtail, ady)
                bt_alarge, bt_a2, bt_a1, bt_a0 = Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0)
                bt_a = (bt_a0, bt_a1, bt_a2, bt_alarge)
                bt_alen = 4
            end
        end
        if iszero(cdxtail)
            if iszero(cdytail)
                ct_a = (zero(T), zero(T), zero(T), zero(T))
                ct_alen = 1
                ct_b = (zero(T), zero(T), zero(T), zero(T))
                ct_blen = 1
            else
                negate = -cdytail
                ct_alarge, ct_a0 = Two_Product(negate, adx)
                ct_a = (ct_a0, ct_alarge, zero(T), zero(T))
                ct_alen = 2
                ct_blarge, ct_b0 = Two_Product(cdytail, bdx)
                ct_b = (ct_b0, ct_blarge, zero(T), zero(T))
                ct_blen = 2
            end
        else
            if iszero(cdytail)
                ct_alarge, ct_a0 = Two_Product(cdxtail, ady)
                ct_a = (ct_a0, ct_alarge, zero(T), zero(T))
                ct_alen = 2
                negate = -cdxtail
                ct_blarge, ct_b0 = Two_Product(negate, bdy)
                ct_b = (ct_b0, ct_blarge, zero(T), zero(T))
                ct_blen = 2
            else
                cdxt_ady1, cdxt_ady0 = Two_Product(cdxtail, ady)
                cdyt_adx1, cdyt_adx0 = Two_Product(cdytail, adx)
                ct_alarge, ct_a2, ct_a1, ct_a0 = Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0)
                ct_a = (ct_a0, ct_a1, ct_a2, ct_alarge)
                ct_alen = 4
                cdyt_bdx1, cdyt_bdx0 = Two_Product(cdytail, bdx)
                cdxt_bdy1, cdxt_bdy0 = Two_Product(cdxtail, bdy)
                ct_blarge, ct_b2, ct_b1, ct_b0 = Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0)
                ct_b = (ct_b0, ct_b1, ct_b2, ct_blarge)
                ct_blen = 4
            end
        end
        bct, bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, h8)
        w, wlength = scale_expansion_zeroelim(bctlen, bct, adz, h16)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
        finnow, finother = finother, finnow

        cat, catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, h8)
        w, wlength = scale_expansion_zeroelim(catlen, cat, bdz, w)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
        finnow, finother = finother, finnow

        abt, abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, h8)
        w, wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
        finnow, finother = finother, finnow

        if !iszero(adztail)
            v, vlength = scale_expansion_zeroelim(4, bc, adztail, h12)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdztail)
            v, vlength = scale_expansion_zeroelim(4, ca, bdztail, h12)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdztail)
            v, vlength = scale_expansion_zeroelim(4, ab, cdztail, h12)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother)
            finnow, finother = finother, finnow
        end

        if !iszero(adxtail)
            if !iszero(bdytail)
                adxt_bdyt1, adxt_bdyt0 = Two_Product(adxtail, bdytail)
                u3, u2, u1, u0 = Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz)
                u = (u0, u1, u2, u3)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(cdztail)
                    u3, u2, u1, u0 = Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail)
                    u = (u0, u1, u2, u3)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
            if !iszero(cdytail)
                negate = -adxtail
                adxt_cdyt1, adxt_cdyt0 = Two_Product(negate, cdytail)
                u3, u2, u1, u0 = Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz)
                u = (u0, u1, u2, u3)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(bdztail)
                    u3, u2, u1, u0 = Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail)
                    u = (u0, u1, u2, u3)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
        end
        if !iszero(bdxtail)
            if !iszero(cdytail)
                bdxt_cdyt1, bdxt_cdyt0 = Two_Product(bdxtail, cdytail)
                u3, u2, u1, u0 = Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz)
                u = (u0, u1, u2, u3)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(adztail)
                    u3, u2, u1, u0 = Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail)
                    u = (u0, u1, u2, u3)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
            if !iszero(adytail)
                negate = -bdxtail
                bdxt_adyt1, bdxt_adyt0 = Two_Product(negate, adytail)
                u3, u2, u1, u0 = Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz)
                u = (u0, u1, u2, u3)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(cdztail)
                    u3, u2, u1, u0 = Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail)
                    u = (u0, u1, u2, u3)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
        end
        if !iszero(cdxtail)
            if !iszero(adytail)
                cdxt_adyt1, cdxt_adyt0 = Two_Product(cdxtail, adytail)
                u3, u2, u1, u0 = Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz)
                u = (u0, u1, u2, u3)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(bdztail)
                    u3, u2, u1, u0 = Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail)
                    u = (u0, u1, u2, u3)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
            if !iszero(bdytail)
                negate = -cdxtail
                cdxt_bdyt1, cdxt_bdyt0 = Two_Product(negate, bdytail)
                u3, u2, u1, u0 = Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz)
                u = (u0, u1, u2, u3)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                finnow, finother = finother, finnow
                if !iszero(adztail)
                    u3, u2, u1, u0 = Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail)
                    u = (u0, u1, u2, u3)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother)
                    finnow, finother = finother, finnow
                end
            end
        end

        if !iszero(adztail)
            w, wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdztail)
            w, wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdztail)
            w, wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother)
            finnow, finother = finother, finnow
        end

        return finnow[finlength]
    end
end

function orient3(pa, pb, pc, pd, cache=nothing)
    @inbounds begin
        T = eltype(pa)
        adx = pa[1] - pd[1]
        bdx = pb[1] - pd[1]
        cdx = pc[1] - pd[1]
        ady = pa[2] - pd[2]
        bdy = pb[2] - pd[2]
        cdy = pc[2] - pd[2]
        adz = pa[3] - pd[3]
        bdz = pb[3] - pd[3]
        cdz = pc[3] - pd[3]

        bdxcdy = bdx * cdy
        cdxbdy = cdx * bdy

        cdxady = cdx * ady
        adxcdy = adx * cdy

        adxbdy = adx * bdy
        bdxady = bdx * ady

        det = adz * (bdxcdy - cdxbdy) +
              bdz * (cdxady - adxcdy) +
              cdz * (adxbdy - bdxady)

        permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz) +
                    (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz) +
                    (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz)
        errbound = o3derrboundA(T) * permanent
        if (det > errbound) || (-det > errbound)
            return det
        end

        return orient3adapt(pa, pb, pc, pd, permanent, cache)
    end
end

orient3exact_cache(_, cache) = cache
orient3exact_cache(::Type{T}, ::Nothing) where {T} = orient3exact_cache(T)
function orient3exact_cache(::Type{T}) where {T}
    cache = Vec{T}(undef, 256)          # cache_size = (64 + 64) + (128)
    h48_1 = view(cache, 1:48)           # 0 .+ (1:48)
    h48_2 = view(cache, 65:112)         # 64 .+ (1:48)
    h96 = view(cache, 129:224)          # 128 .+ (1:96)
    return h48_1, h48_2, h96
end
orient3slow_cache(_, cache) = cache
orient3slow_cache(::Type{T}, ::Nothing) where {T} = orient3slow_cache(T)
function orient3slow_cache(::Type{T}) where {T}
    cache = Vec{T}(undef, 512)          # cache_size = (64 + 64 + 64) + (128) + (192)
    h64_1 = view(cache, 1:64)           # 0 .+ (1:64)
    h64_2 = view(cache, 65:128)         # 64 .+ (1:64)
    h64_3 = view(cache, 129:192)        # 128 .+ (1:64)
    h128 = view(cache, 193:320)         # 192 .+ (1:128)
    h192_1 = view(cache, 321:512)       # 320 .+ (1:192)
    return h64_1, h64_2, h64_3, h128, h192_1
end
orient3adapt_cache(_, cache) = cache
orient3adapt_cache(::Type{T}, ::Nothing) where {T} = orient3adapt_cache(T)
function orient3adapt_cache(::Type{T}) where {T}
    cache = Vec{T}(undef, 384)          # cache_size = (192 + 192)
    h192_1 = view(cache, 1:192)         # 0 .+ (1:192)
    h192_2 = view(cache, 193:384)       # 192 .+ (1:192)
    return h192_1, h192_2
end
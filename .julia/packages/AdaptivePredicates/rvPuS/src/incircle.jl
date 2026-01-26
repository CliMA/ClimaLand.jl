function incirclefast(pa, pb, pc, pd)
    @inbounds begin
        adx = pa[1] - pd[1]
        ady = pa[2] - pd[2]
        bdx = pb[1] - pd[1]
        bdy = pb[2] - pd[2]
        cdx = pc[1] - pd[1]
        cdy = pc[2] - pd[2]

        abdet = adx * bdy - bdx * ady
        bcdet = bdx * cdy - cdx * bdy
        cadet = cdx * ady - adx * cdy
        alift = adx * adx + ady * ady
        blift = bdx * bdx + bdy * bdy
        clift = cdx * cdx + cdy * cdy

        return alift * bcdet + blift * cadet + clift * abdet
    end
end

function incircleexact(pa, pb, pc, pd, cache=nothing)
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
        temp8, templen = fast_expansion_sum_zeroelim(4, da, 4, ab, temp8)
        dab, dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, h12)
        bd = (-bd[1], -bd[2], -bd[3], -bd[4])
        ac = (-ac[1], -ac[2], -ac[3], -ac[4])
        temp8, templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, temp8)
        abc, abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, h12)
        temp8, templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, temp8)
        bcd, bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, h12)

        args = pa, pb, pc, pd, cda, cdalen, dab, dablen, abc, abclen, bcd, bcdlen, cache

        det24x, xlen = scale_expansion_zeroelim(bcdlen, bcd, pa[1], h24)
        @check_length _incircleexact args det48x, xlen = scale_expansion_zeroelim(xlen, det24x, pa[1], h32) # det48x, xlen = scale_expansion_zeroelim(xlen, det24x, pa[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(bcdlen, bcd, pa[2], h24)
        @check_length _incircleexact args det48y, ylen = scale_expansion_zeroelim(ylen, det24y, pa[2], h32) # det48y, ylen = scale_expansion_zeroelim(ylen, det24y, pa[2], h48_2)
        @check_length _incircleexact args adet, alen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h32) # adet, alen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_1)

        det24x, xlen = scale_expansion_zeroelim(cdalen, cda, pb[1], h24)
        @check_length _incircleexact args det48x, xlen = scale_expansion_zeroelim(xlen, det24x, -pb[1], h32) # det48x, xlen = scale_expansion_zeroelim(xlen, det24x, -pb[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(cdalen, cda, pb[2], h24)
        @check_length _incircleexact args det48y, ylen = scale_expansion_zeroelim(ylen, det24y, -pb[2], h32) # det48y, ylen = scale_expansion_zeroelim(ylen, det24y, -pb[2], h48_2)
        @check_length _incircleexact args bdet, blen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h32) # bdet, blen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_2)

        det24x, xlen = scale_expansion_zeroelim(dablen, dab, pc[1], h24)
        @check_length _incircleexact args det48x, xlen = scale_expansion_zeroelim(xlen, det24x, pc[1], h32) # det48x, xlen = scale_expansion_zeroelim(xlen, det24x, pc[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(dablen, dab, pc[2], h24)
        @check_length _incircleexact args det48y, ylen = scale_expansion_zeroelim(ylen, det24y, pc[2], h32) # det48y, ylen = scale_expansion_zeroelim(ylen, det24y, pc[2], h48_2)
        @check_length _incircleexact args cdet, clen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h32) # cdet, clen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_3)

        det24x, xlen = scale_expansion_zeroelim(abclen, abc, pd[1], h24)
        @check_length _incircleexact args det48x, xlen = scale_expansion_zeroelim(xlen, det24x, -pd[1], h32) # det48x, xlen = scale_expansion_zeroelim(xlen, det24x, -pd[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(abclen, abc, pd[2], h24)
        @check_length _incircleexact args det48y, ylen = scale_expansion_zeroelim(ylen, det24y, -pd[2], h32) # det48y, ylen = scale_expansion_zeroelim(ylen, det24y, -pd[2], h48_2)
        @check_length _incircleexact args ddet, dlen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h32) # ddet, dlen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_4)

        @check_length _incircleexact args abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h32) # abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h192_1)
        @check_length _incircleexact args cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, h32) # cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, h192_2)
        @check_length _incircleexact args deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, h32) # deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, h384_1)

        return deter[deterlen]
    end
end
function _incircleexact(pa, pb, pc, pd, cda, cdalen, dab, dablen, abc, abclen, bcd, bcdlen, cache)
    @inbounds begin
        T = eltype(pa)
        h24 = ntuple(_ -> zero(T), Val(24))
        h48_1, h48_2, h96_1, h96_2, h96_3, h96_4, h192_1, h192_2, h384_1 = incircleexact_cache(T, cache)

        det24x, xlen = scale_expansion_zeroelim(bcdlen, bcd, pa[1], h24)
        det48x, xlen = scale_expansion_zeroelim(xlen, det24x, pa[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(bcdlen, bcd, pa[2], h24)
        det48y, ylen = scale_expansion_zeroelim(ylen, det24y, pa[2], h48_2)
        adet, alen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_1)

        det24x, xlen = scale_expansion_zeroelim(cdalen, cda, pb[1], h24)
        det48x, xlen = scale_expansion_zeroelim(xlen, det24x, -pb[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(cdalen, cda, pb[2], h24)
        det48y, ylen = scale_expansion_zeroelim(ylen, det24y, -pb[2], h48_2)
        bdet, blen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_2)

        det24x, xlen = scale_expansion_zeroelim(dablen, dab, pc[1], h24)
        det48x, xlen = scale_expansion_zeroelim(xlen, det24x, pc[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(dablen, dab, pc[2], h24)
        det48y, ylen = scale_expansion_zeroelim(ylen, det24y, pc[2], h48_2)
        cdet, clen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_3)

        det24x, xlen = scale_expansion_zeroelim(abclen, abc, pd[1], h24)
        det48x, xlen = scale_expansion_zeroelim(xlen, det24x, -pd[1], h48_1)
        det24y, ylen = scale_expansion_zeroelim(abclen, abc, pd[2], h24)
        det48y, ylen = scale_expansion_zeroelim(ylen, det24y, -pd[2], h48_2)
        ddet, dlen = fast_expansion_sum_zeroelim(xlen, det48x, ylen, det48y, h96_4)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h192_1)
        cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, h192_2)
        deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, h384_1)

        return deter[deterlen]
    end
end

function incircleslow(pa, pb, pc, pd, cache=nothing)
    @inbounds begin
        T = eltype(pa)
        h16 = ntuple(_ -> zero(T), Val(16))
        h32 = ntuple(_ -> zero(T), Val(32))

        adx, adxtail = Two_Diff(pa[1], pd[1])
        ady, adytail = Two_Diff(pa[2], pd[2])
        bdx, bdxtail = Two_Diff(pb[1], pd[1])
        bdy, bdytail = Two_Diff(pb[2], pd[2])
        cdx, cdxtail = Two_Diff(pc[1], pd[1])
        cdy, cdytail = Two_Diff(pc[2], pd[2])

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

        args = pa, adx, adxtail, ady, adytail, bdx, bdxtail, bdy, bdytail, cdx, cdxtail, cdy, cdytail, axby, bxay, bxcy, cxby, cxay, axcy, cache

        temp16, temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, h16)

        detx, xlen = scale_expansion_zeroelim(temp16len, temp16, adx, h32)
        @check_length _incircleslow args detxx, xxlen = scale_expansion_zeroelim(xlen, detx, adx, h32) # detxx, xxlen = scale_expansion_zeroelim(xlen, detx, adx, h64_1)
        detxt, xtlen = scale_expansion_zeroelim(temp16len, temp16, adxtail, h32)
        @check_length _incircleslow args detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, adx, h32) # detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, adx, h64_2)
        detxxt = @. T(2.0 * detxxt)
        @check_length _incircleslow args detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, adxtail, h32) # detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, adxtail, h64_3)
        @check_length _incircleslow args x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, h32) # x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, h128_1)
        @check_length _incircleslow args x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, h32) # x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, h192_1)

        dety, ylen = scale_expansion_zeroelim(temp16len, temp16, ady, h32)
        @check_length _incircleslow args detyy, yylen = scale_expansion_zeroelim(ylen, dety, ady, h32) # detyy, yylen = scale_expansion_zeroelim(ylen, dety, ady, h64_4)
        detyt, ytlen = scale_expansion_zeroelim(temp16len, temp16, adytail, h32)
        @check_length _incircleslow args detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, ady, h32) # detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, ady, h64_5)
        detyyt = @. T(2.0 * detyyt)
        @check_length _incircleslow args detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, adytail, h32) # detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, adytail, h64_6)
        @check_length _incircleslow args y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, h32) # y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, h128_2)
        @check_length _incircleslow args y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, h32) # y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, h192_2)

        @check_length _incircleslow args adet, alen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h32) # adet, alen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h384_1)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16)

        detx, xlen = scale_expansion_zeroelim(temp16len, temp16, bdx, detx)
        @check_length _incircleslow args detxx, xxlen = scale_expansion_zeroelim(xlen, detx, bdx, detxx)
        detxt, xtlen = scale_expansion_zeroelim(temp16len, temp16, bdxtail, detxt)
        @check_length _incircleslow args detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, bdx, detxxt)
        detxxt = @. T(2.0 * detxxt)
        @check_length _incircleslow args detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, bdxtail, detxtxt)
        @check_length _incircleslow args x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1)
        @check_length _incircleslow args x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2)

        dety, ylen = scale_expansion_zeroelim(temp16len, temp16, bdy, dety)
        @check_length _incircleslow args detyy, yylen = scale_expansion_zeroelim(ylen, dety, bdy, detyy)
        detyt, ytlen = scale_expansion_zeroelim(temp16len, temp16, bdytail, detyt)
        @check_length _incircleslow args detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, bdy, detyyt)
        detyyt = @. T(2.0 * detyyt)
        @check_length _incircleslow args detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, bdytail, detytyt)
        @check_length _incircleslow args y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1)
        @check_length _incircleslow args y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2)

        @check_length _incircleslow args bdet, blen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h32) # bdet, blen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h384_2)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16)

        detx, xlen = scale_expansion_zeroelim(temp16len, temp16, cdx, detx)
        @check_length _incircleslow args detxx, xxlen = scale_expansion_zeroelim(xlen, detx, cdx, detxx)
        @check_length _incircleslow args detxt, xtlen = scale_expansion_zeroelim(temp16len, temp16, cdxtail, detxt)
        @check_length _incircleslow args detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, cdx, detxxt)

        detxxt = @. T(2.0 * detxxt)
        @check_length _incircleslow args detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, cdxtail, detxtxt)
        @check_length _incircleslow args x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1)
        @check_length _incircleslow args x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2)

        dety, ylen = scale_expansion_zeroelim(temp16len, temp16, cdy, dety)
        @check_length _incircleslow args detyy, yylen = scale_expansion_zeroelim(ylen, dety, cdy, detyy)
        @check_length _incircleslow args detyt, ytlen = scale_expansion_zeroelim(temp16len, temp16, cdytail, detyt)
        @check_length _incircleslow args detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, cdy, detyyt)
        detyyt = @. T(2.0 * detyyt)
        @check_length _incircleslow args detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, cdytail, detytyt)
        @check_length _incircleslow args y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1)
        @check_length _incircleslow args y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2)

        @check_length _incircleslow args cdet, clen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h32) # cdet, clen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h384_3)

        @check_length _incircleslow args abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h32) # abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h768)
        @check_length _incircleslow args deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h32) # deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h1152_1)

        return deter[deterlen]
    end
end
function _incircleslow(pa, adx, adxtail, ady, adytail, bdx, bdxtail, bdy, bdytail, cdx, cdxtail, cdy, cdytail, axby, bxay, bxcy, cxby, cxay, axcy, cache)
    @inbounds begin
        T = eltype(pa)
        h16 = ntuple(_ -> zero(T), Val(16))
        h32 = ntuple(_ -> zero(T), Val(32))
        h64_1, h64_2, h64_3, h64_4, h64_5, h64_6, h128_1, h128_2, h192_1, h192_2, h384_1, h384_2, h384_3, h768, h1152_1 = incircleslow_cache(T, cache)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, h16)

        detx, xlen = scale_expansion_zeroelim(temp16len, temp16, adx, h32)
        detxx, xxlen = scale_expansion_zeroelim(xlen, detx, adx, h64_1)
        detxt, xtlen = scale_expansion_zeroelim(temp16len, temp16, adxtail, h32)
        detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, adx, h64_2)
        for i in 1:xxtlen
            detxxt[i] *= 2.0
        end
        detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, adxtail, h64_3)
        x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, h128_1)
        x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, h192_1)

        dety, ylen = scale_expansion_zeroelim(temp16len, temp16, ady, h32)
        detyy, yylen = scale_expansion_zeroelim(ylen, dety, ady, h64_4)
        detyt, ytlen = scale_expansion_zeroelim(temp16len, temp16, adytail, h32)
        detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, ady, h64_5)
        for i in 1:yytlen
            detyyt[i] *= 2.0
        end
        detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, adytail, h64_6)
        y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, h128_2)
        y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, h192_2)

        adet, alen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h384_1)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16)

        detx, xlen = scale_expansion_zeroelim(temp16len, temp16, bdx, detx)
        detxx, xxlen = scale_expansion_zeroelim(xlen, detx, bdx, detxx)
        detxt, xtlen = scale_expansion_zeroelim(temp16len, temp16, bdxtail, detxt)
        detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, bdx, detxxt)
        for i in 1:xxtlen
            detxxt[i] *= 2.0
        end
        detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, bdxtail, detxtxt)
        x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1)
        x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2)

        dety, ylen = scale_expansion_zeroelim(temp16len, temp16, bdy, dety)
        detyy, yylen = scale_expansion_zeroelim(ylen, dety, bdy, detyy)
        detyt, ytlen = scale_expansion_zeroelim(temp16len, temp16, bdytail, detyt)
        detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, bdy, detyyt)
        for i in 1:yytlen
            detyyt[i] *= 2.0
        end
        detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, bdytail, detytyt)
        y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1)
        y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2)

        bdet, blen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h384_2)

        temp16, temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16)

        detx, xlen = scale_expansion_zeroelim(temp16len, temp16, cdx, detx)
        detxx, xxlen = scale_expansion_zeroelim(xlen, detx, cdx, detxx)
        detxt, xtlen = scale_expansion_zeroelim(temp16len, temp16, cdxtail, detxt)
        detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, cdx, detxxt)
        for i in 1:xxtlen
            detxxt[i] *= 2.0
        end
        detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, cdxtail, detxtxt)
        x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1)
        x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2)

        dety, ylen = scale_expansion_zeroelim(temp16len, temp16, cdy, dety)
        detyy, yylen = scale_expansion_zeroelim(ylen, dety, cdy, detyy)
        detyt, ytlen = scale_expansion_zeroelim(temp16len, temp16, cdytail, detyt)
        detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, cdy, detyyt)
        for i in 1:yytlen
            detyyt[i] *= 2.0
        end
        detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, cdytail, detytyt)
        y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1)
        y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2)

        cdet, clen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, h384_3)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h768)
        deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h1152_1)

        return deter[deterlen]
    end
end

function incircleadapt(pa, pb, pc, pd, permanent, cache=nothing)
    @inbounds begin
        T = eltype(pa)
        h8 = ntuple(_ -> zero(T), Val(8))
        h16 = ntuple(_ -> zero(T), Val(16))
        h32 = ntuple(_ -> zero(T), Val(32))

        adx = pa[1] - pd[1]
        bdx = pb[1] - pd[1]
        cdx = pc[1] - pd[1]
        ady = pa[2] - pd[2]
        bdy = pb[2] - pd[2]
        cdy = pc[2] - pd[2]

        bdxcdy1, bdxcdy0 = Two_Product(bdx, cdy)
        cdxbdy1, cdxbdy0 = Two_Product(cdx, bdy)
        bc3, bc2, bc1, bc0 = Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0)
        bc = (bc0, bc1, bc2, bc3)
        axbc, axbclen = scale_expansion_zeroelim(4, bc, adx, h8)
        axxbc, axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, h16)
        aybc, aybclen = scale_expansion_zeroelim(4, bc, ady, h8)
        ayybc, ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, h16)
        adet, alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, h32)

        cdxady1, cdxady0 = Two_Product(cdx, ady)
        adxcdy1, adxcdy0 = Two_Product(adx, cdy)
        ca3, ca2, ca1, ca0 = Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0)
        ca = (ca0, ca1, ca2, ca3)
        bxca, bxcalen = scale_expansion_zeroelim(4, ca, bdx, h8)
        bxxca, bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, h16)
        byca, bycalen = scale_expansion_zeroelim(4, ca, bdy, h8)
        byyca, byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, h16)
        bdet, blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, h32)

        adxbdy1, adxbdy0 = Two_Product(adx, bdy)
        bdxady1, bdxady0 = Two_Product(bdx, ady)
        ab3, ab2, ab1, ab0 = Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0)
        ab = (ab0, ab1, ab2, ab3)
        cxab, cxablen = scale_expansion_zeroelim(4, ab, cdx, h8)
        cxxab, cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, h16)
        cyab, cyablen = scale_expansion_zeroelim(4, ab, cdy, h8)
        cyyab, cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, h16)
        cdet, clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, h32)

        args = pa, pb, pc, pd, alen, adet, blen, bdet, clen, cdet, permanent, adx, bdx, cdx, ady, bdy, cdy, bc, ca, ab, cache

        @check_length _incircleadapt args abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h32) # abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h64_1)
        @check_length _incircleadapt args fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h32) # fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h1152_1)

        det = estimate(finlength, fin1)
        errbound = iccerrboundB(T) * permanent
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        adxtail = Two_Diff_Tail(pa[1], pd[1], adx)
        adytail = Two_Diff_Tail(pa[2], pd[2], ady)
        bdxtail = Two_Diff_Tail(pb[1], pd[1], bdx)
        bdytail = Two_Diff_Tail(pb[2], pd[2], bdy)
        cdxtail = Two_Diff_Tail(pc[1], pd[1], cdx)
        cdytail = Two_Diff_Tail(pc[2], pd[2], cdy)

        if iszero(adxtail) && iszero(bdxtail) && iszero(cdxtail) &&
           iszero(adytail) && iszero(bdytail) && iszero(cdytail)
            return det
        end

        errbound = iccerrboundC(T) * permanent + resulterrbound(T) * Absolute(det)
        detadd = ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) -
                                             (bdy * cdxtail + cdx * bdytail)) +
                  2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
                 ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) -
                                             (cdy * adxtail + adx * cdytail)) +
                  2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
                 ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) -
                                             (ady * bdxtail + bdx * adytail)) +
                  2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx))
        det = T(det + detadd) # Had to change this to match how C handles the 2.0 multiplication with Float32
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        finnow = fin1
        finother = h32 # finother = h1152_2

        if !iszero(bdxtail) || !iszero(bdytail) || !iszero(cdxtail) || !iszero(cdytail)
            adxadx1, adxadx0 = Square(adx)
            adyady1, adyady0 = Square(ady)
            aa3, aa2, aa1, aa0 = Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0)
            aa = (aa0, aa1, aa2, aa3)
        end
        if !iszero(cdxtail) || !iszero(cdytail) || !iszero(adxtail) || !iszero(adytail)
            bdxbdx1, bdxbdx0 = Square(bdx)
            bdybdy1, bdybdy0 = Square(bdy)
            bb3, bb2, bb1, bb0 = Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0)
            bb = (bb0, bb1, bb2, bb3)
        end
        if !iszero(adxtail) || !iszero(adytail) || !iszero(bdxtail) || !iszero(bdytail)
            cdxcdx1, cdxcdx0 = Square(cdx)
            cdycdy1, cdycdy0 = Square(cdy)
            cc3, cc2, cc1, cc0 = Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0)
            cc = (cc0, cc1, cc2, cc3)
        end

        if !iszero(adxtail)
            axtbc, axtbclen = scale_expansion_zeroelim(4, bc, adxtail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2adx, h16)

            axtcc, axtcclen = scale_expansion_zeroelim(4, cc, adxtail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, h16)

            axtbb, axtbblen = scale_expansion_zeroelim(4, bb, adxtail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(adytail)
            aytbc, aytbclen = scale_expansion_zeroelim(4, bc, adytail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2 * ady, h16)

            aytbb, aytbblen = scale_expansion_zeroelim(4, bb, adytail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, h16)

            aytcc, aytcclen = scale_expansion_zeroelim(4, cc, adytail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdxtail)
            bxtca, bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2 * bdx, h16)

            bxtaa, bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, h16)

            bxtcc, bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdytail)
            bytca, bytcalen = scale_expansion_zeroelim(4, ca, bdytail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2 * bdy, h16)

            bytcc, bytcclen = scale_expansion_zeroelim(4, cc, bdytail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, h16)

            bytaa, bytaalen = scale_expansion_zeroelim(4, aa, bdytail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdxtail)
            cxtab, cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2 * cdx, h16)

            cxtbb, cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, h16)

            cxtaa, cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdytail)
            cytab, cytablen = scale_expansion_zeroelim(4, ab, cdytail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2 * cdy, h16)

            cytaa, cytaalen = scale_expansion_zeroelim(4, aa, cdytail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, h16)

            cytbb, cytbblen = scale_expansion_zeroelim(4, bb, cdytail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end

        if !iszero(adxtail) || !iszero(adytail)
            if !iszero(bdxtail) || !iszero(bdytail) || !iszero(cdxtail) || !iszero(cdytail)
                ti1, ti0 = Two_Product(bdxtail, cdy)
                tj1, tj0 = Two_Product(bdx, cdytail)
                u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                u = (u0, u1, u2, u3)
                negate = -bdy
                ti1, ti0 = Two_Product(cdxtail, negate)
                negate = -bdytail
                tj1, tj0 = Two_Product(cdx, negate)
                v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                v = (v0, v1, v2, v3)
                bct, bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, h8)

                ti1, ti0 = Two_Product(bdxtail, cdytail)
                tj1, tj0 = Two_Product(cdxtail, bdytail)
                bctt3, bctt2, bctt1, bctt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
                bctt = (bctt0, bctt1, bctt2, bctt3)
                bcttlen = 4
            else
                bct = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
                bctlen = 1
                bctt = (zero(T), zero(T), zero(T), zero(T))
                bcttlen = 1
            end

            if !iszero(adxtail)
                temp16a, temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, h16)
                axtbct, axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2adx, h32)
                @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow
                if !iszero(bdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, cc, adxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, h16)
                    @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end
                if !iszero(cdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, bb, -adxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, h16)
                    @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end

                temp32a, temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail, h32)
                axtbctt, axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2adx, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                @check_length _incircleadapt args temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h32) # temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(adytail)
                temp16a, temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, h16)
                aytbct, aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2ady, h32)
                @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow

                temp32a, temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail, h32)
                aytbctt, aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2ady, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                @check_length _incircleadapt args temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h32) # temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
        end

        if !iszero(bdxtail) || !iszero(bdytail)
            if !iszero(cdxtail) || !iszero(cdytail) || !iszero(adxtail) || !iszero(adytail)
                ti1, ti0 = Two_Product(cdxtail, ady)
                tj1, tj0 = Two_Product(cdx, adytail)
                u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                u = (u0, u1, u2, u3)
                negate = -cdy
                ti1, ti0 = Two_Product(adxtail, negate)
                negate = -cdytail
                tj1, tj0 = Two_Product(adx, negate)
                v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                v = (v0, v1, v2, v3)
                cat, catlen = fast_expansion_sum_zeroelim(4, u, 4, v, h8)

                ti1, ti0 = Two_Product(cdxtail, adytail)
                tj1, tj0 = Two_Product(adxtail, cdytail)
                catt3, catt2, catt1, catt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
                catt = (catt0, catt1, catt2, catt3)
                cattlen = 4
            else
                cat = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
                catlen = 1
                catt = (zero(T), zero(T), zero(T), zero(T))
                cattlen = 1
            end

            if !iszero(bdxtail)
                temp16a, temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, h16)
                bxtcat, bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2bdx, h32)
                @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow
                if !iszero(cdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, aa, bdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, h16)
                    @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end
                if !iszero(adytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, h16)
                    @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end

                temp32a, temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail, h32)
                bxtcatt, bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2bdx, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                @check_length _incircleadapt args temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h32) # temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(bdytail)
                temp16a, temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, h16)
                bytcat, bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2bdy, h32)
                @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow

                temp32a, temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail, h32)
                bytcatt, bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2bdy, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                @check_length _incircleadapt args temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h32) # temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
        end

        if !iszero(cdxtail) || !iszero(cdytail)
            if !iszero(adxtail) || !iszero(adytail) || !iszero(bdxtail) || !iszero(bdytail)
                ti1, ti0 = Two_Product(adxtail, bdy)
                tj1, tj0 = Two_Product(adx, bdytail)
                u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                u = (u0, u1, u2, u3)
                negate = -ady
                ti1, ti0 = Two_Product(bdxtail, negate)
                negate = -adytail
                tj1, tj0 = Two_Product(bdx, negate)
                v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                v = (v0, v1, v2, v3)
                abt, abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, h8)

                ti1, ti0 = Two_Product(adxtail, bdytail)
                tj1, tj0 = Two_Product(bdxtail, adytail)
                abtt3, abtt2, abtt1, abtt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
                abtt = (abtt0, abtt1, abtt2, abtt3)
                abttlen = 4
            else
                abt = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
                abtlen = 1
                abtt = (zero(T), zero(T), zero(T), zero(T))
                abttlen = 1
            end

            if !iszero(cdxtail)
                temp16a, temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, h16)
                cxtabt, cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2cdx, h32)
                @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow
                if !iszero(adytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, bb, cdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, h16)
                    @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end
                if !iszero(bdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, h16)
                    @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end

                temp32a, temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail, h32)
                cxtabtt, cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2cdx, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                @check_length _incircleadapt args temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h32) # temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(cdytail)
                temp16a, temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, h16)
                cytabt, cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2cdy, h32)
                @check_length _incircleadapt args temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h32) # temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow

                temp32a, temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail, h32)
                cytabtt, cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2cdy, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                @check_length _incircleadapt args temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h32) # temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                @check_length _incircleadapt args finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
        end
        return finnow[finlength]
    end
end
function _incircleadapt(pa, pb, pc, pd, alen, adet, blen, bdet, clen, cdet, permanent, adx, bdx, cdx, ady, bdy, cdy, bc, ca, ab, cache)
    @inbounds begin
        T = eltype(pa)
        h8 = ntuple(_ -> zero(T), Val(8))
        h16 = ntuple(_ -> zero(T), Val(16))
        h32 = ntuple(_ -> zero(T), Val(32))
        h48_1, h64_1, h1152_1, h1152_2 = incircleadapt_cache(T, cache)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h64_1)
        fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, h1152_1)

        det = estimate(finlength, fin1)
        errbound = iccerrboundB(T) * permanent
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        adxtail = Two_Diff_Tail(pa[1], pd[1], adx)
        adytail = Two_Diff_Tail(pa[2], pd[2], ady)
        bdxtail = Two_Diff_Tail(pb[1], pd[1], bdx)
        bdytail = Two_Diff_Tail(pb[2], pd[2], bdy)
        cdxtail = Two_Diff_Tail(pc[1], pd[1], cdx)
        cdytail = Two_Diff_Tail(pc[2], pd[2], cdy)

        if iszero(adxtail) && iszero(bdxtail) && iszero(cdxtail) &&
           iszero(adytail) && iszero(bdytail) && iszero(cdytail)
            return det
        end

        errbound = iccerrboundC(T) * permanent + resulterrbound(T) * Absolute(det)
        detadd = ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) -
                                             (bdy * cdxtail + cdx * bdytail)) +
                  2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
                 ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) -
                                             (cdy * adxtail + adx * cdytail)) +
                  2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
                 ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) -
                                             (ady * bdxtail + bdx * adytail)) +
                  2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx))
        det = T(det + detadd) # Had to change this to match how C handles the 2.0 multiplication with Float32
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        finnow = fin1
        finother = h1152_2

        if !iszero(bdxtail) || !iszero(bdytail) || !iszero(cdxtail) || !iszero(cdytail)
            adxadx1, adxadx0 = Square(adx)
            adyady1, adyady0 = Square(ady)
            aa3, aa2, aa1, aa0 = Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0)
            aa = (aa0, aa1, aa2, aa3)
        end
        if !iszero(cdxtail) || !iszero(cdytail) || !iszero(adxtail) || !iszero(adytail)
            bdxbdx1, bdxbdx0 = Square(bdx)
            bdybdy1, bdybdy0 = Square(bdy)
            bb3, bb2, bb1, bb0 = Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0)
            bb = (bb0, bb1, bb2, bb3)
        end
        if !iszero(adxtail) || !iszero(adytail) || !iszero(bdxtail) || !iszero(bdytail)
            cdxcdx1, cdxcdx0 = Square(cdx)
            cdycdy1, cdycdy0 = Square(cdy)
            cc3, cc2, cc1, cc0 = Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0)
            cc = (cc0, cc1, cc2, cc3)
        end

        if !iszero(adxtail)
            axtbc, axtbclen = scale_expansion_zeroelim(4, bc, adxtail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2adx, h16)

            axtcc, axtcclen = scale_expansion_zeroelim(4, cc, adxtail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, h16)

            axtbb, axtbblen = scale_expansion_zeroelim(4, bb, adxtail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(adytail)
            aytbc, aytbclen = scale_expansion_zeroelim(4, bc, adytail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2 * ady, h16)

            aytbb, aytbblen = scale_expansion_zeroelim(4, bb, adytail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, h16)

            aytcc, aytcclen = scale_expansion_zeroelim(4, cc, adytail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdxtail)
            bxtca, bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2 * bdx, h16)

            bxtaa, bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, h16)

            bxtcc, bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdytail)
            bytca, bytcalen = scale_expansion_zeroelim(4, ca, bdytail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2 * bdy, h16)

            bytcc, bytcclen = scale_expansion_zeroelim(4, cc, bdytail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, h16)

            bytaa, bytaalen = scale_expansion_zeroelim(4, aa, bdytail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdxtail)
            cxtab, cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2 * cdx, h16)

            cxtbb, cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, h16)

            cxtaa, cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdytail)
            cytab, cytablen = scale_expansion_zeroelim(4, ab, cdytail, h8)
            temp16a, temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2 * cdy, h16)

            cytaa, cytaalen = scale_expansion_zeroelim(4, aa, cdytail, h8)
            temp16b, temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, h16)

            cytbb, cytbblen = scale_expansion_zeroelim(4, bb, cdytail, h8)
            temp16c, temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, h16)

            temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48_1)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
        end

        if !iszero(adxtail) || !iszero(adytail)
            if !iszero(bdxtail) || !iszero(bdytail) || !iszero(cdxtail) || !iszero(cdytail)
                ti1, ti0 = Two_Product(bdxtail, cdy)
                tj1, tj0 = Two_Product(bdx, cdytail)
                u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                u = (u0, u1, u2, u3)
                negate = -bdy
                ti1, ti0 = Two_Product(cdxtail, negate)
                negate = -bdytail
                tj1, tj0 = Two_Product(cdx, negate)
                v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                v = (v0, v1, v2, v3)
                bct, bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, h8)

                ti1, ti0 = Two_Product(bdxtail, cdytail)
                tj1, tj0 = Two_Product(cdxtail, bdytail)
                bctt3, bctt2, bctt1, bctt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
                bctt = (bctt0, bctt1, bctt2, bctt3)
                bcttlen = 4
            else
                bct = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
                bctlen = 1
                bctt = (zero(T), zero(T), zero(T), zero(T))
                bcttlen = 1
            end

            if !iszero(adxtail)
                temp16a, temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, h16)
                axtbct, axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2adx, h32)
                temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow
                if !iszero(bdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, cc, adxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, h16)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end
                if !iszero(cdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, bb, -adxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, h16)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end

                temp32a, temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail, h32)
                axtbctt, axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2adx, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(adytail)
                temp16a, temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, h16)
                aytbct, aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2ady, h32)
                temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow

                temp32a, temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail, h32)
                aytbctt, aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2ady, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
        end

        if !iszero(bdxtail) || !iszero(bdytail)
            if !iszero(cdxtail) || !iszero(cdytail) || !iszero(adxtail) || !iszero(adytail)
                ti1, ti0 = Two_Product(cdxtail, ady)
                tj1, tj0 = Two_Product(cdx, adytail)
                u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                u = (u0, u1, u2, u3)
                negate = -cdy
                ti1, ti0 = Two_Product(adxtail, negate)
                negate = -cdytail
                tj1, tj0 = Two_Product(adx, negate)
                v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                v = (v0, v1, v2, v3)
                cat, catlen = fast_expansion_sum_zeroelim(4, u, 4, v, h8)

                ti1, ti0 = Two_Product(cdxtail, adytail)
                tj1, tj0 = Two_Product(adxtail, cdytail)
                catt3, catt2, catt1, catt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
                catt = (catt0, catt1, catt2, catt3)
                cattlen = 4
            else
                cat = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
                catlen = 1
                catt = (zero(T), zero(T), zero(T), zero(T))
                cattlen = 1
            end 

            if !iszero(bdxtail)
                temp16a, temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, h16)
                bxtcat, bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2bdx, h32)
                temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow
                if !iszero(cdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, aa, bdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, h16)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end
                if !iszero(adytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, h16)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end

                temp32a, temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail, h32)
                bxtcatt, bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2bdx, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(bdytail)
                temp16a, temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, h16)
                bytcat, bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2bdy, h32)
                temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow

                temp32a, temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail, h32)
                bytcatt, bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2bdy, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
        end

        if !iszero(cdxtail) || !iszero(cdytail)
            if !iszero(adxtail) || !iszero(adytail) || !iszero(bdxtail) || !iszero(bdytail)
                ti1, ti0 = Two_Product(adxtail, bdy)
                tj1, tj0 = Two_Product(adx, bdytail)
                u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                u = (u0, u1, u2, u3)
                negate = -ady
                ti1, ti0 = Two_Product(bdxtail, negate)
                negate = -adytail
                tj1, tj0 = Two_Product(bdx, negate)
                v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
                v = (v0, v1, v2, v3)
                abt, abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, h8)

                ti1, ti0 = Two_Product(adxtail, bdytail)
                tj1, tj0 = Two_Product(bdxtail, adytail)
                abtt3, abtt2, abtt1, abtt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
                abtt = (abtt0, abtt1, abtt2, abtt3)
                abttlen = 4
            else
                abt = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
                abtlen = 1
                abtt = (zero(T), zero(T), zero(T), zero(T))
                abttlen = 1
            end

            if !iszero(cdxtail)
                temp16a, temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, h16)
                cxtabt, cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2cdx, h32)
                temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow
                if !iszero(adytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, bb, cdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, h16)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end
                if !iszero(bdytail)
                    temp8, temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, h8)
                    temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, h16)
                    finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                    finnow, finother = finother, finnow
                end

                temp32a, temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail, h32)
                cxtabtt, cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2cdx, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(cdytail)
                temp16a, temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, h16)
                cytabt, cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, h16)
                temp32a, temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2cdy, h32)
                temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
                finnow, finother = finother, finnow

                temp32a, temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail, h32)
                cytabtt, cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, h8)
                temp16a, temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2cdy, h16)
                temp16b, temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail, h16)
                temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, h32)
                temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64_1)
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
                finnow, finother = finother, finnow
            end
        end

        return finnow[finlength]
    end
end

function incircle(pa, pb, pc, pd, cache=nothing)
    @inbounds begin
        T = eltype(pa)
        adx = pa[1] - pd[1]
        bdx = pb[1] - pd[1]
        cdx = pc[1] - pd[1]
        ady = pa[2] - pd[2]
        bdy = pb[2] - pd[2]
        cdy = pc[2] - pd[2]

        bdxcdy = bdx * cdy
        cdxbdy = cdx * bdy
        alift = adx * adx + ady * ady

        cdxady = cdx * ady
        adxcdy = adx * cdy
        blift = bdx * bdx + bdy * bdy

        adxbdy = adx * bdy
        bdxady = bdx * ady
        clift = cdx * cdx + cdy * cdy

        det = alift * (bdxcdy - cdxbdy) +
              blift * (cdxady - adxcdy) +
              clift * (adxbdy - bdxady)

        permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift +
                    (Absolute(cdxady) + Absolute(adxcdy)) * blift +
                    (Absolute(adxbdy) + Absolute(bdxady)) * clift
        errbound = iccerrboundA(T) * permanent
        if (det > errbound) || (-det > errbound)
            return det
        end

        return incircleadapt(pa, pb, pc, pd, permanent, cache)
    end
end

incircleexact_cache(_, cache) = cache
incircleexact_cache(::Type{T}, ::Nothing) where {T} = incircleexact_cache(T)
function incircleexact_cache(::Type{T}) where {T}
    cache = Vec{T}(undef, 1408)         # cache_size = (64 + 64) + (128 + 128 + 128 + 128) + (192 + 192) + (384)
    h48_1 = view(cache, 1:48)           # 0 .+ (1:48)
    h48_2 = view(cache, 65:112)         # 64 .+ (1:48)
    h96_1 = view(cache, 129:224)        # 128 .+ (1:96)
    h96_2 = view(cache, 257:352)        # 256 .+ (1:96)
    h96_3 = view(cache, 385:480)        # 384 .+ (1:96)
    h96_4 = view(cache, 513:608)        # 512 .+ (1:96)
    h192_1 = view(cache, 641:832)       # 640 .+ (1:192)
    h192_2 = view(cache, 833:1024)      # 832 .+ (1:192)
    h384_1 = view(cache, 1025:1408)     # 1024 .+ (1:384)
    return h48_1, h48_2, h96_1, h96_2, h96_3, h96_4, h192_1, h192_2, h384_1
end
incircleslow_cache(_, cache) = cache
incircleslow_cache(::Type{T}, ::Nothing) where {T} = incircleslow_cache(T)
function incircleslow_cache(::Type{T}) where {T}
    cache = Vec{T}(undef, 4096)         # cache_size = (64 + 64 + 64 + 64 + 64 + 64) + (128 + 128) + (192 + 192) + (384 + 384 + 384) + (768) + (1152)
    h64_1 = view(cache, 1:64)           # 0 .+ (1:64)
    h64_2 = view(cache, 65:128)         # 64 .+ (1:64)
    h64_3 = view(cache, 129:192)        # 128 .+ (1:64)
    h64_4 = view(cache, 193:256)        # 192 .+ (1:64)
    h64_5 = view(cache, 257:320)        # 256 .+ (1:64)
    h64_6 = view(cache, 321:384)        # 320 .+ (1:64)
    h128_1 = view(cache, 385:512)       # 384 .+ (1:128)
    h128_2 = view(cache, 513:640)       # 512 .+ (1:128)
    h192_1 = view(cache, 641:832)       # 640 .+ (1:192)
    h192_2 = view(cache, 833:1024)      # 832 .+ (1:192)
    h384_1 = view(cache, 1025:1408)     # 1024 .+ (1:384)
    h384_2 = view(cache, 1409:1792)     # 1408 .+ (1:384)
    h384_3 = view(cache, 1793:2176)     # 1792 .+ (1:384)
    h768 = view(cache, 2177:2944)       # 2176 .+ (1:768)
    h1152_1 = view(cache, 2945:4096)    # 2944 .+ (1:1152)
    return h64_1, h64_2, h64_3, h64_4, h64_5, h64_6, h128_1, h128_2, h192_1, h192_2, h384_1, h384_2, h384_3, h768, h1152_1
end
incircleadapt_cache(_, cache) = cache
incircleadapt_cache(::Type{T}, ::Nothing) where {T} = incircleadapt_cache(T)
function incircleadapt_cache(::Type{T}) where {T}
    cache = Vec{T}(undef, 2432)         # cache_size = (64) + (64) + (1152 + 1152)
    h48_1 = view(cache, 1:48)           # 0 .+ (1:48)
    h64_1 = view(cache, 65:128)         # 64 .+ (1:64)
    h1152_1 = view(cache, 129:1280)     # 128 .+ (1:1152)
    h1152_2 = view(cache, 1281:2432)    # 1280 .+ (1:1152)
    return h48_1, h64_1, h1152_1, h1152_2
end
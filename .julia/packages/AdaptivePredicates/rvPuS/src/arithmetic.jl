@inline setindex!!(tup::Tuple, value, index) = @inbounds index > length(tup) ? tup : Base.setindex(tup, value, index)
@inline setindex!!(vec::AbstractVector, value, index) = @inbounds Base.setindex!(vec, value, index)
@inline safe_getindex(e, eindex, elen) = (eindex ≤ elen && eindex ≤ length(e)) ? @inbounds(e[eindex]) : zero(eltype(e)) # Shewchuk's code is relying on undefined behaviour from out-of-bounds access where we call this. We need to be careful.

function grow_expansion(elen, e, b, h)
    @inbounds begin
        Q = b
        for eindex in 1:elen
            enow = e[eindex]
            Q, hh = Two_Sum(Q, enow)
            h = setindex!!(h, hh, eindex)
        end
        h = setindex!!(h, Q, elen + 1)
        return h, elen + 1
    end
end

function grow_expansion_zeroelim(elen, e, b, h)
    @inbounds begin
        hindex = 1
        Q = b
        for eindex in 1:elen
            enow = e[eindex]
            Q, hh = Two_Sum(Q, enow)
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
        end
        if !iszero(Q) || isone(hindex)
            h = setindex!!(h, Q, hindex)
            hindex += 1
        end
        return h, hindex - 1
    end
end

function expansion_sum(elen, e, flen, f, h)
    @inbounds begin
        Q = f[1]
        for hindex in 1:elen
            hnow = e[hindex]
            Q, hh = Two_Sum(Q, hnow)
            h = setindex!!(h, hh, hindex)
        end
        hindex = elen + 1
        h = setindex!!(h, Q, hindex)
        hlast = hindex
        for findex in 2:flen
            Q = f[findex]
            for hindex in findex:hlast
                hnow = h[hindex]
                Q, hh = Two_Sum(Q, hnow)
                h = setindex!!(h, hh, hindex)
            end
            hlast += 1
            h = setindex!!(h, Q, hlast)
        end
        return h, hlast
    end
end

function expansion_sum_zeroelim1(elen, e, flen, f, h)
    @inbounds begin
        Q = f[1]
        for hindex in 1:elen
            hnow = e[hindex]
            Q, hh = Two_Sum(Q, hnow)
            h = setindex!!(h, hh, hindex)
        end
        hindex = elen + 1
        h = setindex!!(h, Q, hindex)
        hlast = hindex
        for findex in 2:flen
            Q = f[findex]
            for hindex in findex:hlast
                hnow = h[hindex]
                Q, hh = Two_Sum(Q, hnow)
                h = setindex!!(h, hh, hindex)
            end
            hlast += 1
            h = setindex!!(h, Q, hlast)
        end
        hindex = 0
        for index in 1:hlast
            hnow = h[index]
            if !iszero(hnow)
                hindex += 1
                h = setindex!!(h, hnow, hindex)
            end
        end
        if iszero(hindex)
            return h, 1
        else
            return h, hindex
        end
    end
end

function expansion_sum_zeroelim2(elen, e, flen, f, h)
    @inbounds begin
        hindex = 1
        Q = f[1]
        for eindex in 1:elen
            enow = e[eindex]
            Q, hh = Two_Sum(Q, enow)
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
        end
        h = setindex!!(h, Q, hindex)
        hlast = hindex
        for findex in 2:flen
            hindex = 1
            Q = f[findex]
            for eindex in 1:hlast
                enow = h[eindex]
                Q, hh = Two_Sum(Q, enow)
                if !iszero(hh)
                    h = setindex!!(h, hh, hindex)
                    hindex += 1
                end
            end
            h = setindex!!(h, Q, hindex)
            hlast = hindex
        end
        return h, hlast
    end
end

function fast_expansion_sum(elen, e, flen, f, h)
    @inbounds begin
        enow = e[1]
        fnow = f[1]
        eindex = findex = 1
        if (fnow > enow) == (fnow > -enow)
            Q = enow
            eindex += 1
            enow = safe_getindex(e, eindex, elen)
        else
            Q = fnow
            findex += 1
            fnow = safe_getindex(f, findex, flen)
        end
        hindex = 1
        if (eindex ≤ elen) && (findex ≤ flen)
            if (fnow > enow) == (fnow > -enow)
                Q, hh = Fast_Two_Sum(enow, Q)
                h = setindex!!(h, hh, hindex)
                eindex += 1
                enow = safe_getindex(e, eindex, elen)
            else
                Q, hh = Fast_Two_Sum(fnow, Q)
                h = setindex!!(h, hh, hindex)
                findex += 1
                fnow = safe_getindex(f, findex, flen)
            end
            hindex += 1
            while (eindex ≤ elen) && (findex ≤ flen)
                if (fnow > enow) == (fnow > -enow)
                    Q, hh = Two_Sum(Q, enow)
                    h = setindex!!(h, hh, hindex)
                    eindex += 1
                    enow = safe_getindex(e, eindex, elen)
                else
                    Q, hh = Two_Sum(Q, fnow)
                    h = setindex!!(h, hh, hindex)
                    findex += 1
                    fnow = safe_getindex(f, findex, flen)
                end
                hindex += 1
            end
        end
        while eindex ≤ elen
            Q, hh = Two_Sum(Q, enow)
            h = setindex!!(h, hh, hindex)
            eindex += 1
            enow = safe_getindex(e, eindex, elen)
            hindex += 1
        end
        while findex ≤ flen
            Q, hh = Two_Sum(Q, fnow)
            h = setindex!!(h, hh, hindex)
            findex += 1
            fnow = safe_getindex(f, findex, flen)
            hindex += 1
        end
        h = setindex!!(h, Q, hindex)
        return h, hindex
    end
end

function fast_expansion_sum_zeroelim(elen, e, flen, f, h)
    @inbounds begin
        enow = e[1]
        fnow = f[1]
        eindex = findex = 1
        if (fnow > enow) == (fnow > -enow)
            Q = enow 
            eindex += 1 
            enow = safe_getindex(e, eindex, elen)
        else
            Q = fnow
            findex += 1
            fnow = safe_getindex(f, findex, flen)
        end
        hindex = 1
        if (eindex ≤ elen) && (findex ≤ flen)
            if (fnow > enow) == (fnow > -enow)
                Q, hh = Fast_Two_Sum(enow, Q)
                eindex += 1
                enow = safe_getindex(e, eindex, elen)
            else
                Q, hh = Fast_Two_Sum(fnow, Q)
                findex += 1
                fnow = safe_getindex(f, findex, flen)
            end
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
            while (eindex ≤ elen) && (findex ≤ flen)
                if (fnow > enow) == (fnow > -enow)
                    Q, hh = Two_Sum(Q, enow)
                    eindex += 1
                    enow = safe_getindex(e, eindex, elen)
                else
                    Q, hh = Two_Sum(Q, fnow)
                    findex += 1
                    fnow = safe_getindex(f, findex, flen)
                end
                if !iszero(hh)
                    h = setindex!!(h, hh, hindex)
                    hindex += 1
                end
            end
        end
        while eindex ≤ elen
            Q, hh = Two_Sum(Q, enow)
            eindex += 1
            enow = safe_getindex(e, eindex, elen)
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
        end
        while findex ≤ flen
            Q, hh = Two_Sum(Q, fnow)
            findex += 1
            fnow = safe_getindex(f, findex, flen)
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
        end
        if !iszero(Q) || isone(hindex)
            h = setindex!!(h, Q, hindex)
            hindex += 1
        end
        return h, hindex - 1
    end
end

function linear_expansion_sum(elen, e, flen, f, h)
    @inbounds begin
        enow = e[1]
        fnow = f[1]
        eindex = findex = 1
        if (fnow > enow) == (fnow > -enow)
            g0 = enow
            eindex += 1
            enow = safe_getindex(e, eindex, elen)
        else
            g0 = fnow
            findex += 1
            fnow = safe_getindex(f, findex, flen)
        end
        if (eindex ≤ elen) && ((findex > flen) || ((fnow > enow) == (fnow > -enow)))
            Q, q = Fast_Two_Sum(enow, g0)
            eindex += 1
            enow = safe_getindex(e, eindex, elen)
        else
            Q, q = Fast_Two_Sum(fnow, g0)
            findex += 1
            fnow = safe_getindex(f, findex, flen)
        end
        local hindex = 0
        for outer hindex in 1:(elen+flen-2)
            if (eindex ≤ elen) && ((findex > flen) || ((fnow > enow) == (fnow > -enow)))
                R, hh = Fast_Two_Sum(enow, q)
                h = setindex!!(h, hh, hindex)
                eindex += 1
                enow = safe_getindex(e, eindex, elen)
            else
                R, hh = Fast_Two_Sum(fnow, q)
                h = setindex!!(h, hh, hindex)
                findex += 1
                fnow = safe_getindex(f, findex, flen)
            end
            Q, q = Two_Sum(Q, R)
        end
        h = setindex!!(h, q, hindex + 1)
        h = setindex!!(h, Q, hindex + 2)
        return h, hindex + 2
    end
end

function linear_expansion_sum_zeroelim(elen, e, flen, f, h)
    @inbounds begin
        enow = e[1]
        fnow = f[1]
        eindex = findex = hindex = 1
        if (fnow > enow) == (fnow > -enow)
            g0 = enow
            eindex += 1
            enow = safe_getindex(e, eindex, elen)
        else
            g0 = fnow
            findex += 1
            fnow = safe_getindex(f, findex, flen)
        end
        if (eindex ≤ elen) && ((findex > flen) || ((fnow > enow) == (fnow > -enow)))
            Q, q = Fast_Two_Sum(enow, g0)
            eindex += 1
            enow = safe_getindex(e, eindex, elen)
        else
            Q, q = Fast_Two_Sum(fnow, g0)
            findex += 1
            fnow = safe_getindex(f, findex, flen)
        end
        for _ in 3:(elen+flen)
            if (eindex ≤ elen) && ((findex > flen) || ((fnow > enow) == (fnow > -enow)))
                R, hh = Fast_Two_Sum(enow, q)
                eindex += 1
                enow = safe_getindex(e, eindex, elen)
            else
                R, hh = Fast_Two_Sum(fnow, q)
                findex += 1
                fnow = safe_getindex(f, findex, flen)
            end
            Q, q = Two_Sum(Q, R)
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
        end
        if !iszero(q)
            h = setindex!!(h, q, hindex)
            hindex += 1
        end
        if !iszero(Q) || isone(hindex)
            h = setindex!!(h, Q, hindex)
            hindex += 1
        end
        return h, hindex - 1
    end
end

function scale_expansion(elen, e, b, h)
    @inbounds begin
        bhi, blo = Split(b)
        Q, hh = Two_Product_Presplit(e[1], b, bhi, blo)
        h = setindex!!(h, hh, 1)
        hindex = 2
        for eindex in 2:elen
            enow = e[eindex]
            product1, product0 = Two_Product_Presplit(enow, b, bhi, blo)
            sum, hh = Two_Sum(Q, product0)
            h = setindex!!(h, hh, hindex)
            hindex += 1
            Q, hh = Two_Sum(product1, sum)
            h = setindex!!(h, hh, hindex)
            hindex += 1
        end
        h = setindex!!(h, Q, hindex)
        return h, 2elen
    end
end

function scale_expansion_zeroelim(elen, e, b, h)
    @inbounds begin
        bhi, blo = Split(b)
        Q, hh = Two_Product_Presplit(e[1], b, bhi, blo)
        hindex = 1
        if !iszero(hh)
            h = setindex!!(h, hh, hindex)
            hindex += 1
        end
        for eindex in 2:elen
            enow = e[eindex]
            product1, product0 = Two_Product_Presplit(enow, b, bhi, blo)
            sum, hh = Two_Sum(Q, product0)
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
            Q, hh = Fast_Two_Sum(product1, sum)
            if !iszero(hh)
                h = setindex!!(h, hh, hindex)
                hindex += 1
            end
        end
        if !iszero(Q) || isone(hindex)
            h = setindex!!(h, Q, hindex)
            hindex += 1
        end
        return h, hindex - 1
    end
end

function compress(elen, e, h)
    @inbounds begin
        bottom = elen
        Q = e[bottom]
        for eindex in (elen-1):-1:1
            enow = e[eindex]
            Qnew, q = Fast_Two_Sum(Q, enow)
            if !iszero(q)
                h = setindex!!(h, Qnew, bottom)
                bottom -= 1
                Q = q
            else
                Q = Qnew
            end
        end
        top = 1
        for hindex in (bottom+1):elen
            hnow = h[hindex]
            Qnew, q = Fast_Two_Sum(hnow, Q)
            if !iszero(q)
                h = setindex!!(h, q, top)
                top += 1
            end
            Q = Qnew
        end
        h = setindex!!(h, Q, top)
        return h, top
    end
end

function estimate(elen, e)
    @inbounds begin
        Q = e[1]
        for i in 2:elen
            Q += e[i]
        end
        return Q
    end
end
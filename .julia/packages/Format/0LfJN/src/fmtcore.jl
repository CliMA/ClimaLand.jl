# core formatting functions

### auxiliary functions

### print char n times

function _repprint(out::IO, c::AbstractChar, n::Int)
    wid = textwidth(c)
    n < wid && return
    while n > 0
        print(out, c)
        n -= wid
    end
end

### print string or char

function _truncstr(s::AbstractString, slen, prec)
    prec == 0 && return ("", 0)
    i, n = 0, 1
    siz = ncodeunits(s)
    while n <= siz
	(prec -= textwidth(s[n])) < 0 && break
	i = n
        n = nextind(s, i)
    end
    str = SubString(s, 1, i)
    return (str, textwidth(str))
end

_truncstr(s::AbstractChar, slen, prec) = ("", 0)

function _pfmt_s(out::IO, fs::FormatSpec, s::Union{AbstractString,AbstractChar})
    slen = textwidth(s)
    str, slen = 0 <= fs.prec < slen ? _truncstr(s, slen, fs.prec) : (s, slen)
    prepad = postpad = 0
    pad = fs.width - slen
    if pad > 0
        if fs.align == '<'
            postpad = pad
        elseif fs.align == '^'
            prepad, postpad = pad>>1, (pad+1)>>1
        else
            prepad = pad
        end
    end
    # left padding
    prepad == 0 || _repprint(out, fs.fill, prepad)
    # print string
    print(out, str)
    # right padding
    postpad == 0 || _repprint(out, fs.fill, postpad)
end

_unsigned_abs(x::Signed) = unsigned(abs(x))
_unsigned_abs(x::Bool) = UInt(x)
_unsigned_abs(x::Unsigned) = x
# Special case because unsigned fails for BigInt
_unsigned_abs(x::BigInt) = abs(x)

_ndigits(x, ::_Dec) = ndigits(x)
_ndigits(x, ::_Bin) = ndigits(x, base=2)
_ndigits(x, ::_Oct) = ndigits(x, base=8)
_ndigits(x, ::Union{_Hex, _HEX}) = ndigits(x, base=16)

_sepcnt(::_Dec) = 3
_sepcnt(::Any) = 4

_mul(x::Integer, ::_Dec) = x * 10
_mul(x::Integer, ::_Bin) = x << 1
_mul(x::Integer, ::_Oct) = x << 3
_mul(x::Integer, ::Union{_Hex, _HEX}) = x << 4

_div(x::Integer, ::_Dec) = div(x, 10)
_div(x::Integer, ::_Bin) = x >> 1
_div(x::Integer, ::_Oct) = x >> 3
_div(x::Integer, ::Union{_Hex, _HEX}) = x >> 4

_ipre(op) = ""
_ipre(::_Oct) = "0o"
_ipre(::_Bin) = "0b"
_ipre(::_Hex) = "0x"
_ipre(::_HEX) = "0X"

_digitchar(x::Integer, ::_Bin) = x == 0 ? '0' : '1'
_digitchar(x::Integer, ::_Dec) = Char(Int('0') + x)
_digitchar(x::Integer, ::_Oct) = Char(Int('0') + x)
_digitchar(x::Integer, ::_Hex) = Char(x < 10 ? '0' + x : 'a' + (x - 10))
_digitchar(x::Integer, ::_HEX) = Char(x < 10 ? '0' + x : 'A' + (x - 10))

_signchar(x::Real, s::AbstractChar) = signbit(x) ? '-' : s == '+' ? '+' : s == ' ' ? ' ' : '\0'

### output integers (with or without separators)

function _outint(out::IO, ax::T, op::Op=_Dec()) where {Op, T<:Integer}
    b_lb = _div(ax, op)
    b = one(T)
    while b <= b_lb
        b = _mul(b, op)
    end
    r = ax
    while b > 0
        (q, r) = divrem(r, b)
        print(out, _digitchar(q, op))
        b = _div(b, op)
    end
end

function _outint(out::IO, ax::T, op::Op, numini, sep) where {Op, T<:Integer}
    b_lb = _div(ax, op)
    b = one(T)
    while b <= b_lb
        b = _mul(b, op)
    end
    r = ax
    while true
        (q, r) = divrem(r, b)
        print(out, _digitchar(q, op))
        b = _div(b, op)
        b == 0 && break
        if numini == 0
            numini = _sepcnt(op)
            print(out, sep)
        end
        numini -= 1
    end
end

# Print integer

function _pfmt_i(out::IO, fs::FormatSpec, x::Integer, op::Op) where {Op}
    # calculate actual length
    ax = _unsigned_abs(x)
    xlen = _ndigits(ax, op)
    numsep, numini = fs.tsep ? divrem(xlen-1, _sepcnt(op)) : (0, 0)
    xlen += numsep

    # sign char
    sch = _signchar(x, fs.sign)
    xlen += (sch != '\0')
    # prefix (e.g. 0x, 0b, 0o)
    ip = ""
    if fs.ipre
        ip = _ipre(op)
        xlen += length(ip)
    end
    # printing
    pad = fs.width - xlen
    prepad = postpad = zpad = 0
    if pad > 0
        if fs.zpad
            zpad = pad
        elseif fs.align == '<'
            postpad = pad
        elseif fs.align == '^'
            prepad, postpad = pad>>1, (pad+1)>>1
        else
            prepad = pad
        end
    end
    # left padding
    prepad == 0 || _repprint(out, fs.fill, prepad)
    # print sign
    sch != '\0' && print(out, sch)
    # print prefix
    !isempty(ip) && print(out, ip)
    # print padding zeros
    zpad > 0 && _repprint(out, '0', zpad)
    # print actual digits
    ax == 0 ? print(out, '0') :
        numsep == 0 ? _outint(out, ax, op) : _outint(out, ax, op, numini, fs.sep)
    # right padding
    postpad == 0 || _repprint(out, fs.fill, postpad)
end

function _truncval(v)
    try
        return trunc(Integer, v)
    catch e
        e isa InexactError || rethrow(e)
    end
    try
        return trunc(Int128, v)
    catch e
        e isa InexactError || rethrow(e)
    end
    trunc(BigInt, v)
end

### print floating point numbers

function _pfmt_f(out::IO, fs::FormatSpec, x::AbstractFloat)
    # Handle %
    percentflag = (fs.typ == '%')
    percentflag && (x *= 100)
    # separate sign, integer, and decimal part
    prec = fs.prec
    rax = round(abs(x); digits = prec)
    sch = _signchar(x, fs.sign)
    intv = _truncval(rax)
    decv = rax - intv

    # calculate length
    xlen = ndigits(intv)
    numsep, numini = fs.tsep ? divrem(xlen - 1, 3) : (0, 0)
    xlen += ifelse(prec > 0, prec + 1, 0) + (sch != '\0') + numsep + percentflag

    # calculate padding needed
    pad = fs.width - xlen
    prepad = postpad = zpad = 0
    if pad > 0
        if fs.zpad
            zpad = pad
        elseif fs.align == '<'
            postpad = pad
        elseif fs.align == '^'
            prepad, postpad = pad>>1, (pad+1)>>1
        else
            prepad = pad
        end
    end
    # left padding
    prepad == 0 || _repprint(out, fs.fill, prepad)

    # print sign
    sch != '\0' && print(out, sch)

    # print padding zeros
    zpad > 0 && _repprint(out, '0', zpad)

    # print integer part
    intv == 0 ? print(out, '0') :
        numsep == 0 ? _outint(out, intv) : _outint(out, intv, _Dec(), numini, fs.sep)

    # print decimal part
    if prec > 0
        print(out, '.')
        idecv = round(Integer, decv * exp10(prec))
        nd = ndigits(idecv)
        nd < prec && _repprint(out, '0', prec - nd)
        _outint(out, idecv)
    end

    # print trailing percent sign
    percentflag && print(out, '%')

    # right padding
    postpad == 0 || _repprint(out, fs.fill, postpad)
end

function _pfmt_e(out::IO, fs::FormatSpec, x::AbstractFloat)
    # extract sign, significand, and exponent
    prec = fs.prec
    ax = abs(x)
    sch = _signchar(x, fs.sign)
    if ax == 0.0
        e = 0
        u = zero(x)
    else
        rax = round(ax; sigdigits = prec + 1)
        e = floor(Integer, log10(rax))  # exponent
        u = round(rax * exp10(-e); sigdigits = prec + 1)  # significand
        i = 0
        v10 = 1
        while isinf(u)
            i += 1
            i > 18 && (u = 0.0; e = 0; break)
            v10 *= 10
            u = round(v10 * rax * exp10(-e - i); sigdigits = prec + 1)
        end
    end

    # calculate length
    xlen = 6 + prec + (sch != '\0') + (abs(e) > 99 ? ndigits(abs(e)) - 2 : 0)

    # calculate padding
    pad = fs.width - xlen
    prepad = postpad = zpad = 0
    if pad > 0
        if fs.zpad
            zpad = pad
        elseif fs.align == '<'
            postpad = pad
        elseif fs.align == '^'
            prepad, postpad = pad>>1, (pad+1)>>1
        else
            prepad = pad
        end
    end

    # left padding
    prepad == 0 || _repprint(out, fs.fill, prepad)

    # print sign
    sch != '\0' && print(out, sch)

    # print padding zeros
    zpad > 0 && _repprint(out, '0', zpad)

    # print actual digits
    intv = trunc(Integer, u)
    decv = u - intv
    if intv == 0 && decv != 0
        intv = 1
        decv -= 1
    end

    # print integer part (should always be 0-9)
    print(out, Char(Int('0') + intv))

    # print decimal part
    if prec > 0
        print(out, '.')
        idecv = round(Integer, decv * exp10(prec))
        nd = ndigits(idecv)
        nd < prec && _repprint(out, '0', prec - nd)
        _outint(out, idecv)
    end

    # print exponent
    print(out, isuppercase(fs.typ) ? 'E' : 'e')
    if e >= 0
        print(out, '+')
    else
        print(out, '-')
        e = -e
    end
    e < 10 && print(out, '0')
    _outint(out, e)
    # right padding
    postpad == 0 || _repprint(out, fs.fill, postpad)
end

function _pfmt_g(out::IO, fs::FormatSpec, x::AbstractFloat)
    # number decomposition
    ax = abs(x)
    if 1.0e-4 <= ax < 1.0e6
        _pfmt_f(out, fs, x)
    else
        _pfmt_e(out, fs, x)
    end
end

function _pfmt_specialf(out::IO, fs::FormatSpec, x::AbstractFloat)
    if isinf(x)
        _pfmt_s(out, fs, x > 0 ? "Inf" : "-Inf")
    else
        @assert isnan(x)
        _pfmt_s(out, fs, "NaN")
    end
end

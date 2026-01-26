# test format spec parsing

# default spec
@testset "Default spec" begin
    fs = FormatSpec("")
    @test fs.typ == 's'
    @test fs.fill == ' '
    @test fs.align == '<'
    @test fs.sign == '-'
    @test fs.width == -1
    @test fs.prec == -1
    @test fs.ipre == false
    @test fs.zpad == false
    @test fs.tsep == false
end

_contains(s, r) = occursin(r, s)

@testset "Show" begin
    x = FormatSpec("#8,d")
    io = IOBuffer()
    show(io, x)
    str = String(take!(io))
    @test _contains(str, "width = 8")
end

@testset "Literal incorrect" begin
    @test_throws ErrorException FormatSpec("Z")
end

# more cases

@testset "FormatSpec(\"d\")" begin
    fs = FormatSpec("d")
    @test fs == FormatSpec('d')
    @test fs.align == '>'

    @test FormatSpec("8x") == FormatSpec('x'; width=8)
    @test FormatSpec("08b") == FormatSpec('b'; width=8, zpad=true)
    @test FormatSpec("12f") == FormatSpec('f'; width=12, prec=6)
    @test FormatSpec("12.7f") == FormatSpec('f'; width=12, prec=7)
    @test FormatSpec("+08o") == FormatSpec('o'; width=8, zpad=true, sign='+')

    @test FormatSpec("8") == FormatSpec('s'; width=8)
    @test FormatSpec(".6f") == FormatSpec('f'; prec=6)
    @test FormatSpec("<8d") == FormatSpec('d'; width=8, align='<')
    @test FormatSpec("#<8d") == FormatSpec('d'; width=8, fill='#', align='<')
    @test FormatSpec("⋆<8d") == FormatSpec('d'; width=8, fill='⋆', align='<')
    @test FormatSpec("#^8d") == FormatSpec('d'; width=8, fill='#', align='^')
    @test FormatSpec("#8,d") == FormatSpec('d'; width=8, ipre=true, tsep=true)
end

@testset "Format prefix" begin
    @test pyfmt("#b", 6) == "0b110"
    @test pyfmt("#o", 6) == "0o6"
    @test pyfmt("#x", 6) == "0x6"
    @test pyfmt("#X", 10) == "0XA"
end

@testset "Format string" begin
    @test pyfmt("", "abc") == "abc"
    @test pyfmt("", "αβγ") == "αβγ"
    @test pyfmt("s", "abc") == "abc"
    @test pyfmt("s", "αβγ") == "αβγ"
    @test pyfmt("2s", "abc") == "abc"
    @test pyfmt("2s", "αβγ") == "αβγ"
    @test pyfmt("5s", "abc") == "abc  "
    @test pyfmt("5s", "αβγ") == "αβγ  "
    @test pyfmt(">5s", "abc") == "  abc"
    @test pyfmt(">5s", "αβγ") == "  αβγ"
    @test pyfmt("^5s", "abc") == " abc "
    @test pyfmt("^5s", "αβγ") == " αβγ "
    @test pyfmt("*>5s", "abc") == "**abc"
    @test pyfmt("⋆>5s", "αβγ") == "⋆⋆αβγ"
    @test pyfmt("*^5s", "abc") == "*abc*"
    @test pyfmt("⋆^5s", "αβγ") == "⋆αβγ⋆"
    @test pyfmt("*<5s", "abc") == "abc**"
    @test pyfmt("⋆<5s", "αβγ") == "αβγ⋆⋆"

    # Issue #83 (in Formatting.jl)
    @test pyfmt(".3s", "\U1f595\U1f596") == "🖕"
    @test pyfmt(".4s", "\U1f595\U1f596") == "🖕🖖"
    @test pyfmt(".5s", "\U1f595\U1f596") == "🖕🖖"

    @test pyfmt("6.3s", "\U1f595\U1f596") == "🖕    "
    @test pyfmt("6.4s", "\U1f595\U1f596") == "🖕🖖  "
    @test pyfmt("6.5s", "\U1f595\U1f596") == "🖕🖖  "
end

@testset "Format Symbol" begin
    @test pyfmt("", :abc) == "abc"
    @test pyfmt("s", :abc) == "abc"
end

@testset "Format Char" begin
    @test pyfmt("", 'c') == "c"
    @test pyfmt("", 'γ') == "γ"
    @test pyfmt("c", 'c') == "c"
    @test pyfmt("c", 'γ') == "γ"
    @test pyfmt("3c", 'c') == "c  "
    @test pyfmt("3c", 'γ') == "γ  "
    @test pyfmt(">3c", 'c') == "  c"
    @test pyfmt(">3c", 'γ') == "  γ"
    @test pyfmt("^3c", 'c') == " c "
    @test pyfmt("^3c", 'γ') == " γ "
    @test pyfmt("*>3c", 'c') == "**c"
    @test pyfmt("⋆>3c", 'γ') == "⋆⋆γ"
    @test pyfmt("*^3c", 'c') == "*c*"
    @test pyfmt("⋆^3c", 'γ') == "⋆γ⋆"
    @test pyfmt("*<3c", 'c') == "c**"
    @test pyfmt("⋆<3c", 'γ') == "γ⋆⋆"
end

@testset "Format integer" begin

    @test pyfmt("", 1234) == "1234"
    @test pyfmt("d", 1234) == "1234"
    @test pyfmt("n", 1234) == "1234"
    @test pyfmt("x", 0x2ab) == "2ab"
    @test pyfmt("X", 0x2ab) == "2AB"
    @test pyfmt("o", 0o123) == "123"
    @test pyfmt("b", 0b1101) == "1101"

    @test pyfmt("d", 0) == "0"
    @test pyfmt("d", 9) == "9"
    @test pyfmt("d", 10) == "10"
    @test pyfmt("d", 99) == "99"
    @test pyfmt("d", 100) == "100"
    @test pyfmt("d", 1000) == "1000"

    @test pyfmt("06d", 123) == "000123"
    @test pyfmt("+6d", 123) == "  +123"
    @test pyfmt("+06d", 123) == "+00123"
    @test pyfmt(" d", 123) == " 123"
    @test pyfmt(" 6d", 123) == "   123"
    @test pyfmt("<6d", 123) == "123   "
    @test pyfmt(">6d", 123) == "   123"
    @test pyfmt("^6d", 123) == " 123  "
    @test pyfmt("*<6d", 123) == "123***"
    @test pyfmt("⋆<6d", 123) == "123⋆⋆⋆"
    @test pyfmt("*^6d", 123) == "*123**"
    @test pyfmt("⋆^6d", 123) == "⋆123⋆⋆"
    @test pyfmt("*>6d", 123) == "***123"
    @test pyfmt("⋆>6d", 123) == "⋆⋆⋆123"
    @test pyfmt("< 6d", 123) == " 123  "
    @test pyfmt("<+6d", 123) == "+123  "
    @test pyfmt("^ 6d", 123) == "  123 "
    @test pyfmt("^+6d", 123) == " +123 "
    @test pyfmt("> 6d", 123) == "   123"
    @test pyfmt(">+6d", 123) == "  +123"

    @test pyfmt("+d", -123) == "-123"
    @test pyfmt("-d", -123) == "-123"
    @test pyfmt(" d", -123) == "-123"
    @test pyfmt("06d", -123) == "-00123"
    @test pyfmt("<6d", -123) == "-123  "
    @test pyfmt("^6d", -123) == " -123 "
    @test pyfmt(">6d", -123) == "  -123"

    # Issue #110 (in Formatting.jl)
    f = FormatExpr("{:+d}")
    for T in (Int8, Int16, Int32, Int64, Int128)
        @test format(f, typemin(T)) == string(typemin(T))
    end
end

@testset "Python separators (, and _)" begin

    # May need to change so that 'd' is a default, if none specified for integer value
    @test pyfmt(",d", 1) == "1"
    @test pyfmt(",d", 12) == "12"
    @test pyfmt(",d", 123) == "123"
    @test pyfmt(",d", 1234) == "1,234"
    @test pyfmt(",d", 12345) == "12,345"
    @test pyfmt(",d", 123456) == "123,456"
    @test pyfmt(",d", 1234567) == "1,234,567"
    @test pyfmt(",d", 12345678) == "12,345,678"
    @test pyfmt(",d", 123456789) == "123,456,789"
    @test pyfmt(",d", 1234567890) == "1,234,567,890"

    @test pyfmt("_d", 1) == "1"
    @test pyfmt("_d", 12) == "12"
    @test pyfmt("_d", 123) == "123"
    @test pyfmt("_d", 1234) == "1_234"
    @test pyfmt("_d", 12345) == "12_345"
    @test pyfmt("_d", 123456) == "123_456"
    @test pyfmt("_d", 1234567) == "1_234_567"
    @test pyfmt("_d", 12345678) == "12_345_678"
    @test pyfmt("_d", 123456789) == "123_456_789"
    @test pyfmt("_d", 1234567890) == "1_234_567_890"

    @test pyfmt("_x", 0x1) == "1"
    @test pyfmt("_x", 0x12) == "12"
    @test pyfmt("_x", 0x123) == "123"
    @test pyfmt("_x", 0x1234) == "1234"
    @test pyfmt("_x", 0x12345) == "1_2345"
    @test pyfmt("_x", 0x123456) == "12_3456"
    @test pyfmt("_x", 0x1234567) == "123_4567"
    @test pyfmt("_x", 0x12345678) == "1234_5678"
    @test pyfmt("_x", 0x123456789) == "1_2345_6789"
    @test pyfmt("_x", 0x123456789a) == "12_3456_789a"

    @test pyfmt("_b", 0b1) == "1"
    @test pyfmt("_b", 0b10) == "10"
    @test pyfmt("_b", 0b101) == "101"
    @test pyfmt("_b", 0b1010) == "1010"
    @test pyfmt("_b", 0b10101) == "1_0101"
    @test pyfmt("_b", 0b101010) == "10_1010"

    @test pyfmt("_o", 0o1) == "1"
    @test pyfmt("_o", 0o12) == "12"
    @test pyfmt("_o", 0o123) == "123"
    @test pyfmt("_o", 0o1234) == "1234"
    @test pyfmt("_o", 0o12345) == "1_2345"
    @test pyfmt("_o", 0o123456) == "12_3456"

    @test pyfmt("15,d",  123456789) == "    123,456,789"
    @test pyfmt("15,d", -123456789) == "   -123,456,789"

    @test pyfmt("06,d",  1234) == "01,234"
    @test pyfmt("07,d",  1234) == "001,234"

    # May need to change how separators are handled with zero padding,
    # to produce the same results as Python

#    @test pyfmt("08,d",  1234) == "0,001,234"
#    @test pyfmt("09,d",  1234) == "0,001,234"
#    @test pyfmt("010,d", 1234) == "00,001,234"

    @test pyfmt("+#012_b", 42) == "+0b0010_1010"
#    @test pyfmt("+#013_b", 42) == "+0b0_0010_1010"
#    @test pyfmt("+#014_b", 42) == "+0b0_0010_1010"

end

@testset "Format floating point (f)" begin

    @test pyfmt("", 0.125) == "0.125"
    @test pyfmt("f", 0.0) == "0.000000"
    @test pyfmt("f", -0.0) == "-0.000000"
    @test pyfmt("f", 0.001) == "0.001000"
    @test pyfmt("f", 0.125) == "0.125000"
    @test pyfmt("f", 1.0/3) == "0.333333"
    @test pyfmt("f", 1.0/6) == "0.166667"
    @test pyfmt("f", -0.125) == "-0.125000"
    @test pyfmt("f", -1.0/3) == "-0.333333"
    @test pyfmt("f", -1.0/6) == "-0.166667"
    @test pyfmt("f", 1234.5678) == "1234.567800"
    @test pyfmt("8f", 1234.5678) == "1234.567800"

    @test pyfmt("8.2f", 8.376) == "    8.38"
    @test pyfmt("<8.2f", 8.376) == "8.38    "
    @test pyfmt("^8.2f", 8.376) == "  8.38  "
    @test pyfmt(">8.2f", 8.376) == "    8.38"
    @test pyfmt("8.2f", -8.376) == "   -8.38"
    @test pyfmt("<8.2f", -8.376) == "-8.38   "
    @test pyfmt("^8.2f", -8.376) == " -8.38  "
    @test pyfmt(">8.2f", -8.376) == "   -8.38"
    @test pyfmt(".0f", 8.376) == "8"

    # Note: these act differently than Python, but it looks like Python might be wrong
    # in at least some of these cases (zero padding should be *after* the sign, IMO)
    @test pyfmt("<08.2f", 8.376) == "00008.38"
    @test pyfmt("^08.2f", 8.376) == "00008.38"
    @test pyfmt(">08.2f", 8.376) == "00008.38"
    @test pyfmt("<08.2f", -8.376) == "-0008.38"
    @test pyfmt("^08.2f", -8.376) == "-0008.38"
    @test pyfmt(">08.2f", -8.376) == "-0008.38"
    @test pyfmt("*<8.2f", 8.376) == "8.38****"
    @test pyfmt("⋆<8.2f", 8.376) == "8.38⋆⋆⋆⋆"
    @test pyfmt("*^8.2f", 8.376) == "**8.38**"
    @test pyfmt("⋆^8.2f", 8.376) == "⋆⋆8.38⋆⋆"
    @test pyfmt("*>8.2f", 8.376) == "****8.38"
    @test pyfmt("⋆>8.2f", 8.376) == "⋆⋆⋆⋆8.38"
    @test pyfmt("*<8.2f", -8.376) == "-8.38***"
    @test pyfmt("⋆<8.2f", -8.376) == "-8.38⋆⋆⋆"
    @test pyfmt("*^8.2f", -8.376) == "*-8.38**"
    @test pyfmt("⋆^8.2f", -8.376) == "⋆-8.38⋆⋆"
    @test pyfmt("*>8.2f", -8.376) == "***-8.38"
    @test pyfmt("⋆>8.2f", -8.376) == "⋆⋆⋆-8.38"

    @test pyfmt(".2f", 0.999) == "1.00"
    @test pyfmt(".2f", 0.996) == "1.00"
    @test pyfmt("6.2f", 9.999) == " 10.00"
    # Floating point error can upset this one (i.e. 0.99500000 or 0.994999999)
    @test (pyfmt(".2f", 0.995) == "1.00" || pyfmt(".2f", 0.995) == "0.99")
    @test pyfmt(".2f", 0.994) == "0.99"

    # issue #102 (from Formatting.jl)
    @test pyfmt("15.6f", 1e20) == "100000000000000000000.000000"
end

@testset "Format floating point (e)" begin

    @test pyfmt("E", 0.0) == "0.000000E+00"
    @test pyfmt("e", 0.0) == "0.000000e+00"
    @test pyfmt("e", 0.001) == "1.000000e-03"
    @test pyfmt("e", 0.125) == "1.250000e-01"
    @test pyfmt("e", 100/3) == "3.333333e+01"
    @test pyfmt("e", -0.125) == "-1.250000e-01"
    @test pyfmt("e", -100/6) == "-1.666667e+01"
    @test pyfmt("e", 1234.5678) == "1.234568e+03"
    @test pyfmt("8e", 1234.5678) == "1.234568e+03"

    @test pyfmt("<12.2e", 13.89) == "1.39e+01    "
    @test pyfmt("^12.2e", 13.89) == "  1.39e+01  "
    @test pyfmt(">12.2e", 13.89) == "    1.39e+01"
    @test pyfmt("*<12.2e", 13.89) == "1.39e+01****"
    @test pyfmt("⋆<12.2e", 13.89) == "1.39e+01⋆⋆⋆⋆"
    @test pyfmt("*^12.2e", 13.89) == "**1.39e+01**"
    @test pyfmt("⋆^12.2e", 13.89) == "⋆⋆1.39e+01⋆⋆"
    @test pyfmt("*>12.2e", 13.89) == "****1.39e+01"
    @test pyfmt("⋆>12.2e", 13.89) == "⋆⋆⋆⋆1.39e+01"
    @test pyfmt("012.2e", 13.89) == "00001.39e+01"
    @test pyfmt("012.2e", -13.89) == "-0001.39e+01"
    @test pyfmt("+012.2e", 13.89) == "+0001.39e+01"

    @test pyfmt(".1e", 0.999) == "1.0e+00"
    @test pyfmt(".1e", 0.996) == "1.0e+00"
    # Floating point error can upset this one (i.e. 0.99500000 or 0.994999999)
    @test (pyfmt(".1e", 0.995) == "1.0e+00" || pyfmt(".1e", 0.995) == "9.9e-01")
    @test pyfmt(".1e", 0.994) == "9.9e-01"
    @test pyfmt(".1e", 0.6) == "6.0e-01"
    @test pyfmt(".1e", 0.9) == "9.0e-01"

    # issue #61 (from Formatting.jl)
    @test pyfmt("1.0e", 1e-21) == "1e-21"
    @test pyfmt("1.1e", 1e-21) == "1.0e-21"

    @test pyfmt("10.2e", 1.2e100) == " 1.20e+100"
    @test pyfmt("11.2e", BigFloat("1.2e1000")) == " 1.20e+1000"
    @test pyfmt("11.2e", BigFloat("1.2e-1000")) == " 1.20e-1000"
    @test pyfmt("9.2e", 9.999e9) == " 1.00e+10"
    @test pyfmt("10.2e", 9.999e99) == " 1.00e+100"
    @test pyfmt("11.2e", BigFloat("9.999e999")) == " 1.00e+1000"
    @test pyfmt("10.2e", -9.999e-100) == " -1.00e-99"

    # issue #84 (from Formatting.jl)
    @test pyfmt("+11.3e", 1.0e-309) == "+1.000e-309"
    @test pyfmt("+11.3e", 1.0e-313) == "+1.000e-313"

    # issue #108 (from Formatting.jl)
    @test pyfmt(".1e", 0.0003) == "3.0e-04"
    @test pyfmt(".1e", 0.0006) == "6.0e-04"
end

@testset "Format percentage (%)" begin
    @test pyfmt("8.2%", 0.123456)  == "  12.35%"
    @test pyfmt("<8.2%", 0.123456) == "12.35%  "
    @test pyfmt("^8.2%", 0.123456) == " 12.35% "
    @test pyfmt(">8.2%", 0.123456) == "  12.35%"
end

@testset "Format special floating point value" begin

    @test pyfmt("f", NaN) == "NaN"
    @test pyfmt("e", NaN) == "NaN"
    @test pyfmt("f", NaN32) == "NaN"
    @test pyfmt("e", NaN32) == "NaN"

    @test pyfmt("f", Inf) == "Inf"
    @test pyfmt("e", Inf) == "Inf"
    @test pyfmt("f", Inf32) == "Inf"
    @test pyfmt("e", Inf32) == "Inf"

    @test pyfmt("f", -Inf) == "-Inf"
    @test pyfmt("e", -Inf) == "-Inf"
    @test pyfmt("f", -Inf32) == "-Inf"
    @test pyfmt("e", -Inf32) == "-Inf"

    @test pyfmt("<5f", Inf) == "Inf  "
    @test pyfmt("^5f", Inf) == " Inf "
    @test pyfmt(">5f", Inf) == "  Inf"

    @test pyfmt("*<5f", Inf) == "Inf**"
    @test pyfmt("⋆<5f", Inf) == "Inf⋆⋆"
    @test pyfmt("*^5f", Inf) == "*Inf*"
    @test pyfmt("⋆^5f", Inf) == "⋆Inf⋆"
    @test pyfmt("*>5f", Inf) == "**Inf"
    @test pyfmt("⋆>5f", Inf) == "⋆⋆Inf"
end

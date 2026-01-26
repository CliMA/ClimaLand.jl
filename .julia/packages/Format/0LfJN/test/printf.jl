# Based on stdlib Printf.jl runtests.jl

using Format
using Test

@testset "cfmt" begin

@testset "%p" begin

    # pointers
    if Sys.WORD_SIZE == 64
        @test cfmt("%20p", 0) == "  0x0000000000000000"
        @test cfmt("%-20p", 0) == "0x0000000000000000  "
        @test cfmt("%20p", C_NULL) == "  0x0000000000000000"
        @test cfmt("%-20p", C_NULL) == "0x0000000000000000  "
    elseif Sys.WORD_SIZE == 32
        @test cfmt("%20p", 0) == "          0x00000000"
        @test cfmt("%-20p", 0) == "0x00000000          "
        @test cfmt("%20p", C_NULL) == "          0x00000000"
        @test cfmt("%-20p", C_NULL) == "0x00000000          "
    end

    # Printf.jl #40318
    @test cfmt("%p", 0xfffffffffffe0000) == "0xfffffffffffe0000"

end

@testset "%a" begin

    # hex float
    @test cfmt("%a", 0.0) == "0x0p+0"
    @test cfmt("%a", -0.0) == "-0x0p+0"
    @test cfmt("%.3a", 0.0) == "0x0.000p+0"
    @test cfmt("%a", 1.5) == "0x1.8p+0"
    @test cfmt("%a", 1.5f0) == "0x1.8p+0"
    @test cfmt("%a", big"1.5") == "0x1.8p+0"
    @test cfmt("%#.0a", 1.5) == "0x2.p+0"
    @test cfmt("%+30a", 1/3) == "         +0x1.5555555555555p-2"

    @test cfmt("%a", 1.5) == "0x1.8p+0"
    @test cfmt("%a", 3.14) == "0x1.91eb851eb851fp+1"
    @test cfmt("%.0a", 3.14) == "0x2p+1"
    @test cfmt("%.1a", 3.14) == "0x1.9p+1"
    @test cfmt("%.2a", 3.14) == "0x1.92p+1"
    @test cfmt("%#a", 3.14) == "0x1.91eb851eb851fp+1"
    @test cfmt("%#.0a", 3.14) == "0x2.p+1"
    @test cfmt("%#.1a", 3.14) == "0x1.9p+1"
    @test cfmt("%#.2a", 3.14) == "0x1.92p+1"
    @test cfmt("%.6a", 1.5) == "0x1.800000p+0"

end

@testset "%g" begin

    # %g
    for (val, res) in ((12345678., "1.23457e+07"),
                    (1234567.8, "1.23457e+06"),
                    (123456.78, "123457"),
                    (12345.678, "12345.7"),
                    (12340000.0, "1.234e+07"))
        @test cfmt("%.6g", val) == res
    end
    for (val, res) in ((big"12345678.", "1.23457e+07"),
                    (big"1234567.8", "1.23457e+06"),
                    (big"123456.78", "123457"),
                    (big"12345.678", "12345.7"))
        @test cfmt("%.6g", val) == res
    end
    for (fmt, val) in (("%10.5g", "     123.4"),
                    ("%+10.5g", "    +123.4"),
                    ("% 10.5g", "     123.4"),
                    ("%#10.5g", "    123.40"),
                    ("%-10.5g", "123.4     "),
                    ("%-+10.5g", "+123.4    "),
                    ("%010.5g", "00000123.4")),
        num in (123.4, big"123.4")
        @test cfmt(Format.FmtSpec(fmt), num) == val
    end
    @test cfmt( "%10.5g", -123.4 ) == "    -123.4"
    @test cfmt( "%010.5g", -123.4 ) == "-0000123.4"
    @test cfmt( "%.6g", 12340000.0 ) == "1.234e+07"
    @test cfmt( "%#.6g", 12340000.0 ) == "1.23400e+07"
    @test cfmt( "%10.5g", big"-123.4" ) == "    -123.4"
    @test cfmt( "%010.5g", big"-123.4" ) == "-0000123.4"
    @test cfmt( "%.6g", big"12340000.0" ) == "1.234e+07"
    @test cfmt( "%#.6g", big"12340000.0") == "1.23400e+07"

    # %g regression gh #14331
    @test cfmt( "%.5g", 42) == "42"
    @test cfmt( "%#.2g", 42) == "42."
    @test cfmt( "%#.5g", 42) == "42.000"

    @test cfmt("%g", 0.00012) == "0.00012"
    @test cfmt("%g", 0.000012) == "1.2e-05"
    @test cfmt("%g", 123456.7) == "123457"
    @test cfmt("%g", 1234567.8) == "1.23457e+06"

    # %g regression gh #41631
    for (val, res) in ((Inf, "Inf"),
                       (-Inf, "-Inf"),
                       (NaN, "NaN"),
                       (-NaN, "NaN"))
        @test cfmt("%g", val) == res
        @test cfmt("%G", val) == res
    end

    # zeros
    @test cfmt("%.15g", 0) == "0"
    @test cfmt("%#.15g", 0) == "0.00000000000000"

end

@testset "%f" begin

    # Inf / NaN handling
    @test cfmt("%f", Inf) == "Inf"
    @test cfmt("%+f", Inf) == "+Inf"
    @test cfmt("% f", Inf) == " Inf"
    @test cfmt("% #f", Inf) == " Inf"
    @test cfmt("%07f", Inf) == "    Inf"
    @test cfmt("%f", -Inf) == "-Inf"
    @test cfmt("%+f", -Inf) == "-Inf"
    @test cfmt("%07f", -Inf) == "   -Inf"
    @test cfmt("%f", NaN) == "NaN"
    @test cfmt("%+f", NaN) == "+NaN"
    @test cfmt("% f", NaN) == " NaN"
    @test cfmt("% #f", NaN) == " NaN"
    @test cfmt("%07f", NaN) == "    NaN"
    @test cfmt("%e", big"Inf") == "Inf"
    @test cfmt("%e", big"NaN") == "NaN"

    @test cfmt("%.0f", 3e142) == "29999999999999997463140672961703247153805615792184250659629251954072073858354858644285983761764971823910371920726635399393477049701891710124032"

    @test cfmt("%f", 1.234) == "1.234000"
    @test cfmt("%F", 1.234) == "1.234000"
    @test cfmt("%+f", 1.234) == "+1.234000"
    @test cfmt("% f", 1.234) == " 1.234000"
    @test cfmt("%f", -1.234) == "-1.234000"
    @test cfmt("%+f", -1.234) == "-1.234000"
    @test cfmt("% f", -1.234) == "-1.234000"
    @test cfmt("%#f", 1.234) == "1.234000"
    @test cfmt("%.2f", 1.234) == "1.23"
    @test cfmt("%.2f", 1.235) == "1.24"
    @test cfmt("%.2f", 0.235) == "0.23"
    @test cfmt("%4.1f", 1.234) == " 1.2"
    @test cfmt("%8.1f", 1.234) == "     1.2"
    @test cfmt("%+8.1f", 1.234) == "    +1.2"
    @test cfmt("% 8.1f", 1.234) == "     1.2"
    @test cfmt("% 7.1f", 1.234) == "    1.2"
    @test cfmt("% 08.1f", 1.234) == " 00001.2"
    @test cfmt("%08.1f", 1.234) == "000001.2"
    @test cfmt("%-08.1f", 1.234) == "1.2     "
    @test cfmt("%-8.1f", 1.234) == "1.2     "
    @test cfmt("%08.1f", -1.234) == "-00001.2"
    @test cfmt("%09.1f", -1.234) == "-000001.2"
    @test cfmt("%09.1f", 1.234) == "0000001.2"
    @test cfmt("%+09.1f", 1.234) == "+000001.2"
    @test cfmt("% 09.1f", 1.234) == " 000001.2"
    @test cfmt("%+ 09.1f", 1.234) == "+000001.2"
    @test cfmt("%+ 09.1f", 1.234) == "+000001.2"
    @test cfmt("%+ 09.0f", 1.234) == "+00000001"
    @test cfmt("%+ #09.0f", 1.234) == "+0000001."

    # Printf.jl issue #40303
    @test cfmt("%+7.1f", 9.96) == "  +10.0"
    @test cfmt("% 7.1f", 9.96) == "   10.0"
end

@testset "%e" begin

    # Inf / NaN handling
    @test cfmt("%e", Inf) == "Inf"
    @test cfmt("%+e", Inf) == "+Inf"
    @test cfmt("% e", Inf) == " Inf"
    @test cfmt("% #e", Inf) == " Inf"
    @test cfmt("%07e", Inf) == "    Inf"
    @test cfmt("%e", -Inf) == "-Inf"
    @test cfmt("%+e", -Inf) == "-Inf"
    @test cfmt("%07e", -Inf) == "   -Inf"
    @test cfmt("%e", NaN) == "NaN"
    @test cfmt("%+e", NaN) == "+NaN"
    @test cfmt("% e", NaN) == " NaN"
    @test cfmt("% #e", NaN) == " NaN"
    @test cfmt("%07e", NaN) == "    NaN"
    @test cfmt("%e", big"Inf") == "Inf"
    @test cfmt("%e", big"NaN") == "NaN"

    # scientific notation
    @test cfmt("%.0e", 3e142) == "3e+142"
    @test cfmt("%#.0e", 3e142) == "3.e+142"
    @test cfmt("%.0e", big"3e142") == "3e+142"
    @test cfmt("%#.0e", big"3e142") == "3.e+142"

    @test cfmt("%.0e", big"3e1042") == "3e+1042"

    @test cfmt("%e", 3e42) == "3.000000e+42"
    @test cfmt("%E", 3e42) == "3.000000E+42"
    @test cfmt("%e", 3e-42) == "3.000000e-42"
    @test cfmt("%E", 3e-42) == "3.000000E-42"

    @test cfmt("%e", 1.234) == "1.234000e+00"
    @test cfmt("%E", 1.234) == "1.234000E+00"
    @test cfmt("%+e", 1.234) == "+1.234000e+00"
    @test cfmt("% e", 1.234) == " 1.234000e+00"
    @test cfmt("%e", -1.234) == "-1.234000e+00"
    @test cfmt("%+e", -1.234) == "-1.234000e+00"
    @test cfmt("% e", -1.234) == "-1.234000e+00"
    @test cfmt("%#e", 1.234) == "1.234000e+00"
    @test cfmt("%.2e", 1.234) == "1.23e+00"
    @test cfmt("%.2e", 1.235) == "1.24e+00"
    @test cfmt("%.2e", 0.235) == "2.35e-01"
    @test cfmt("%4.1e", 1.234) == "1.2e+00"
    @test cfmt("%8.1e", 1.234) == " 1.2e+00"
    @test cfmt("%+8.1e", 1.234) == "+1.2e+00"
    @test cfmt("% 8.1e", 1.234) == " 1.2e+00"
    @test cfmt("% 7.1e", 1.234) == " 1.2e+00"
    @test cfmt("% 08.1e", 1.234) == " 1.2e+00"
    @test cfmt("%08.1e", 1.234) == "01.2e+00"
    @test cfmt("%-08.1e", 1.234) == "1.2e+00 "
    @test cfmt("%-8.1e", 1.234) == "1.2e+00 "
    @test cfmt("%-8.1e", 1.234) == "1.2e+00 "
    @test cfmt("%08.1e", -1.234) == "-1.2e+00"
    @test cfmt("%09.1e", -1.234) == "-01.2e+00"
    @test cfmt("%09.1e", 1.234) == "001.2e+00"
    @test cfmt("%+09.1e", 1.234) == "+01.2e+00"
    @test cfmt("% 09.1e", 1.234) == " 01.2e+00"
    @test cfmt("%+ 09.1e", 1.234) == "+01.2e+00"
    @test cfmt("%+ 09.1e", 1.234) == "+01.2e+00"
    @test cfmt("%+ 09.0e", 1.234) == "+0001e+00"
    @test cfmt("%+ #09.0e", 1.234) == "+001.e+00"

    # Printf.jl #40303
    @test cfmt("%+9.1e", 9.96) == " +1.0e+01"
    @test cfmt("% 9.1e", 9.96) == "  1.0e+01"
end

@testset "strings" begin

    @test cfmt("%.1s", "foo") == "f"
    @test cfmt("%s", "%%%%") == "%%%%"
    @test cfmt("%s", "Hallo heimur") == "Hallo heimur"
    @test cfmt("%+s", "Hallo heimur") == "Hallo heimur"
    @test cfmt("% s", "Hallo heimur") == "Hallo heimur"
    @test cfmt("%+ s", "Hallo heimur") == "Hallo heimur"
    @test cfmt("%1s", "Hallo heimur") == "Hallo heimur"
    @test cfmt("%20s", "Hallo") == "               Hallo"
    @test cfmt("%-20s", "Hallo") == "Hallo               "
    @test cfmt("%0-20s", "Hallo") == "Hallo               "
    @test cfmt("%.20s", "Hallo heimur") == "Hallo heimur"
    @test cfmt("%20.5s", "Hallo heimur") == "               Hallo"
    @test cfmt("%.0s", "Hallo heimur") == ""
    @test cfmt("%20.0s", "Hallo heimur") == "                    "
    @test cfmt("%.s", "Hallo heimur") == ""
    @test cfmt("%20.s", "Hallo heimur") == "                    "
    @test cfmt("%s", "test") == "test"
    @test cfmt("%s", "tést") == "tést"

    @test cfmt("%8s", "test") == "    test"
    @test cfmt("%-8s", "test") == "test    "

    @test cfmt("%s", :test) == "test"
    @test cfmt("%#s", :test) == ":test"
    @test cfmt("%#8s", :test) == "   :test"
    @test cfmt("%#-8s", :test) == ":test   "

    @test cfmt("%8.3s", "test") == "     tes"
    @test cfmt("%#8.3s", "test") == "     \"te"
    @test cfmt("%-8.3s", "test") == "tes     "
    @test cfmt("%#-8.3s", "test") == "\"te     "
    @test cfmt("%.3s", "test") == "tes"
    @test cfmt("%#.3s", "test") == "\"te"
    @test cfmt("%-.3s", "test") == "tes"
    @test cfmt("%#-.3s", "test") == "\"te"

    # Printf.jl issue #41068
    @test cfmt("%.2s", "föó") == "fö"
    @test cfmt("%5s", "föó") == "  föó"
    @test cfmt("%6s", "😍🍕") == "  😍🍕"
    @test cfmt("%2c", '🍕') == "🍕"
    @test cfmt("%3c", '🍕') == " 🍕"

end

@testset "chars" begin

    @test cfmt("%c", 'a') == "a"
    @test cfmt("%c",  32) == " "
    @test cfmt("%c",  36) == "\$"
    @test cfmt("%3c", 'a') == "  a"
    @test cfmt( "%c", 'x') == "x"
    @test cfmt("%+c", 'x') == "x"
    @test cfmt("% c", 'x') == "x"
    @test cfmt("%+ c", 'x') == "x"
    @test cfmt("%1c", 'x') == "x"
    @test cfmt("%20c", 'x') == "                   x"
    @test cfmt("%-20c", 'x') == "x                   "
    @test cfmt("%-020c", 'x') == "x                   "
    @test cfmt("%c", 65) == "A"
    @test cfmt("%c", 'A') == "A"
    @test cfmt("%3c", 'A') == "  A"
    @test cfmt("%-3c", 'A') == "A  "
    @test cfmt("%c", 248) == "ø"
    @test cfmt("%c", 'ø') == "ø"
    @test cfmt("%c", "ø") == "ø"
    @test cfmt("%c", '𐀀') == "𐀀"

end

function _test_flags(val, vflag::AbstractString, fmt::AbstractString, res::AbstractString, prefix::AbstractString)
    vflag = string("%", vflag)
    space_fmt = string(length(res) + length(prefix) + 3, fmt)
    fsign = string((val < 0 ? "-" : "+"), prefix)
    nsign = string((val < 0 ? "-" : " "), prefix)
    osign = val < 0 ? string("-", prefix) : string(prefix, "0")
    esign = string(val < 0 ? "-" : "", prefix)
    esignend = val < 0 ? "" : " "

    for (flag::AbstractString, ans::AbstractString) in (
            ("", string("  ", nsign, res)),
            ("+", string("  ", fsign, res)),
            (" ", string("  ", nsign, res)),
            ("0", string(osign, "00", res)),
            ("-", string(esign, res, "  ", esignend)),
            ("0+", string(fsign, "00", res)),
            ("0 ", string(nsign, "00", res)),
            ("-+", string(fsign, res, "  ")),
            ("- ", string(nsign, res, "  ")),
        )
        fmt_string = string(vflag, flag, space_fmt)
        @test cfmt(fmt_string, val) == ans
    end
end

@testset "basics" begin

    @test cfmt("%10.5d", 4) == "     00004"
    @test cfmt("%d", typemax(Int64)) == "9223372036854775807"

    for (fmt, val) in (("%7.2f", "   1.23"),
                   ("%-7.2f", "1.23   "),
                   ("%07.2f", "0001.23"),
                   ("%.0f", "1"),
                   ("%#.0f", "1."),
                   ("%.4e", "1.2345e+00"),
                   ("%.4E", "1.2345E+00"),
                   ("%.2a", "0x1.3cp+0"),
                   ("%.2A", "0X1.3CP+0")),
        num in (1.2345, big"1.2345")
        @test cfmt(Format.FmtSpec(fmt), num) == val
    end

    for (fmt, val) in (("%i", "42"),
                   ("%u", "42"),
                   #("Test: %i", "Test: 42"),
                   ("%#x", "0x2a"),
                   ("%x", "2a"),
                   ("%X", "2A"),
                   ("% i", " 42"),
                   ("%+i", "+42"),
                   ("%4i", "  42"),
                   ("%-4i", "42  "),
                   ("%f", "42.000000"),
                   ("%g", "42"),
                   ("%e", "4.200000e+01")),
        num in (UInt16(42), UInt32(42), UInt64(42), UInt128(42),
                Int16(42), Int32(42), Int64(42), Int128(42)) #, big"42") #Fix!
        @test cfmt(Format.FmtSpec(fmt), num) == val
    end

    for i in (
            (42, "", "i", "42", ""),
            (42, "", "d", "42", ""),

            (42, "", "u", "42", ""),
            (42, "", "x", "2a", ""),
            (42, "", "X", "2A", ""),
            (42, "", "o", "52", ""),

            (42, "#", "x", "2a", "0x"),
            (42, "#", "X", "2A", "0X"),
            (42, "#", "o", "052", ""),

            (1.2345, "", ".2f", "1.23", ""),
            (1.2345, "", ".2e", "1.23e+00", ""),
            (1.2345, "", ".2E", "1.23E+00", ""),

            (1.2345, "#", ".0f", "1.", ""),
            (1.2345, "#", ".0e", "1.e+00", ""),
            (1.2345, "#", ".0E", "1.E+00", ""),

            (1.2345, "", ".2a", "1.3cp+0", "0x"),
            (1.2345, "", ".2A", "1.3CP+0", "0X"),
        )
        _test_flags(i...)
        _test_flags(-i[1], i[2:5]...)
    end

    # Check bug with trailing nul printing BigFloat
    @test cfmt("%.330f", BigFloat(1))[end] != '\0'

    # Check bugs with truncated output printing BigFloat
    @test cfmt("%f", parse(BigFloat, "1e400")) ==
           "10000000000000000000000000000000000000000000000000000000000000000000000000000025262527574416492004687051900140830217136998040684679611623086405387447100385714565637522507383770691831689647535911648520404034824470543643098638520633064715221151920028135130764414460468236314621044034960475540018328999334468948008954289495190631358190153259681118693204411689043999084305348398480210026863210192871358464.000000"

    # Check that Printf does not attempt to output more than 8KiB worth of digits
    #@test_throws ArgumentError cfmt("%f", parse(BigFloat, "1e99999"))

    # Check bug with precision > length of string
    @test cfmt("%4.2s", "a") == "   a"

    # Printf.jl issue #29662
    @test cfmt("%12.3e", pi*1e100) == "  3.142e+100"

    #@test string(Printf.Format("%a").formats[1]) == "%a"
    #@test string(Printf.Format("%a").formats[1]; modifier="R") == "%Ra"

    @test cfmt("%d", 3.14) == "3"
    @test cfmt("%2d", 3.14) == " 3"
    @test_broken cfmt("%2d", big(3.14)) == " 3"
    @test cfmt("%s", 1) == "1"
    @test cfmt("%f", 1) == "1.000000"
    @test cfmt("%e", 1) == "1.000000e+00"
    @test cfmt("%g", 1) == "1"

    # Printf.jl issue #39748
    @test cfmt("%.16g", 194.4778127560983) == "194.4778127560983"
    @test cfmt("%.17g", 194.4778127560983) == "194.4778127560983"
    @test cfmt("%.18g", 194.4778127560983) == "194.477812756098302"
    @test cfmt("%.1g", 1.7976931348623157e308) == "2e+308"
    @test cfmt("%.2g", 1.7976931348623157e308) == "1.8e+308"
    @test cfmt("%.3g", 1.7976931348623157e308) == "1.8e+308"

    # print float as %d uses round(x)
    @test cfmt("%d", 25.5) == "26"

    # Printf.jl issue #37539
    #@test @sprintf(" %.1e\n", 0.999) == " 1.0e+00\n"
    #@test @sprintf("   %.1f", 9.999) == "   10.0"

    # Printf.jl issue #37552
    @test cfmt("%d", 1.0e100) == "10000000000000000159028911097599180468360808563945281389781327557747838772170381060813469985856815104"
    @test_broken cfmt("%d", 3//1) == "3"
    @test cfmt("%d", Inf) == "Inf"
end
    
@testset "integers" begin

    @test cfmt("% d",  42) == " 42"
    @test cfmt("% d", -42) == "-42"
    @test cfmt("% 5d",  42) == "   42"
    @test cfmt("% 5d", -42) == "  -42"
    @test cfmt("% 15d",  42) == "             42"
    @test cfmt("% 15d", -42) == "            -42"
    @test cfmt("%+d",  42) == "+42"
    @test cfmt("%+d", -42) == "-42"
    @test cfmt("%+5d",  42) == "  +42"
    @test cfmt("%+5d", -42) == "  -42"
    @test cfmt("%+15d",  42) == "            +42"
    @test cfmt("%+15d", -42) == "            -42"
    @test cfmt("%0d",  42) == "42"
    @test cfmt("%0d", -42) == "-42"
    @test cfmt("%05d",  42) == "00042"
    @test cfmt("%05d", -42) == "-0042"
    @test cfmt("%015d",  42) == "000000000000042"
    @test cfmt("%015d", -42) == "-00000000000042"
    @test cfmt("%-d",  42) == "42"
    @test cfmt("%-d", -42) == "-42"
    @test cfmt("%-5d",  42) == "42   "
    @test cfmt("%-5d", -42) == "-42  "
    @test cfmt("%-15d",  42) == "42             "
    @test cfmt("%-15d", -42) == "-42            "
    @test cfmt("%-0d",  42) == "42"
    @test cfmt("%-0d", -42) == "-42"
    @test cfmt("%-05d",  42) == "42   "
    @test cfmt("%-05d", -42) == "-42  "
    @test cfmt("%-015d",  42) == "42             "
    @test cfmt("%-015d", -42) == "-42            "
    @test cfmt("%0-d",  42) == "42"
    @test cfmt("%0-d", -42) == "-42"
    @test cfmt("%0-5d",  42) == "42   "
    @test cfmt("%0-5d", -42) == "-42  "
    @test cfmt("%0-15d",  42) == "42             "
    @test cfmt("%0-15d", -42) == "-42            "

    @test cfmt("%lld", 18446744065119617025) == "18446744065119617025"
    @test cfmt("%+8lld", 100) == "    +100"
    @test cfmt("%+.8lld", 100) == "+00000100"
    @test cfmt("%+10.8lld", 100) == " +00000100"
    #@test_throws Printf.InvalidFormatStringError Printf.Format("%_1lld")
    @test cfmt("%-1.5lld", -100) == "-00100"
    @test cfmt("%5lld", 100) == "  100"
    @test cfmt("%5lld", -100) == " -100"
    @test cfmt("%-5lld", 100) == "100  "
    @test cfmt("%-5lld", -100) == "-100 "
    @test cfmt("%-.5lld", 100) == "00100"
    @test cfmt("%-.5lld", -100) == "-00100"
    @test cfmt("%-8.5lld", 100) == "00100   "
    @test cfmt("%-8.5lld", -100) == "-00100  "
    @test cfmt("%05lld", 100) == "00100"
    @test cfmt("%05lld", -100) == "-0100"
    @test cfmt("% lld", 100) == " 100"
    @test cfmt("% lld", -100) == "-100"
    @test cfmt("% 5lld", 100) == "  100"
    @test cfmt("% 5lld", -100) == " -100"
    @test cfmt("% .5lld", 100) == " 00100"
    @test cfmt("% .5lld", -100) == "-00100"
    @test cfmt("% 8.5lld", 100) == "   00100"
    @test cfmt("% 8.5lld", -100) == "  -00100"
    @test cfmt("%.0lld", 0) == "0"
    @test cfmt("%#+21.18llx", -100) == "-0x000000000000000064"
    @test cfmt("%#.25llo", -100) == "-00000000000000000000000144"
    @test cfmt("%#+24.20llo", -100) == "  -000000000000000000144"
    @test cfmt("%#+18.21llX", -100) == "-0X000000000000000000064"
    @test cfmt("%#+20.24llo", -100) == "-0000000000000000000000144"
    @test cfmt("%#+25.22llu", -1) == "  -0000000000000000000001"
    @test cfmt("%#+25.22llu", -1) == "  -0000000000000000000001"
    @test cfmt("%#+30.25llu", -1) == "    -0000000000000000000000001"
    @test cfmt("%+#25.22lld", -1) == "  -0000000000000000000001"
    @test cfmt("%#-8.5llo", 100) == "000144  "
    @test cfmt("%#-+ 08.5lld", 100) == "+00100  "
    @test cfmt("%#-+ 08.5lld", 100) == "+00100  "
    @test cfmt("%.40lld",  1) == "0000000000000000000000000000000000000001"
    @test cfmt("% .40lld",  1) == " 0000000000000000000000000000000000000001"
    @test cfmt("% .40d",  1) == " 0000000000000000000000000000000000000001"
    @test cfmt("%lld",  18446744065119617025) == "18446744065119617025"

    @test cfmt("%#012x",  1) == "0x0000000001"
    @test cfmt("%#04.8x",  1) == "0x00000001"

    @test cfmt("%#-08.2x",  1) == "0x01    "
    @test cfmt("%#08o",  1) == "00000001"
    @test cfmt("%d",  1024) == "1024"
    @test cfmt("%d", -1024) == "-1024"
    @test cfmt("%i",  1024) == "1024"
    @test cfmt("%i", -1024) == "-1024"
    @test cfmt("%u",  1024) == "1024"
    @test cfmt("%u",  UInt(4294966272)) == "4294966272"
    @test cfmt("%o",  511) == "777"
    @test cfmt("%o",  UInt(4294966785)) == "37777777001"
    @test cfmt("%x",  305441741) == "1234abcd"
    @test cfmt("%x",  UInt(3989525555)) == "edcb5433"
    @test cfmt("%X",  305441741) == "1234ABCD"
    @test cfmt("%X",  UInt(3989525555)) == "EDCB5433"
    @test cfmt("%+d",  1024) == "+1024"
    @test cfmt("%+d", -1024) == "-1024"
    @test cfmt("%+i",  1024) == "+1024"
    @test cfmt("%+i", -1024) == "-1024"
    @test cfmt("%+u",  1024) == "+1024"
    @test cfmt("%+u",  UInt(4294966272)) == "+4294966272"
    @test cfmt("%+o",  511) == "+777"
    @test cfmt("%+o",  UInt(4294966785)) == "+37777777001"
    @test cfmt("%+x",  305441741) == "+1234abcd"
    @test cfmt("%+x",  UInt(3989525555)) == "+edcb5433"
    @test cfmt("%+X",  305441741) == "+1234ABCD"
    @test cfmt("%+X",  UInt(3989525555)) == "+EDCB5433"
    @test cfmt("% d",  1024) == " 1024"
    @test cfmt("% d", -1024) == "-1024"
    @test cfmt("% i",  1024) == " 1024"
    @test cfmt("% i", -1024) == "-1024"
    @test cfmt("% u",  1024) == " 1024"
    @test cfmt("% u",  UInt(4294966272)) == " 4294966272"
    @test cfmt("% o",  511) == " 777"
    @test cfmt("% o",  UInt(4294966785)) == " 37777777001"
    @test cfmt("% x",  305441741) == " 1234abcd"
    @test cfmt("% x",  UInt(3989525555)) == " edcb5433"
    @test cfmt("% X",  305441741) == " 1234ABCD"
    @test cfmt("% X",  UInt(3989525555)) == " EDCB5433"
    @test cfmt("%+ d",  1024) == "+1024"
    @test cfmt("%+ d", -1024) == "-1024"
    @test cfmt("%+ i",  1024) == "+1024"
    @test cfmt("%+ i", -1024) == "-1024"
    @test cfmt("%+ u",  1024) == "+1024"
    @test cfmt("%+ u",  UInt(4294966272)) == "+4294966272"
    @test cfmt("%+ o",  511) == "+777"
    @test cfmt("%+ o",  UInt(4294966785)) == "+37777777001"
    @test cfmt("%+ x",  305441741) == "+1234abcd"
    @test cfmt("%+ x",  UInt(3989525555)) == "+edcb5433"
    @test cfmt("%+ X",  305441741) == "+1234ABCD"
    @test cfmt("%+ X",  UInt(3989525555)) == "+EDCB5433"
    @test cfmt("%#o",  511) == "0777"
    @test cfmt("%#o",  UInt(4294966785)) == "037777777001"
    @test cfmt("%#x",  305441741) == "0x1234abcd"
    @test cfmt("%#x",  UInt(3989525555)) == "0xedcb5433"
    @test cfmt("%#X",  305441741) == "0X1234ABCD"
    @test cfmt("%#X",  UInt(3989525555)) == "0XEDCB5433"
    @test cfmt("%#o",  UInt(0)) == "00"
    @test cfmt("%#x",  UInt(0)) == "0x0"
    @test cfmt("%#X",  UInt(0)) == "0X0"
    @test cfmt("%1d",  1024) == "1024"
    @test cfmt("%1d", -1024) == "-1024"
    @test cfmt("%1i",  1024) == "1024"
    @test cfmt("%1i", -1024) == "-1024"
    @test cfmt("%1u",  1024) == "1024"
    @test cfmt("%1u",  UInt(4294966272)) == "4294966272"
    @test cfmt("%1o",  511) == "777"
    @test cfmt("%1o",  UInt(4294966785)) == "37777777001"
    @test cfmt("%1x",  305441741) == "1234abcd"
    @test cfmt("%1x",  UInt(3989525555)) == "edcb5433"
    @test cfmt("%1X",  305441741) == "1234ABCD"
    @test cfmt("%1X",  UInt(3989525555)) == "EDCB5433"
    @test cfmt("%20d",  1024) == "                1024"
    @test cfmt("%20d", -1024) == "               -1024"
    @test cfmt("%20i",  1024) == "                1024"
    @test cfmt("%20i", -1024) == "               -1024"
    @test cfmt("%20u",  1024) == "                1024"
    @test cfmt("%20u",  UInt(4294966272)) == "          4294966272"
    @test cfmt("%20o",  511) == "                 777"
    @test cfmt("%20o",  UInt(4294966785)) == "         37777777001"
    @test cfmt("%20x",  305441741) == "            1234abcd"
    @test cfmt("%20x",  UInt(3989525555)) == "            edcb5433"
    @test cfmt("%20X",  305441741) == "            1234ABCD"
    @test cfmt("%20X",  UInt(3989525555)) == "            EDCB5433"
    @test cfmt("%-20d",  1024) == "1024                "
    @test cfmt("%-20d", -1024) == "-1024               "
    @test cfmt("%-20i",  1024) == "1024                "
    @test cfmt("%-20i", -1024) == "-1024               "
    @test cfmt("%-20u",  1024) == "1024                "
    @test cfmt("%-20u",  UInt(4294966272)) == "4294966272          "
    @test cfmt("%-20o",  511) == "777                 "
    @test cfmt("%-20o",  UInt(4294966785)) == "37777777001         "
    @test cfmt("%-20x",  305441741) == "1234abcd            "
    @test cfmt("%-20x",  UInt(3989525555)) == "edcb5433            "
    @test cfmt("%-20X",  305441741) == "1234ABCD            "
    @test cfmt("%-20X",  UInt(3989525555)) == "EDCB5433            "
    @test cfmt("%020d",  1024) == "00000000000000001024"
    @test cfmt("%020d", -1024) == "-0000000000000001024"
    @test cfmt("%020i",  1024) == "00000000000000001024"
    @test cfmt("%020i", -1024) == "-0000000000000001024"
    @test cfmt("%020u",  1024) == "00000000000000001024"
    @test cfmt("%020u",  UInt(4294966272)) == "00000000004294966272"
    @test cfmt("%020o",  511) == "00000000000000000777"
    @test cfmt("%020o",  UInt(4294966785)) == "00000000037777777001"
    @test cfmt("%020x",  305441741) == "0000000000001234abcd"
    @test cfmt("%020x",  UInt(3989525555)) == "000000000000edcb5433"
    @test cfmt("%020X",  305441741) == "0000000000001234ABCD"
    @test cfmt("%020X",  UInt(3989525555)) == "000000000000EDCB5433"
    @test cfmt("%#20o",  511) == "                0777"
    @test cfmt("%#20o",  UInt(4294966785)) == "        037777777001"
    @test cfmt("%#20x",  305441741) == "          0x1234abcd"
    @test cfmt("%#20x",  UInt(3989525555)) == "          0xedcb5433"
    @test cfmt("%#20X",  305441741) == "          0X1234ABCD"
    @test cfmt("%#20X",  UInt(3989525555)) == "          0XEDCB5433"
    @test cfmt("%#020o",  511) == "00000000000000000777"
    @test cfmt("%#020o",  UInt(4294966785)) == "00000000037777777001"
    @test cfmt("%#020x",  305441741) == "0x00000000001234abcd"
    @test cfmt("%#020x",  UInt(3989525555)) == "0x0000000000edcb5433"
    @test cfmt("%#020X",  305441741) == "0X00000000001234ABCD"
    @test cfmt("%#020X",  UInt(3989525555)) == "0X0000000000EDCB5433"
    @test cfmt("%0-20d",  1024) == "1024                "
    @test cfmt("%0-20d", -1024) == "-1024               "
    @test cfmt("%0-20i",  1024) == "1024                "
    @test cfmt("%0-20i", -1024) == "-1024               "
    @test cfmt("%0-20u",  1024) == "1024                "
    @test cfmt("%0-20u",  UInt(4294966272)) == "4294966272          "
    @test cfmt("%-020o",  511) == "777                 "
    @test cfmt("%-020o",  UInt(4294966785)) == "37777777001         "
    @test cfmt("%-020x",  305441741) == "1234abcd            "
    @test cfmt("%-020x",  UInt(3989525555)) == "edcb5433            "
    @test cfmt("%-020X",  305441741) == "1234ABCD            "
    @test cfmt("%-020X",  UInt(3989525555)) == "EDCB5433            "
    @test cfmt("%.20d",  1024) == "00000000000000001024"
    @test cfmt("%.20d", -1024) == "-00000000000000001024"
    @test cfmt("%.20i",  1024) == "00000000000000001024"
    @test cfmt("%.20i", -1024) == "-00000000000000001024"
    @test cfmt("%.20u",  1024) == "00000000000000001024"
    @test cfmt("%.20u",  UInt(4294966272)) == "00000000004294966272"
    @test cfmt("%.20o",  511) == "00000000000000000777"
    @test cfmt("%.20o",  UInt(4294966785)) == "00000000037777777001"
    @test cfmt("%.20x",  305441741) == "0000000000001234abcd"
    @test cfmt("%.20x",  UInt(3989525555)) == "000000000000edcb5433"
    @test cfmt("%.20X",  305441741) == "0000000000001234ABCD"
    @test cfmt("%.20X",  UInt(3989525555)) == "000000000000EDCB5433"
    @test cfmt("%20.5d",  1024) == "               01024"
    @test cfmt("%20.5d", -1024) == "              -01024"
    @test cfmt("%20.5i",  1024) == "               01024"
    @test cfmt("%20.5i", -1024) == "              -01024"
    @test cfmt("%20.5u",  1024) == "               01024"
    @test cfmt("%20.5u",  UInt(4294966272)) == "          4294966272"
    @test cfmt("%20.5o",  511) == "               00777"
    @test cfmt("%20.5o",  UInt(4294966785)) == "         37777777001"
    @test cfmt("%20.5x",  305441741) == "            1234abcd"
    @test cfmt("%20.10x",  UInt(3989525555)) == "          00edcb5433"
    @test cfmt("%20.5X",  305441741) == "            1234ABCD"
    @test cfmt("%20.10X",  UInt(3989525555)) == "          00EDCB5433"
    @test cfmt("%020.5d",  1024) == "               01024"
    @test cfmt("%020.5d", -1024) == "              -01024"
    @test cfmt("%020.5i",  1024) == "               01024"
    @test cfmt("%020.5i", -1024) == "              -01024"
    @test cfmt("%020.5u",  1024) == "               01024"
    @test cfmt("%020.5u",  UInt(4294966272)) == "          4294966272"
    @test cfmt("%020.5o",  511) == "               00777"
    @test cfmt("%020.5o",  UInt(4294966785)) == "         37777777001"
    @test cfmt("%020.5x",  305441741) == "            1234abcd"
    @test cfmt("%020.10x",  UInt(3989525555)) == "          00edcb5433"
    @test cfmt("%020.5X",  305441741) == "            1234ABCD"
    @test cfmt("%020.10X",  UInt(3989525555)) == "          00EDCB5433"
    @test cfmt("%20.0d",  1024) == "                1024"
    @test cfmt("%20.d", -1024) == "               -1024"
    @test cfmt("%20.d",  0) == "                   0"
    @test cfmt("%20.0i",  1024) == "                1024"
    @test cfmt("%20.i", -1024) == "               -1024"
    @test cfmt("%20.i",  0) == "                   0"
    @test cfmt("%20.u",  1024) == "                1024"
    @test cfmt("%20.0u",  UInt(4294966272)) == "          4294966272"
    @test cfmt("%20.u",  UInt(0)) == "                   0"
    @test cfmt("%20.o",  511) == "                 777"
    @test cfmt("%20.0o",  UInt(4294966785)) == "         37777777001"
    @test cfmt("%20.o",  UInt(0)) == "                   0"
    @test cfmt("%20.x",  305441741) == "            1234abcd"
    @test cfmt("%20.0x",  UInt(3989525555)) == "            edcb5433"
    @test cfmt("%20.x",  UInt(0)) == "                   0"
    @test cfmt("%20.X",  305441741) == "            1234ABCD"
    @test cfmt("%20.0X",  UInt(3989525555)) == "            EDCB5433"
    @test cfmt("%20.X",  UInt(0)) == "                   0"

    # Printf.jl issue #41971
    @test cfmt("%4d", typemin(Int8)) == "-128"
    @test cfmt("%4d", typemax(Int8)) == " 127"
    @test cfmt("%6d", typemin(Int16)) == "-32768"
    @test cfmt("%6d", typemax(Int16)) == " 32767"
    @test cfmt("%11d", typemin(Int32)) == "-2147483648"
    @test cfmt("%11d", typemax(Int32)) == " 2147483647"
    @test cfmt("%20d", typemin(Int64)) == "-9223372036854775808"
    @test cfmt("%20d", typemax(Int64)) == " 9223372036854775807"
    @test cfmt("%40d", typemin(Int128)) == "-170141183460469231731687303715884105728"
    @test cfmt("%40d", typemax(Int128)) == " 170141183460469231731687303715884105727"
end

# Printf.jl issue #52749
@test cfmt("%.160g", 1.38e-23) == "1.380000000000000060010582465734078799297660966782642624395399644741944111814291318296454846858978271484375e-23"

end # @testset "cfmt"

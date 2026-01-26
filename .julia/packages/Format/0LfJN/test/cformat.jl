_erfinv(z) = sqrt(π) * Base.Math.@horner(z, 0, 1, 0, π/12, 0, 7π^2/480, 0, 127π^3/40320, 0,
                                         4369π^4/5806080, 0, 34807π^5/182476800) / 2

set_seed!(seed) = Random.seed!(seed)

function test_equality()
    println( "test cformat equality...")
    set_seed!(10)
    fmts = [ (x->@sprintf("%10.4f",x), "%10.4f"),
             (x->@sprintf("%f", x),    "%f"),
             (x->@sprintf("%e", x),    "%e"),
             (x->@sprintf("%10f", x),  "%10f"),
             (x->@sprintf("%.3f", x),  "%.3f"),
             (x->@sprintf("%.3e", x),  "%.3e")]
    for (mfmtr,fmt) in fmts
        for i in 1:10000
            n = _erfinv( rand() * 1.99 - 1.99/2.0 )
            expect = mfmtr( n )
            actual = cfmt( fmt, n )
            @test expect == actual
        end
    end

    fmts = [ (x->@sprintf("%d",x),    "%d"),
             (x->@sprintf("%10d",x),  "%10d"),
             (x->@sprintf("%010d",x), "%010d"),
             (x->@sprintf("%-10d",x), "%-10d")]
    for (mfmtr,fmt) in fmts
        for i in 1:10000
            j = round(Int, _erfinv( rand() * 1.99 - 1.99/2.0 ) * 100000 )
            expect = mfmtr( j )
            actual = cfmt( fmt, j )
            @test expect == actual
        end
    end
    println( "...Done" )
end

@time test_equality()

include("speedtest.jl")

@testset "test commas..." begin
    @test cfmt( "%'d", 1000 ) == "1,000"
    @test cfmt( "%'d", -1000 ) == "-1,000"
    @test cfmt( "%'d", 100 ) == "100"
    @test cfmt( "%'d", -100 ) == "-100"
    @test cfmt( "%'f", Inf ) == "Inf"
    @test cfmt( "%'f", -Inf ) == "-Inf"
    @test cfmt( "%'s", 1000.0 ) == "1,000.0"
    @test cfmt( "%'s", 1234567.0 ) == "1.234567e6"
    @test cfmt( "%'g", 1000.0 ) == "1,000"
end

@testset "Test bug introduced by stdlib/Printf rewrite" begin
    @test cfmt( "%4.2s", "a" ) == "   a"
end

@testset "Test precision with wide characters" begin
    @test cfmt( "%10.9s", "\U1f355" ^ 6) == "  " * "\U1f355" ^ 4
    @test cfmt( "%10.10s", "\U1f355" ^ 6) == "\U1f355" ^ 5
end

@testset "test format..." begin
    @test format( 10 ) == "10"
    @test format( 10.0 ) == "10"
    @test format( 10.0, precision=2 ) == "10.00"
    @test format( 111//100, precision=2 ) == "1.11"
    @test format( 111//100 ) == "111/100"
    @test format( 1234, commas=true ) == "1,234"
    @test format( 1234, conversion="f", precision=2 ) == "1234.00"
    @test format( 1.23, precision=3 ) == "1.230"
    @test format( 1.23, precision=3, stripzeros=true ) == "1.23"
    @test format( 1.00, precision=3, stripzeros=true ) == "1"

    @test format( 1.0, conversion="e", stripzeros=true ) == "1e+00"
    @test format( 1.0, conversion="e", precision=4 ) == "1.0000e+00"
    @test format( 1.0, signed=true ) == "+1"
    @test format( 1.0, positivespace=true ) == " 1"
    @test_throws ErrorException format( 1234.56, signed=true, commas=true )
    @test format( 1.0, width=6, precision=4, stripzeros=true, leftjustified=true) == "1     "
end

@testset "hex output" begin
    @test format( 1118, conversion="x" ) == "45e"
    @test format( 1118, width=4, conversion="x" ) == " 45e"
    @test format( 1118, width=4, zeropadding=true, conversion="x" ) == "045e"
    @test format( 1118, alternative=true, conversion="x" ) == "0x45e"
    @test format( 1118, width=4, alternative=true, conversion="x" ) == "0x45e"
    @test format( 1118, width=6, alternative=true, conversion="x", zeropadding=true ) == "0x045e"
end

@testset "mixed fractions" begin
    @test format( 3//2, mixedfraction=true ) == "1_1/2"
    @test format( -3//2, mixedfraction=true ) == "-1_1/2"
    @test format( 3//100, mixedfraction=true ) == "3/100"
    @test format( -3//100, mixedfraction=true ) == "-3/100"
    @test format( 307//100, mixedfraction=true ) == "3_7/100"
    @test format( -307//100, mixedfraction=true ) == "-3_7/100"
    @test format( 307//100, mixedfraction=true, fractionwidth=6 ) == "3_07/100"
    @test format( -307//100, mixedfraction=true, fractionwidth=6 ) == "-3_07/100"
    @test format( -302//100, mixedfraction=true ) == "-3_1/50"
    # try to make the denominator 100
    @test format( -302//100, mixedfraction=true,tryden = 100 ) == "-3_2/100"
    @test format( -302//30, mixedfraction=true,tryden = 100 ) == "-10_1/15" # lose precision otherwise
    @test format( -302//100, mixedfraction=true,tryden = 100,fractionwidth=6 ) == "-3_02/100" # lose precision otherwise
end

@testset "commas" begin
    @test format( 12345678, width=10, commas=true ) == "12,345,678"
    # it would try to squeeze out the commas
    @test format( 12345678, width=9, commas=true ) == "12345,678"
    # until it can't anymore
    @test format( 12345678, width=8, commas=true ) == "12345678"
    @test format( 12345678, width=7, commas=true ) == "12345678"

    # only the numerator would have commas
    @test format( 1111//1000, commas=true ) == "1,111/1000"

    # this shows how, with enough space, parens line up with empty spaces
    @test format(  12345678, width=12, commas=true, parens=true )== " 12,345,678 "
    @test format( -12345678, width=12, commas=true, parens=true )== "(12,345,678)"
    # same with unspecified width
    @test format(  12345678, commas=true, parens=true )== " 12,345,678 "
    @test format( -12345678, commas=true, parens=true )== "(12,345,678)"
end

@testset "autoscale" begin
    @test format( 1.2e9, autoscale = :metric ) == "1.2G"
    @test format( 1.2e6, autoscale = :metric ) == "1.2M"
    @test format( 1.2e3, autoscale = :metric ) == "1.2k"
    @test format( 1.2e-6, autoscale = :metric ) == "1.2μ"
    @test format( 1.2e-9, autoscale = :metric ) == "1.2n"
    @test format( 1.2e-12, autoscale = :metric ) == "1.2p"

    @test format( 1.2e9, autoscale = :finance ) == "1.2b"
    @test format( 1.2e6, autoscale = :finance ) == "1.2m"
    @test format( 1.2e3, autoscale = :finance ) == "1.2k"

    @test format( 0x40000000, autoscale = :binary ) == "1Gi"
    @test format( 0x100000, autoscale = :binary ) == "1Mi"
    @test format( 0x800, autoscale = :binary ) == "2Ki"
    @test format( 0x400, autoscale = :binary ) == "1Ki"
end

@testset "suffix" begin
    @test format( 100.00, precision=2, suffix="%" ) == "100.00%"
    @test format( 100, precision=2, suffix="%" ) == "100%"
    @test format( 100, precision=2, suffix="%", conversion="f" ) == "100.00%"
end


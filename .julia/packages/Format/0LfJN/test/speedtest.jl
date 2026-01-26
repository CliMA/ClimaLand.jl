println( "\nTest speed" )

const BENCH_REP = 200_000

function native_int()
    for i in 1:BENCH_REP
        @sprintf( "%10d", i )
    end
end
@static if VERSION >= v"1.6"
function format_int()
    fmt = Printf.Format( "%10d" )
    for i in 1:BENCH_REP
        Printf.format( fmt, i )
    end
end
end
function runtime_int()
    for i in 1:BENCH_REP
        cfmt( "%10d", i )
    end
end
function runtime_int_spec()
    fs = Format.FmtSpec("%10d")
    for i in 1:BENCH_REP
        cfmt( fs, i )
    end
end
function runtime_int_bypass()
    f = generate_formatter( "%10d" )
    for i in 1:BENCH_REP
        f( i )
    end
end

println( "integer @sprintf speed")
@time native_int()
@static if VERSION >= v"1.6"
    println( "integer format speed")
    @time format_int()
end
println( "integer cfmt speed")
@time runtime_int()
println( "integer cfmt spec speed")
@time runtime_int_spec()
println( "integer cfmt speed, bypass repeated lookup")
@time runtime_int_bypass()

set_seed!(10)
const testflts = [ _erfinv( rand() ) for i in 1:200_000 ]

function native_float()
    for v in testflts
        @sprintf( "%10.4f", v)
    end
end
@static if VERSION >= v"1.6"
function format_float()
    fmt = Printf.Format( "%10.4f" )
    for v in testflts
        Printf.format( fmt, v )
    end
end
end
function runtime_float()
    for v in testflts
        cfmt( "%10.4f", v)
    end
end
function runtime_float_spec()
    fs = Format.FmtSpec("%10.4f")
    for v in testflts
        cfmt( fs, v )
    end
end
function runtime_float_bypass()
    f = generate_formatter( "%10.4f" )
    for v in testflts
        f( v )
    end
end

println()
println( "float64 @sprintf speed")
@time native_float()
@static if VERSION >= v"1.6"
println( "float64 format speed")
@time format_float()
end
println( "float64 cfmt speed")
@time runtime_float()
println( "float64 cfmt spec speed")
@time runtime_float_spec()
println( "float64 cfmt speed, bypass repeated lookup")
@time runtime_float_bypass()

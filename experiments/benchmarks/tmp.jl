using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--profiler"
        help = "Profiler option: nsight or flamegraph"
        arg_type = String
        default = "flamegraph"
    end
    return parse_args(s)
end
parsed_args = parse_commandline()

@info parsed_args["profiler"]

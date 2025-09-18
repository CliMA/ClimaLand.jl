

# (filename::String, is_stdout_log::Bool)
file_list = map(f -> (f, !isnothing(match(r"scale\.o\d*", basename(f)))), readdir())

scaling_for_res = Dict{Float64, Array}()

for (fname, is_stdout_log) in file_list
    if is_stdout_log
        @show fname
        filestr = read(fname, String)
        @show filestr
        rx = r"("
        for match in eachmatch(r"\(\d\.\d* , \d* , \d*\.\d*\)", filestr)
        end
    end
end

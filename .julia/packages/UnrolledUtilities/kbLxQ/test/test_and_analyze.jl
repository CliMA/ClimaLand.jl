using Test
using JET
using OrderedCollections
using PrettyTables
using InteractiveUtils

using UnrolledUtilities
include("recursively_unrolled_functions.jl")

comparison_table_dicts = OrderedDict()

function print_comparison_table(title, comparison_table_dict, io = stdout)
    table_data =
        mapreduce(vcat, collect(comparison_table_dict)) do (key, entries)
            stack(entry -> (key..., entry...), entries; dims = 1)
        end

    writing_to_docs = io isa IOStream

    color(color_str) =
        writing_to_docs ? HtmlDecoration(; color = color_str) :
        Crayon(; foreground = Symbol(color_str))
    highlighter_color(optimization, run_time, compile_time, allocs) =
        if contains(optimization, "better") ||
           contains(optimization, "fewer allocs") &&
           !contains(run_time, "more") ||
           contains(optimization, "similar") && contains(run_time, "less")
            # better performance
            if !contains(run_time, "more") &&
               !contains(compile_time, "more") &&
               !contains(allocs, "more")
                # similar or better run time, compilation, and total allocations
                if contains(optimization, "better")
                    # better optimization
                    color(writing_to_docs ? "darkturquoise" : "cyan")
                else
                    # faster run time or fewer allocations at run time
                    color(writing_to_docs ? "mediumseagreen" : "green")
                end
            else
                # worse run time, compilation, or total allocations
                if contains(optimization, "better")
                    # better optimization
                    color(writing_to_docs ? "royalblue" : "blue")
                else
                    # faster run time or fewer allocations at run time
                    color(writing_to_docs ? "khaki" : "yellow")
                end
            end
        elseif contains(optimization, "similar") &&
               contains(run_time, "similar")
            # similar performance
            if contains(compile_time, "less") && !contains(allocs, "more") ||
               !contains(compile_time, "more") && contains(allocs, "less")
                # better compilation or total allocations
                color(writing_to_docs ? "mediumorchid" : "magenta")
            elseif contains(compile_time, "less") && contains(allocs, "more") ||
                   contains(compile_time, "more") && contains(allocs, "less")
                # mixed compilation and total allocations
                color(writing_to_docs ? "silver" : "light_gray")
            elseif contains(compile_time, "similar") &&
                   contains(allocs, "similar")
                # similar compilation and total allocations
                color(writing_to_docs ? "gray" : "dark_gray")
            else
                # worse compilation or total allocations
                color(writing_to_docs ? "indianred" : "red")
            end
        else
            # worse performance
            color(writing_to_docs ? "indianred" : "red")
        end
    highlighter = (writing_to_docs ? HtmlHighlighter : Highlighter)(
        Returns(true),
        (_, data, row, _) -> highlighter_color(data[row, 6:9]...),
    )

    # TODO: Why does Sys.maxrss() always seem to be 0 on Ubuntu systems?
    has_rss = any(contains('['), table_data[:, 9])

    other_kwargs =
        writing_to_docs ?
        (;
            backend = Val(:html),
            table_style = Dict(
                "font-family" => "monospace",
                "font-size" => "70%",
            ),
        ) :
        (;
            title,
            title_alignment = :c,
            title_same_width_as_table = true,
            columns_width = [45, 45, 15, 10, 30, 25, 20, 20, has_rss ? 30 : 20],
            linebreaks = true,
            autowrap = true,
            crop = :none,
        )

    if writing_to_docs
        println(io, "## $title")
        println(io, "```@raw html")
        println(io, "<div style=\"width: max(80vw, 100%)\">") # 80% of viewport
    end
    pretty_table(
        io,
        table_data;
        alignment = :l,
        header = [
            "Unrolled Expression",
            "Reference Expression",
            "Itr Type",
            "Itr Length",
            "Itr Contents",
            "Optimization",
            "Run Time",
            "Compilation Time",
            "Total $(has_rss ? "GC [and RSS] " : "")Allocations",
        ],
        highlighters = highlighter,
        other_kwargs...,
    )
    if writing_to_docs
        println(io, "</div>")
        println(io, "```")
    else
        println(io)
    end
end

function time_string(nanoseconds)
    nanoseconds == 0 && return "$nanoseconds ns"
    n_decimal_digits = floor(Int, log10(nanoseconds) + 1)
    return if n_decimal_digits <= 3
        "$nanoseconds ns"
    elseif n_decimal_digits <= 6
        "$(round(Int, nanoseconds / 10^3)) μs"
    elseif n_decimal_digits <= 9
        "$(round(Int, nanoseconds / 10^6)) ms"
    else
        "$(round(Int, nanoseconds / 10^9)) s"
    end
end

function memory_string(bytes)
    bytes == 0 && return "$bytes B"
    n_binary_digits = floor(Int, log2(bytes) + 1)
    return if n_binary_digits <= 10
        "$bytes B"
    elseif n_binary_digits <= 20
        "$(round(Int, bytes / 2^10)) kB"
    elseif n_binary_digits <= 30
        "$(round(Int, bytes / 2^20)) MB"
    else
        "$(round(Int, bytes / 2^30)) GB"
    end
end

function comparison_string(
    value1,
    value2,
    to_string;
    to_number = identity,
    epsilon = 0,
)
    number1 = to_number(value1)
    number2 = to_number(value2)
    ratio = number1 / number2
    ratio_str =
        if ratio < 2 && inv(ratio) < 2 || number1 < epsilon && number2 < epsilon
            "similar"
        elseif ratio >= 2
            floored_ratio = ratio == Inf ? Inf : floor(Int, ratio)
            "$floored_ratio times more"
        else
            floored_inv_ratio = ratio == 0 ? Inf : floor(Int, inv(ratio))
            "$floored_inv_ratio times less"
        end
    return "$ratio_str ($(to_string(value1)) vs. $(to_string(value2)))"
end

function drop_line_numbers(expr)
    expr isa Expr || return expr
    new_args = map(drop_line_numbers, expr.args)
    expr.head == :block || return Expr(expr.head, new_args...)
    filtered_args = filter(arg -> !(arg isa LineNumberNode), new_args)
    return length(filtered_args) == 1 ? filtered_args[1] :
           Expr(expr.head, filtered_args...)
end

simplified_expression_string(expr) =
    replace(string(drop_line_numbers(expr)), r"#=.+=# @" => '@', r"\s+" => ' ')

function code_instance(f, args...)
    available_methods = methods(f, Tuple{map(typeof, args)...})
    @assert length(available_methods) == 1
    (; specializations) = available_methods[1]
    specTypes = Tuple{typeof(f), map(typeof, args)...}
    return if specializations isa Core.MethodInstance
        @assert specializations.specTypes == specTypes
        specializations.cache
    else
        matching_specialization_indices =
            findall(specializations) do specialization
                !isnothing(specialization) &&
                    specialization.specTypes == specTypes
            end
        @assert length(matching_specialization_indices) == 1
        specializations[matching_specialization_indices[1]].cache
    end
end

macro benchmark(expression)
    return quote
        prev_time = time_ns()
        $(esc(expression))
        new_time = time_ns()
        best_time = new_time - prev_time

        # Benchmark for at most 0.1 s (10^8 ns), ignoring the first call above.
        n_trials = 0
        start_time = new_time
        while n_trials < 10^4 && new_time - start_time < 10^8
            prev_time = time_ns()
            $(esc(expression))
            new_time = time_ns()
            best_time = min(best_time, new_time - prev_time)
            n_trials += 1
        end

        best_time
    end
end

macro test_unrolled(
    args_expr,
    unrolled_expr,
    reference_expr,
    itr_contents_str,
    skip_allocations_test = false,
    skip_type_stability_test = false,
    load_recursively_unrolled_functions = false,
)
    @assert Meta.isexpr(args_expr, :tuple)
    arg_names = args_expr.args
    @assert all(arg_name -> arg_name isa Symbol, arg_names)
    args = map(esc, arg_names)
    unrolled_expr_str = simplified_expression_string(unrolled_expr)
    reference_expr_str = simplified_expression_string(reference_expr)
    contains_str = length(args) == 1 ? " that contains" : "s that each contain"
    quote
        itr_types = map(arg -> typeof(arg).name.wrapper, ($(args...),))
        itr_lengths = map(length, ($(args...),))

        itr_type_str =
            length(unique(itr_types)) == 1 ? string(itr_types[1]) :
            join(itr_types, '/')
        itr_length_str =
            length(unique(itr_lengths)) == 1 ? string(itr_lengths[1]) :
            join(itr_lengths, '/')
        itr_str =
            $(isempty(args)) ? "nothing" :
            "$($(length(args))) $itr_type_str$($contains_str) $itr_length_str \
             $($(esc(itr_contents_str)))"

        @info "Testing $($unrolled_expr_str) with $itr_str"

        unrolled_func($(arg_names...)) = $(esc(unrolled_expr))
        reference_func($(arg_names...)) = $(esc(reference_expr))

        # Test for correctness.
        @test unrolled_func($(args...)) == reference_func($(args...))

        unrolled_func_and_nothing($(arg_names...)) =
            ($(esc(unrolled_expr)); nothing)
        reference_func_and_nothing($(arg_names...)) =
            ($(esc(reference_expr)); nothing)

        unrolled_func_and_nothing($(args...)) # Run once to compile.
        reference_func_and_nothing($(args...))

        # Test for allocations.
        unrolled_run_memory = @allocated unrolled_func_and_nothing($(args...))
        reference_run_memory = @allocated reference_func_and_nothing($(args...))
        $(esc(skip_allocations_test)) || @test unrolled_run_memory == 0

        # Test for type-stability.
        is_unrolled_stable =
            isempty(JET.get_reports(@report_opt unrolled_func($(args...))))
        is_reference_stable =
            isempty(JET.get_reports(@report_opt reference_func($(args...))))
        $(esc(skip_type_stability_test)) || @test_opt unrolled_func($(args...))

        unrolled_code = code_instance(unrolled_func, $(args...))
        reference_code = code_instance(reference_func, $(args...))

        # Test for constant propagation.
        is_unrolled_const =
            isdefined(unrolled_code, :rettype_const) &&
            isbits(unrolled_code.rettype_const)
        is_reference_const =
            isdefined(reference_code, :rettype_const) &&
            isbits(reference_code.rettype_const)

        buffer = IOBuffer()
        args_type = Tuple{map(typeof, ($(args...),))...}
        code_llvm(buffer, unrolled_func, args_type; debuginfo = :none)
        unrolled_llvm = String(take!(buffer))
        code_llvm(buffer, reference_func, args_type; debuginfo = :none)
        reference_llvm = String(take!(buffer))
        code_llvm(buffer, Returns(nothing), Tuple{}; debuginfo = :none)
        no_op_llvm = String(take!(buffer))

        # Test whether the functions are compiled into switch instructions.
        is_unrolled_switch = contains(unrolled_llvm, "switch")
        is_reference_switch = contains(reference_llvm, "switch")

        # Test whether the functions are fully optimized out.
        unrolled_llvm_lines = length(split(unrolled_llvm, '\n'))
        reference_llvm_lines = length(split(reference_llvm, '\n'))
        no_op_llvm_lines = length(split(no_op_llvm, '\n'))
        is_unrolled_no_op = unrolled_llvm_lines == no_op_llvm_lines
        is_reference_no_op = reference_llvm_lines == no_op_llvm_lines

        # Test the overall level of optimization.
        unrolled_opt_str, unrolled_opt_score = if unrolled_run_memory > 0
            "$(memory_string(unrolled_run_memory)) allocs", 1 / unrolled_run_memory
        elseif !is_unrolled_stable
            "type-unstable", 2
        elseif !is_unrolled_const && !is_unrolled_switch
            "$unrolled_llvm_lines LLVM lines", 3
        elseif !is_unrolled_const
            "$unrolled_llvm_lines LLVM lines w/ switch", 3 # same as no switch
        elseif !is_unrolled_no_op
            "$unrolled_llvm_lines LLVM lines w/ Const", 4
        else
            "optimized out", 5
        end
        reference_opt_str, reference_opt_score = if reference_run_memory > 0
            "$(memory_string(reference_run_memory)) allocs",
            1 / reference_run_memory
        elseif !is_reference_stable
            "type-unstable", 2
        elseif !is_reference_const && !is_reference_switch
            "$reference_llvm_lines LLVM lines", 3
        elseif !is_reference_const
            "$reference_llvm_lines LLVM lines w/ switch", 3 # same as no switch
        elseif !is_reference_no_op
            "$reference_llvm_lines LLVM lines w/ Const", 4
        else
            "optimized out", 5
        end
        $(esc(skip_type_stability_test)) ||
            @test unrolled_opt_score >= reference_opt_score

        # Measure the run times.
        unrolled_run_time = @benchmark unrolled_func($(args...))
        reference_run_time = @benchmark reference_func($(args...))

        # Measure the compilation times and memory allocations in separate
        # processes to ensure that they are not under-counted.
        arg_name_strs = ($(map(string, arg_names)...),)
        arg_definition_strs =
            map((name, value) -> "$name = $value", arg_name_strs, ($(args...),))
        arg_definitions_str = join(arg_definition_strs, '\n')
        recursively_unrolled_functions_file_path = escape_string(
            joinpath(@__DIR__, "recursively_unrolled_functions.jl"),
        )
        load_recursively_unrolled_functions_str =
            $load_recursively_unrolled_functions ?
            "include(\"$recursively_unrolled_functions_file_path\")" : ""
        command_str(func_str) = """
            using UnrolledUtilities
            $arg_definitions_str
            $load_recursively_unrolled_functions_str
            Base.cumulative_compile_timing(true)
            nanoseconds1 = Base.cumulative_compile_time_ns()[1]
            rss_bytes_1 = Sys.maxrss()
            Δgc_bytes = @allocated $func_str
            rss_bytes_2 = Sys.maxrss()
            nanoseconds2 = Base.cumulative_compile_time_ns()[1]
            Base.cumulative_compile_timing(false)
            Δnanoseconds = nanoseconds2 - nanoseconds1
            Δrss_bytes = rss_bytes_2 - rss_bytes_1
            print(Δnanoseconds, ", ", Δgc_bytes, ", ", Δrss_bytes)
            """

        unrolled_command_str = command_str($(string(unrolled_expr)))
        run(pipeline(`julia --project -e $unrolled_command_str`, buffer))
        unrolled_compile_time, unrolled_total_memory, unrolled_total_rss =
            parse.((Int, Int, Int), split(String(take!(buffer)), ','))

        # Make a new buffer to avoid a potential data race:
        # discourse.julialang.org/t/iobuffer-becomes-not-writable-after-run/92323/3
        close(buffer)
        buffer = IOBuffer()

        reference_command_str = command_str($(string(reference_expr)))
        run(pipeline(`julia --project -e $reference_command_str`, buffer))
        reference_compile_time, reference_total_memory, reference_total_rss =
            parse.((Int, Int, Int), split(String(take!(buffer)), ','))

        close(buffer)

        optimization_str = if unrolled_opt_score > reference_opt_score
            if unrolled_opt_score <= 1
                "fewer allocs ($unrolled_opt_str vs. $reference_opt_str)"
            else
                "better ($unrolled_opt_str vs. $reference_opt_str)"
            end
        elseif unrolled_opt_score < reference_opt_score
            "worse ($unrolled_opt_str vs. $reference_opt_str)"
        else
            "similar ($unrolled_opt_str)"
        end
        run_time_str = comparison_string(
            unrolled_run_time,
            reference_run_time,
            time_string;
            epsilon = 50, # Ignore differences between times shorter than 50 ns.
        )
        compile_time_str = comparison_string(
            unrolled_compile_time,
            reference_compile_time,
            time_string,
        )
        memory_str = comparison_string(
            (unrolled_total_memory, unrolled_total_rss),
            (reference_total_memory, reference_total_rss),
            ((gc_bytes, rss_bytes),) ->
                rss_bytes == 0 ? memory_string(gc_bytes) :
                "$(memory_string(gc_bytes)) [$(memory_string(rss_bytes))]";
            to_number = first, # Use GC number since RSS might be unavailable.
        )

        dict_key = ($unrolled_expr_str, $reference_expr_str)
        dict_entry = (
            itr_type_str,
            itr_length_str,
            $(esc(itr_contents_str)),
            optimization_str,
            run_time_str,
            compile_time_str,
            memory_str,
        )
        if dict_key in keys(comparison_table_dict)
            push!(comparison_table_dict[dict_key], dict_entry)
        else
            comparison_table_dict[dict_key] = [dict_entry]
        end
    end
end

tuple_of_tuples(num_tuples, min_tuple_length, singleton, identical) =
    ntuple(num_tuples) do index
        tuple_length = min_tuple_length + (identical ? 0 : (index - 1) % 7)
        ntuple(singleton ? Val : identity, tuple_length)
    end
function tuples_of_tuples_contents_str(itrs...)
    str = ""
    all(itr -> length(itr) > 1 && length(unique(itr)) == 1, itrs) &&
        (str *= "identical ")
    all(itr -> length(itr) > 1 && length(unique(itr)) != 1, itrs) &&
        (str *= "distinct ")
    all(itr -> all(isempty, itr), itrs) && (str *= "empty ")
    all(itr -> all(!isempty, itr), itrs) && (str *= "nonempty ")
    all(itr -> any(isempty, itr) && any(!isempty, itr), itrs) &&
        (str *= "empty & nonempty ")
    all(itr -> Base.issingletontype(typeof(itr)), itrs) && (str *= "singleton ")
    all(itr -> !Base.issingletontype(typeof(itr)), itrs) &&
        (str *= "non-singleton ")
    str *= "Tuple"
    all(itr -> length(itr) > 1, itrs) && (str *= "s")
    return str
end

itr_lengths = "fast_mode" in ARGS ? (8,) : (2, 4, 8, 32, 128)

# NOTE: In the tests below, random numbers are meant to emulate values that
# cannot be inferred during compilation.

title = "Isolated Unrolled Functions"
comparison_table_dict = (comparison_table_dicts[title] = OrderedDict())

for itr in (
    tuple_of_tuples(1, 0, true, true),
    tuple_of_tuples(1, 1, true, true),
    tuple_of_tuples(1, 1, false, true),
    map(n -> tuple_of_tuples(n, 1, true, true), itr_lengths)...,
    map(n -> tuple_of_tuples(n, 1, false, true), itr_lengths)...,
    map(n -> tuple_of_tuples(n, 0, true, false), itr_lengths)...,
    map(n -> tuple_of_tuples(n, 1, true, false), itr_lengths)...,
    map(n -> tuple_of_tuples(n, 1, false, false), itr_lengths)...,
)
    str = tuples_of_tuples_contents_str(itr)
    itr_description = "a Tuple that contains $(length(itr)) $str"
    @testset "individual unrolled functions of $itr_description" begin
        @test_unrolled (itr,) unrolled_push(itr, itr[1]) (itr..., itr[1]) str

        @test_unrolled(
            (itr,),
            unrolled_append(itr, Iterators.reverse(itr)),
            (itr..., Iterators.reverse(itr)...),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_prepend(itr, Iterators.reverse(itr)),
            (Iterators.reverse(itr)..., itr...),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_take(itr, Val(length(itr) ÷ 2)),
            itr[1:(length(itr) ÷ 2)],
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_drop(itr, Val(length(itr) ÷ 2)),
            itr[(length(itr) ÷ 2 + 1):end],
            str,
        )

        @test_unrolled (itr,) unrolled_map(length, itr) map(length, itr) str

        @test_unrolled (itr,) unrolled_any(isempty, itr) any(isempty, itr) str
        @test_unrolled(
            (itr,),
            unrolled_any(x -> length(x) == rand(8:10), itr),
            any(x -> length(x) == rand(8:10), itr),
            str,
        )

        @test_unrolled (itr,) unrolled_all(isempty, itr) all(isempty, itr) str
        @test_unrolled(
            (itr,),
            unrolled_all(x -> length(x) == rand(8:10), itr),
            all(x -> length(x) == rand(8:10), itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_foreach(x -> @assert(length(x) <= 7), itr),
            foreach(x -> @assert(length(x) <= 7), itr),
            str,
        )

        @test_unrolled (itr,) unrolled_reduce(tuple, itr) reduce(tuple, itr) str
        @test_unrolled(
            (itr,),
            unrolled_reduce(tuple, itr; init = ()),
            reduce(tuple, itr; init = ()),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_mapreduce(x -> (x, length(x)), tuple, itr),
            mapreduce(x -> (x, length(x)), tuple, itr),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_mapreduce(x -> (x, length(x)), tuple, itr; init = ()),
            mapreduce(x -> (x, length(x)), tuple, itr; init = ()),
            str,
        )

        if length(itr) <= 32
            @test_unrolled(
                (itr,),
                unrolled_accumulate(tuple, itr),
                accumulate(tuple, itr),
                str,
            )
            @test_unrolled(
                (itr,),
                unrolled_accumulate(tuple, itr; init = ()),
                accumulate(tuple, itr; init = ()),
                str,
            )
        end # These can take half a minute to compile when the length is 128.

        @test_unrolled(
            (itr,),
            unrolled_applyat(length, rand(1:7:length(itr)), itr),
            length(itr[rand(1:7:length(itr))]),
            str,
        )

        @test_unrolled (itr,) unrolled_in(nothing, itr) (nothing in itr) str
        @test_unrolled (itr,) unrolled_in(itr[1], itr) (itr[1] in itr) str
        @test_unrolled (itr,) unrolled_in(itr[end], itr) (itr[end] in itr) str

        @test_unrolled(
            (itr,),
            unrolled_unique(length, itr),
            Tuple(unique(length, itr)),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_unique(itr),
            Tuple(unique(itr)),
            str,
            !Base.issingletontype(typeof(itr)),
            !Base.issingletontype(typeof(itr)),
        ) # unrolled_unique is only type-stable when comparing Core.Const data

        if VERSION >= v"1.11"
            @test_unrolled(
                (itr,),
                unrolled_allunique(length, itr),
                allunique(length, itr),
                str,
            )
            @test_unrolled(
                (itr,),
                unrolled_allequal(length, itr),
                allequal(length, itr),
                str,
            )
        else
            @test_unrolled(
                (itr,),
                unrolled_allunique(length, itr),
                allunique(Iterators.map(length, itr)),
                str,
            )
            @test_unrolled(
                (itr,),
                unrolled_allequal(length, itr),
                allequal(Iterators.map(length, itr)),
                str,
            )
        end # allunique(f, itr) and allequal(f, itr) require at least Julia 1.11

        @test_unrolled (itr,) unrolled_sum(length, itr) sum(length, itr) str
        @test_unrolled (itr,) unrolled_prod(length, itr) prod(length, itr) str

        @test_unrolled(
            (itr,),
            unrolled_cumsum(length, itr),
            Tuple(cumsum(Iterators.map(length, itr))),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_cumprod(length, itr),
            Tuple(cumprod(Iterators.map(length, itr))),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_count(!isempty, itr),
            count(!isempty, itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_maximum(length, itr),
            maximum(length, itr),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_minimum(length, itr),
            minimum(length, itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_extrema(length, itr),
            extrema(length, itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_findmax(length, itr),
            findmax(length, itr),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_findmin(length, itr),
            findmin(length, itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_argmax(length, itr),
            argmax(length, itr),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_argmin(length, itr),
            argmin(length, itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_argmax(unrolled_map(length, itr)),
            argmax(Iterators.map(length, itr)),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_argmin(unrolled_map(length, itr)),
            argmin(Iterators.map(length, itr)),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_findfirst(!isempty, itr),
            findfirst(!isempty, itr),
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_findlast(isempty, itr),
            findlast(isempty, itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_argfirst(x -> length(x) >= length(itr[end]), itr),
            itr[findfirst(x -> length(x) >= length(itr[end]), itr)],
            str,
        )
        @test_unrolled(
            (itr,),
            unrolled_arglast(x -> length(x) >= length(itr[1]), itr),
            itr[findlast(x -> length(x) >= length(itr[1]), itr)],
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_filter(!isempty, itr),
            filter(!isempty, itr),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_split(isempty, itr),
            (filter(isempty, itr), filter(!isempty, itr)),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_flatten(itr),
            Tuple(Iterators.flatten(itr)),
            str,
        )

        @test_unrolled(
            (itr,),
            unrolled_flatmap(Iterators.reverse, itr),
            Tuple(Iterators.flatmap(Iterators.reverse, itr)),
            str,
        )

        # TODO: Testing with coverage triggers allocations for unrolled_product!
        if length(itr) <= 32
            @test_unrolled(
                (itr,),
                unrolled_product(itr, itr),
                Tuple(Iterators.product(itr, itr)),
                str,
                "fast_mode" in ARGS,
                "fast_mode" in ARGS,
            )
        end # This can take several minutes to compile when the length is 128.
        if length(itr) <= 8
            @test_unrolled(
                (itr,),
                unrolled_product(itr, itr, itr),
                Tuple(Iterators.product(itr, itr, itr)),
                str,
                "fast_mode" in ARGS,
                "fast_mode" in ARGS,
            )
        end # This can take several minutes to compile when the length is 32.

        if VERSION >= v"1.11"
            @test_unrolled(
                (itr,),
                unrolled_cycle(itr, Val(3)),
                Tuple(Iterators.cycle(itr, 3)),
                str,
            )
        end # Iterators.cycle(itr, n) requires at least Julia 1.11

        @test_unrolled(
            (itr,),
            unrolled_partition(itr, Val(3)),
            Tuple(Iterators.map(Tuple, Iterators.partition(itr, 3))),
            str,
        )
    end
end

title = "Nested Unrolled Functions"
comparison_table_dict = (comparison_table_dicts[title] = OrderedDict())

for (itr1, itr2, itr3) in (
    (
        tuple_of_tuples(1, 0, true, true),
        tuple_of_tuples(1, 1, true, true),
        tuple_of_tuples(1, 1, false, true),
    ),
    zip(
        map(n -> tuple_of_tuples(n, 0, true, true), itr_lengths),
        map(n -> tuple_of_tuples(n, 1, true, true), itr_lengths),
        map(n -> tuple_of_tuples(n, 1, false, true), itr_lengths),
    )...,
    zip(
        map(n -> tuple_of_tuples(n, 0, true, false), itr_lengths),
        map(n -> tuple_of_tuples(n, 1, true, false), itr_lengths),
        map(n -> tuple_of_tuples(n, 1, false, false), itr_lengths),
    )...,
)
    str3 = tuples_of_tuples_contents_str(itr3)
    str12 = tuples_of_tuples_contents_str(itr1, itr2)
    str23 = tuples_of_tuples_contents_str(itr2, itr3)
    str123 = tuples_of_tuples_contents_str(itr1, itr2, itr3)
    itr_description = "Tuples that contain $(length(itr1)) $str123"
    @testset "nested unrolled functions of $itr_description" begin
        @test_unrolled(
            (itr3,),
            unrolled_any(x -> unrolled_sum(x) > 7, itr3),
            any(x -> sum(x) > 7, itr3),
            str3,
        )

        @test_unrolled(
            (itr3,),
            unrolled_mapreduce(unrolled_sum, max, itr3),
            mapreduce(sum, max, itr3),
            str3,
        )

        @test_unrolled(
            (itr3,),
            unrolled_applyat(unrolled_minimum, rand(1:length(itr3)), itr3),
            minimum(itr3[rand(1:length(itr3))]),
            str3,
        )

        @test_unrolled(
            (itr1, itr2),
            unrolled_foreach(
                (x1, x2) -> @assert(x1 == unrolled_take(x2, Val(length(x1)))),
                itr1,
                itr2,
            ),
            foreach((x1, x2) -> @assert(x1 == x2[1:length(x1)]), itr1, itr2),
            str12,
        )
        @test_unrolled(
            (itr2, itr3),
            unrolled_foreach(
                (x2, x3) -> @assert(x2 == unrolled_map(Val, x3)),
                itr2,
                itr3,
            ),
            foreach((x2, x3) -> @assert(x2 == map(Val, x3)), itr2, itr3),
            str23,
        )

        @test_unrolled(
            (itr1, itr2),
            unrolled_applyat(
                (x1, x2) -> @assert(x1 == unrolled_take(x2, Val(length(x1)))),
                rand(1:length(itr1)),
                itr1,
                itr2,
            ),
            let n = rand(1:length(itr1))
                @assert(itr1[n] == itr2[n][1:length(itr1[n])])
            end,
            str12,
        )
        @test_unrolled(
            (itr2, itr3),
            unrolled_applyat(
                (x2, x3) -> @assert(x2 == unrolled_map(Val, x3)),
                rand(1:length(itr2)),
                itr2,
                itr3,
            ),
            let n = rand(1:length(itr2))
                @assert(itr2[n] == map(Val, itr3[n]))
            end,
            str23,
        )
    end
end

nested_iterator(depth, n, inner_n) =
    depth == 1 ? ntuple(identity, n) :
    ntuple(Returns(nested_iterator(depth - 1, n ÷ inner_n, inner_n)), inner_n)

title = "Recursive Unrolled Functions"
comparison_table_dict = (comparison_table_dicts[title] = OrderedDict())

for n in itr_lengths
    @testset "recursive unrolled functions of $n values in nested Tuples" begin
        for depth in 2:2:(Int(log2(n)) + 1)
            itr = nested_iterator(depth, n, 2)
            str =
                depth > 2 ?
                "nested Tuples with $(n ÷ 2) values at depth $(depth - 1)" :
                (n > 2 ? "Tuples with $(n ÷ 2) values" : "Tuples with 1 value")
            # In the following definitions, use var"#self#" to avoid boxing:
            # discourse.julialang.org/t/performant-recursive-anonymous-functions/90984/5
            @test_unrolled(
                (itr,),
                map(
                    x ->
                        eltype(x) <: Tuple ? unrolled_sum(var"#self#", x) :
                        length(x),
                    (itr,),
                )[1],
                map(
                    x -> eltype(x) <: Tuple ? sum(var"#self#", x) : length(x),
                    (itr,),
                )[1],
                str,
            ) # nested iterator length
            @test_unrolled(
                (itr,),
                map(
                    x ->
                        eltype(x) <: Tuple ? unrolled_sum(var"#self#", x) :
                        unrolled_sum(x),
                    (itr,),
                )[1],
                map(
                    x -> eltype(x) <: Tuple ? sum(var"#self#", x) : sum(x),
                    (itr,),
                )[1],
                str,
            ) # nested iterator sum
            @test_unrolled(
                (itr,),
                map(
                    x ->
                        eltype(x) <: Tuple ?
                        unrolled_sum(var"#self#", x) +
                        unrolled_sum(log ∘ var"#self#", x) : unrolled_sum(x),
                    (itr,),
                )[1],
                map(
                    x ->
                        eltype(x) <: Tuple ?
                        sum(var"#self#", x) + sum(log ∘ var"#self#", x) :
                        sum(x),
                    (itr,),
                )[1],
                str,
            ) # recursive function that is improved by unrolling on Julia 1.11
        end
    end
end

title = "Nested Unrolled Closures"
comparison_table_dict = (comparison_table_dicts[title] = OrderedDict())

@testset "nested unrolled closures of Tuples vs. StaticBitVectors" begin
    for (itr, skip_allocations_test) in (
        (ntuple(Returns(true), 32), false),
        (ntuple(Returns(true), 33), true),
        (StaticBitVector{256}(true), false),
        (StaticBitVector{257}(true), true),
    )
        @test_unrolled(
            (itr,),
            unrolled_reduce(
                (itr′, i) -> Base.setindex(itr′, !itr′[i], i),
                StaticOneTo(length(itr));
                init = itr,
            ),
            reduce(
                (itr′, i) -> Base.setindex(itr′, !itr′[i], i),
                StaticOneTo(length(itr));
                init = itr,
            ),
            "Bools",
            skip_allocations_test,
        )
        @test_unrolled(
            (itr,),
            unrolled_reduce(
                (itr′, i) -> unrolled_reduce(
                    (itr′′, j) ->
                        Base.setindex(itr′′, !itr′′[min(i, j)], j),
                    StaticOneTo(length(itr′));
                    init = itr′,
                ),
                StaticOneTo(length(itr));
                init = itr,
            ),
            reduce(
                (itr′, i) -> reduce(
                    (itr′′, j) ->
                        Base.setindex(itr′′, !itr′′[min(i, j)], j),
                    StaticOneTo(length(itr′));
                    init = itr′,
                ),
                StaticOneTo(length(itr));
                init = itr,
            ),
            "Bools",
            skip_allocations_test,
        )
        if length(itr) <= 256
            @test_unrolled(
                (itr,),
                unrolled_reduce(
                    (itr′, i) -> unrolled_reduce(
                        (itr′′, j) -> unrolled_reduce(
                            (itr′′′, k) -> Base.setindex(
                                itr′′′,
                                !itr′′′[min(i, j, k)],
                                k,
                            ),
                            StaticOneTo(length(itr′′));
                            init = itr′′,
                        ),
                        StaticOneTo(length(itr′));
                        init = itr′,
                    ),
                    StaticOneTo(length(itr));
                    init = itr,
                ),
                reduce(
                    (itr′, i) -> reduce(
                        (itr′′, j) -> reduce(
                            (itr′′′, k) -> Base.setindex(
                                itr′′′,
                                !itr′′′[min(i, j, k)],
                                k,
                            ),
                            StaticOneTo(length(itr′′));
                            init = itr′′,
                        ),
                        StaticOneTo(length(itr′));
                        init = itr′,
                    ),
                    StaticOneTo(length(itr));
                    init = itr,
                ),
                "Bools",
                skip_allocations_test,
            )
        end # The StaticBitVector{257} allocates over 2 GB for this test.
    end
end

title = "Empty Iterators"
comparison_table_dict = (comparison_table_dicts[title] = OrderedDict())

@testset "unrolled functions of an empty Tuple" begin
    itr = ()
    str = "nothing"
    @test_unrolled (itr,) unrolled_map(error, itr) map(error, itr) str
    @test_unrolled (itr,) unrolled_any(error, itr) any(error, itr) str
    @test_unrolled (itr,) unrolled_all(error, itr) all(error, itr) str
    @test_unrolled (itr,) unrolled_foreach(error, itr) foreach(error, itr) str
    @test_throws "init" unrolled_reduce(error, itr)
    @test_unrolled(
        (itr,),
        unrolled_reduce(error, itr; init = 0),
        reduce(error, itr; init = 0),
        str,
    )
    @test_unrolled(
        (itr,),
        unrolled_accumulate(error, itr),
        accumulate(error, itr),
        str,
    )
    @test_unrolled(
        (itr,),
        unrolled_accumulate(error, itr; init = 0),
        accumulate(error, itr; init = 0),
        str,
    )
end

title = "Very Long Iterators"
comparison_table_dict = (comparison_table_dicts[title] = OrderedDict())

@testset "unrolled functions of Tuples vs. StaticOneTos" begin
    for itr in (ntuple(identity, 2000), StaticOneTo(2000), StaticOneTo(9000))
        @test_unrolled (itr,) unrolled_sum(itr) sum(itr) "Ints"
        @test_unrolled((itr,), unrolled_sum(log, itr), sum(log, itr), "Ints",)
    end # These can take over a minute to compile for ntuple(identity, 9000).
end

title = "Manual vs. Recursive Unrolling"
comparison_table_dict = (comparison_table_dicts[title] = OrderedDict())

man_vs_rec_itr_lengths = if "fast_mode" in ARGS
    (8,)
elseif VERSION >= v"1.11" && (Sys.iswindows() || Sys.isapple())
    # JET sometimes throws stack overflow errors for 256 elements on Julia 1.11.
    (1:8..., 16, 32, 128)
else
    (1:8..., 16, 32, 128, 256)
end

for itr in (
    tuple_of_tuples(1, 0, true, true),
    tuple_of_tuples(1, 1, true, true),
    tuple_of_tuples(1, 1, false, true),
    map(n -> tuple_of_tuples(n, 0, true, true), man_vs_rec_itr_lengths)...,
    map(n -> tuple_of_tuples(n, 1, true, true), man_vs_rec_itr_lengths)...,
    map(n -> tuple_of_tuples(n, 1, false, true), man_vs_rec_itr_lengths)...,
    map(n -> tuple_of_tuples(n, 0, true, false), man_vs_rec_itr_lengths)...,
    map(n -> tuple_of_tuples(n, 1, true, false), man_vs_rec_itr_lengths)...,
    map(n -> tuple_of_tuples(n, 1, false, false), man_vs_rec_itr_lengths)...,
)
    str = tuples_of_tuples_contents_str(itr)
    itr_description = "a Tuple that contains $(length(itr)) $str"
    @testset "manual vs. recursive unrolling of $itr_description" begin
        @test_unrolled(
            (itr,),
            UnrolledUtilities.unrolled_map_into_tuple(length, itr),
            rec_unrolled_map(length, itr),
            str,
            false,
            false,
            true,
        )

        @test_unrolled(
            (itr,),
            unrolled_any(isempty, itr),
            rec_unrolled_any(isempty, itr),
            str,
            false,
            false,
            true,
        )

        @test_unrolled(
            (itr,),
            unrolled_all(isempty, itr),
            rec_unrolled_all(isempty, itr),
            str,
            false,
            false,
            true,
        )

        @test_unrolled(
            (itr,),
            unrolled_foreach(x -> @assert(length(x) <= 7), itr),
            rec_unrolled_foreach(x -> @assert(length(x) <= 7), itr),
            str,
            false,
            false,
            true,
        )

        if length(itr) <= 32
            @test_unrolled(
                (itr,),
                unrolled_reduce(tuple, itr, ()),
                rec_unrolled_reduce(tuple, itr, ()),
                str,
                false,
                false,
                true,
            )

            @test_unrolled(
                (itr,),
                UnrolledUtilities.unrolled_accumulate_into_tuple(
                    tuple,
                    itr,
                    (),
                ),
                rec_unrolled_accumulate(tuple, itr, ()),
                str,
                false,
                false,
                true,
            )
        end # These can take over a minute to compile when the length is 128.

        @test_unrolled(
            (itr,),
            UnrolledUtilities.unrolled_ifelse(
                x -> length(x) >= length(itr[end]),
                identity,
                error,
                itr,
            ),
            rec_unrolled_ifelse(
                x -> length(x) >= length(itr[end]),
                identity,
                error,
                itr,
            ),
            str,
            false,
            false,
            true,
        ) # unrolled_argfirst(x -> length(x) >= length(itr[end]), itr)

        @test_unrolled(
            (itr,),
            UnrolledUtilities.unrolled_ifelse(
                !isempty,
                identity,
                Returns(nothing),
                itr,
                StaticOneTo(length(itr)),
            ),
            rec_unrolled_ifelse(
                !isempty,
                identity,
                Returns(nothing),
                itr,
                StaticOneTo(length(itr)),
            ),
            str,
            false,
            false,
            true,
        ) # unrolled_findfirst(!isempty, itr)
    end
end

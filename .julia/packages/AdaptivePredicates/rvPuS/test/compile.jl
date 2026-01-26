function define_forward_declarations(pred)
    declarations = raw"""
    double doublerand(void);
    double narrowdoublerand(void);
    double uniformdoublerand(void);
    float floatrand(void);
    float narrowfloatrand(void);
    float uniformfloatrand(void);
    void exactinit(void);
    int grow_expansion(int elen, REAL *e, REAL b, REAL *h);
    int grow_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h);
    int expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h);
    int expansion_sum_zeroelim1(int elen, REAL *e, int flen, REAL *f, REAL *h);
    int expansion_sum_zeroelim2(int elen, REAL *e, int flen, REAL *f, REAL *h);
    int fast_expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h);
    int fast_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h);
    int linear_expansion_sum(int elen, REAL *e, int flen, REAL *f, REAL *h);
    int linear_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h);
    int scale_expansion(int elen, REAL *e, REAL b, REAL *h);
    int scale_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h);
    int compress(int elen, REAL *e, REAL *h);
    REAL estimate(int elen, REAL *e);
    REAL orient2dfast(REAL *pa, REAL *pb, REAL *pc);
    REAL orient2dexact(REAL *pa, REAL *pb, REAL *pc);
    REAL orient2dslow(REAL *pa, REAL *pb, REAL *pc);
    REAL orient2dadapt(REAL *pa, REAL *pb, REAL *pc, REAL detsum);
    REAL orient2d(REAL *pa, REAL *pb, REAL *pc);
    REAL orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
    REAL orient3dexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
    REAL orient3dslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
    REAL orient3dadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL permanent);
    REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
    REAL incirclefast(REAL *pa, REAL*pb, REAL *pc, REAL *pd);
    REAL incircleexact(REAL *pa, REAL*pb, REAL *pc, REAL *pd);
    REAL incircleslow(REAL *pa, REAL*pb, REAL *pc, REAL *pd);
    REAL incircleadapt(REAL *pa, REAL*pb, REAL *pc, REAL *pd, REAL permanent);
    REAL incirclefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
    REAL incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
    REAL inspherefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
    REAL insphereexact(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
    REAL insphereslow(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
    REAL insphereadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe, REAL permanent);
    REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
    """
    # We put the declarations in some arbitrary place (after REAL has been defined, of course)
    return replace(pred,
        "#define UNIFORMRAND uniformdoublerand" => "#define UNIFORMRAND uniformdoublerand\n" * declarations
    )
end

function fix_windows(pred)
    if Sys.iswindows()
        # random(), being a POSIX function, is not available on Windows. We need to use rand().
        return replace(pred, "random()" => "rand()")
    else
        return pred
    end
end

function fix_two_square(pred)
    return replace(pred, "#define Two_Square" => "#define Old_Two_Square")
end

function convert_macros(pred)
    converted_macros = read("macro_defs.txt", String)
    return pred * "\n" * converted_macros
end

function convert_float(pred)
    return replace(
        pred,
        "REAL double" => "REAL float",
        "1.793662034335766e-43" => "1.793662034335766e-17",
        "3.2138760885179806e60" => "3.2138760885179806e12"
    )
end

function compile(pred)
    write("predicates64.c", pred)
    write("predicates32.c", convert_float(pred))
    if !isfile("libpredicates64.so")
        run(`gcc -shared -fPIC -g2 predicates64.c -o libpredicates64.so`)
    end
    if !isfile("libpredicates32.so")
        run(`gcc -shared -fPIC -g2 predicates32.c -o libpredicates32.so`)
    end
end

function test_macros()
    run(`gcc -o macro_tests macro_defs.c`)
    macro32 = read("macro_defs.c", String) |> convert_float
    write("macro_defs32.c", macro32)
    run(`gcc -o macro_tests32 macro_defs32.c`)
    @test success(run(ignorestatus(`./macro_tests`)))
    @test success(run(ignorestatus(`./macro_tests32`)))
end

function cleanup()
    rm("predicates32.c")
    rm("predicates64.c")
    rm("macro_defs32.c")
    rm("./macro_tests" * (Sys.iswindows() ? ".exe" : ""))
    rm("./macro_tests32" * (Sys.iswindows() ? ".exe" : ""))
end

function compile()
    read("predicates.c", String) |>
    define_forward_declarations |>
    fix_windows |>
    fix_two_square |>
    convert_macros |>
    compile
    test_macros()
    cleanup()
end

compile()
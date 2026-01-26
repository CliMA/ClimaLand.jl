let
    while true
        ccall(:jl_generating_output, Cint, ()) == 1 || break
        gl = GridLayout()
        gl2 = GridLayout()
        gl[1, 1] = gl2
        determinedirsize(gl, GridLayoutBase.Row())
        compute_rowcols(gl, GridLayoutBase.suggestedbboxobservable(gl)[])
        update!(gl)
        align_to_bbox!(gl2, GridLayoutBase.suggestedbboxobservable(gl2)[])
        compute_rowcols(gl2, GridLayoutBase.suggestedbboxobservable(gl2)[])
        break
    end
    nothing
end

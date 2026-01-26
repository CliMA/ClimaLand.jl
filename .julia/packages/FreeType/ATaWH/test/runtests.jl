using FreeType
using Pkg
using Test

library = Ref{FT_Library}()
error = FT_Init_FreeType(library)
@test error == 0

@testset "test FT_Outline_Decompose" begin

refface = Ref{FT_Face}()
@test FT_New_Face(library[], joinpath(@__DIR__, "hack_regular.ttf"), 0, refface) == 0

glyph_index = FT_Get_Char_Index(refface[], 'J')
@test glyph_index == 0x00000031
@test FT_Set_Char_Size(refface[], 0, 16*64, 3, 3) == 0
@test FT_Load_Glyph(refface[], glyph_index, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP) == 0

function pos(p::Ptr{FT_Vector})
    v = unsafe_load(p)
    (v.x, v.y)
end

paths = []
function move_to_func(to, user)
    push!(paths, (:move, pos(to)))
    Cint(0)
end
function line_to_func(to, user)
    push!(paths, (:line, pos(to)))
    Cint(0)
end
function conic_to_func(control, to, user)
    push!(paths, (:conic, pos.((control, to))...))
    Cint(0)
end
function cubic_to_func(control1, control2, to, user)
    push!(paths, (:cubic, pos.((control1, control2, to))...))
    Cint(0)
end

move_f = @cfunction $move_to_func Cint (Ptr{FT_Vector}, Ptr{Cvoid})
line_f = @cfunction $line_to_func Cint (Ptr{FT_Vector}, Ptr{Cvoid})
conic_f = @cfunction $conic_to_func Cint (Ptr{FT_Vector}, Ptr{FT_Vector}, Ptr{Cvoid})
cubic_f = @cfunction $cubic_to_func Cint (Ptr{FT_Vector}, Ptr{FT_Vector}, Ptr{FT_Vector}, Ptr{Cvoid})

GC.@preserve move_f line_f conic_f cubic_f begin
face = unsafe_load(refface[])
glyph = unsafe_load(face.glyph)
outline_funcs = FreeType.FT_Outline_Funcs(Base.unsafe_convert.(Ptr{Cvoid}, (move_f, line_f, conic_f, cubic_f))..., 0, 0)
FT_Outline_Decompose(pointer_from_objref.((Ref(glyph.outline), Ref(outline_funcs)))..., C_NULL)
end

@test paths == [(:move, (502, -29)), (:conic, (399, -29), (307, -7)), (:conic, (210, 16), (109, 61)), (:line, (109, 297)), (:conic, (200, 216), (297, 176)), (:conic, (396, 135), (499, 135)), (:conic, (566, 135), (617, 153)), (:conic, (669, 171), (698, 210)), (:conic, (754, 284), (754, 487)), (:line, (754, 1323)), (:line, (373, 1323)), (:line, (373, 1493)), (:line, (956, 1493)), (:line, (956, 487)), (:conic, (956, 343), (930, 246)), (:conic, (905, 150), (850, 89)), (:conic, (795, 28), (709, 0)), (:conic, (624, -29), (502, -29))]

@test FT_Done_FreeType(library[]) == 0

end


# since there are no meaningful tests, please manually do a test for FreeTypeAbstraction
# using FreeTypeAbstraction
# Pkg.test("FreeTypeAbstraction")

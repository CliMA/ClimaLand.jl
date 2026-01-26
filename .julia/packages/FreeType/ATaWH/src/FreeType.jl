module FreeType

using FreeType2_jll
export FreeType2_jll

using CEnum

const ptrdiff_t = Cptrdiff_t
const FILE = Libc.FILE
const SHRT_MAX = 32767
const USHRT_MAX = 65535
const CHAR_BIT = 8
const INT_MAX = typemax(Cint)
const INT_MIN = typemin(Cint)
const UINT_MAX = typemax(Cuint)
const UINT_MIN = typemin(Cuint)
const LONG_MAX = typemax(Clong)
const LONG_MIN = typemin(Clong)
const ULONG_MAX = typemax(Culong)


const FT_Int16 = Cshort

const FT_UInt16 = Cushort

const FT_Int32 = Cint

const FT_UInt32 = Cuint

const FT_Fast = Cint

const FT_UFast = Cuint

const FT_Int64 = Clong

const FT_UInt64 = Culong

# typedef void * ( * FT_Alloc_Func ) ( FT_Memory memory , long size )
const FT_Alloc_Func = Ptr{Cvoid}

# typedef void ( * FT_Free_Func ) ( FT_Memory memory , void * block )
const FT_Free_Func = Ptr{Cvoid}

# typedef void * ( * FT_Realloc_Func ) ( FT_Memory memory , long cur_size , long new_size , void * block )
const FT_Realloc_Func = Ptr{Cvoid}

struct FT_MemoryRec_
    user::Ptr{Cvoid}
    alloc::FT_Alloc_Func
    free::FT_Free_Func
    realloc::FT_Realloc_Func
end

const FT_Memory = Ptr{FT_MemoryRec_}

struct FT_StreamDesc_
    data::NTuple{8, UInt8}
end

function Base.getproperty(x::Ptr{FT_StreamDesc_}, f::Symbol)
    f === :value && return Ptr{Clong}(x + 0)
    f === :pointer && return Ptr{Ptr{Cvoid}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::FT_StreamDesc_, f::Symbol)
    r = Ref{FT_StreamDesc_}(x)
    ptr = Base.unsafe_convert(Ptr{FT_StreamDesc_}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{FT_StreamDesc_}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

const FT_StreamDesc = FT_StreamDesc_

# typedef unsigned long ( * FT_Stream_IoFunc ) ( FT_Stream stream , unsigned long offset , unsigned char * buffer , unsigned long count )
const FT_Stream_IoFunc = Ptr{Cvoid}

# typedef void ( * FT_Stream_CloseFunc ) ( FT_Stream stream )
const FT_Stream_CloseFunc = Ptr{Cvoid}

struct FT_StreamRec_
    base::Ptr{Cuchar}
    size::Culong
    pos::Culong
    descriptor::FT_StreamDesc
    pathname::FT_StreamDesc
    read::FT_Stream_IoFunc
    close::FT_Stream_CloseFunc
    memory::FT_Memory
    cursor::Ptr{Cuchar}
    limit::Ptr{Cuchar}
end

const FT_Stream = Ptr{FT_StreamRec_}

const FT_StreamRec = FT_StreamRec_

const FT_Pos = Clong

struct FT_Vector_
    x::FT_Pos
    y::FT_Pos
end

const FT_Vector = FT_Vector_

struct FT_BBox_
    xMin::FT_Pos
    yMin::FT_Pos
    xMax::FT_Pos
    yMax::FT_Pos
end

const FT_BBox = FT_BBox_

@cenum FT_Pixel_Mode_::UInt32 begin
    FT_PIXEL_MODE_NONE = 0
    FT_PIXEL_MODE_MONO = 1
    FT_PIXEL_MODE_GRAY = 2
    FT_PIXEL_MODE_GRAY2 = 3
    FT_PIXEL_MODE_GRAY4 = 4
    FT_PIXEL_MODE_LCD = 5
    FT_PIXEL_MODE_LCD_V = 6
    FT_PIXEL_MODE_BGRA = 7
    FT_PIXEL_MODE_MAX = 8
end

const FT_Pixel_Mode = FT_Pixel_Mode_

struct FT_Bitmap_
    rows::Cuint
    width::Cuint
    pitch::Cint
    buffer::Ptr{Cuchar}
    num_grays::Cushort
    pixel_mode::Cuchar
    palette_mode::Cuchar
    palette::Ptr{Cvoid}
end

const FT_Bitmap = FT_Bitmap_

struct FT_Outline_
    n_contours::Cshort
    n_points::Cshort
    points::Ptr{FT_Vector}
    tags::Ptr{Cchar}
    contours::Ptr{Cshort}
    flags::Cint
end

const FT_Outline = FT_Outline_

# typedef int ( * FT_Outline_MoveToFunc ) ( const FT_Vector * to , void * user )
const FT_Outline_MoveToFunc = Ptr{Cvoid}

# typedef int ( * FT_Outline_LineToFunc ) ( const FT_Vector * to , void * user )
const FT_Outline_LineToFunc = Ptr{Cvoid}

# typedef int ( * FT_Outline_ConicToFunc ) ( const FT_Vector * control , const FT_Vector * to , void * user )
const FT_Outline_ConicToFunc = Ptr{Cvoid}

# typedef int ( * FT_Outline_CubicToFunc ) ( const FT_Vector * control1 , const FT_Vector * control2 , const FT_Vector * to , void * user )
const FT_Outline_CubicToFunc = Ptr{Cvoid}

struct FT_Outline_Funcs_
    move_to::FT_Outline_MoveToFunc
    line_to::FT_Outline_LineToFunc
    conic_to::FT_Outline_ConicToFunc
    cubic_to::FT_Outline_CubicToFunc
    shift::Cint
    delta::FT_Pos
end

const FT_Outline_Funcs = FT_Outline_Funcs_

@cenum FT_Glyph_Format_::UInt32 begin
    FT_GLYPH_FORMAT_NONE = 0
    FT_GLYPH_FORMAT_COMPOSITE = 1668246896
    FT_GLYPH_FORMAT_BITMAP = 1651078259
    FT_GLYPH_FORMAT_OUTLINE = 1869968492
    FT_GLYPH_FORMAT_PLOTTER = 1886154612
end

const FT_Glyph_Format = FT_Glyph_Format_

mutable struct FT_RasterRec_ end

const FT_Raster = Ptr{FT_RasterRec_}

struct FT_Span_
    x::Cshort
    len::Cushort
    coverage::Cuchar
end

const FT_Span = FT_Span_

# typedef void ( * FT_SpanFunc ) ( int y , int count , const FT_Span * spans , void * user )
const FT_SpanFunc = Ptr{Cvoid}

# typedef int ( * FT_Raster_BitTest_Func ) ( int y , int x , void * user )
const FT_Raster_BitTest_Func = Ptr{Cvoid}

# typedef void ( * FT_Raster_BitSet_Func ) ( int y , int x , void * user )
const FT_Raster_BitSet_Func = Ptr{Cvoid}

struct FT_Raster_Params_
    target::Ptr{FT_Bitmap}
    source::Ptr{Cvoid}
    flags::Cint
    gray_spans::FT_SpanFunc
    black_spans::FT_SpanFunc
    bit_test::FT_Raster_BitTest_Func
    bit_set::FT_Raster_BitSet_Func
    user::Ptr{Cvoid}
    clip_box::FT_BBox
end

const FT_Raster_Params = FT_Raster_Params_

# typedef int ( * FT_Raster_NewFunc ) ( void * memory , FT_Raster * raster )
const FT_Raster_NewFunc = Ptr{Cvoid}

# typedef void ( * FT_Raster_DoneFunc ) ( FT_Raster raster )
const FT_Raster_DoneFunc = Ptr{Cvoid}

# typedef void ( * FT_Raster_ResetFunc ) ( FT_Raster raster , unsigned char * pool_base , unsigned long pool_size )
const FT_Raster_ResetFunc = Ptr{Cvoid}

# typedef int ( * FT_Raster_SetModeFunc ) ( FT_Raster raster , unsigned long mode , void * args )
const FT_Raster_SetModeFunc = Ptr{Cvoid}

# typedef int ( * FT_Raster_RenderFunc ) ( FT_Raster raster , const FT_Raster_Params * params )
const FT_Raster_RenderFunc = Ptr{Cvoid}

struct FT_Raster_Funcs_
    glyph_format::FT_Glyph_Format
    raster_new::FT_Raster_NewFunc
    raster_reset::FT_Raster_ResetFunc
    raster_set_mode::FT_Raster_SetModeFunc
    raster_render::FT_Raster_RenderFunc
    raster_done::FT_Raster_DoneFunc
end

const FT_Raster_Funcs = FT_Raster_Funcs_

const FT_Bool = Cuchar

const FT_FWord = Cshort

const FT_UFWord = Cushort

const FT_Char = Int8

const FT_Byte = Cuchar

const FT_Bytes = Ptr{FT_Byte}

const FT_Tag = FT_UInt32

const FT_String = Cchar

const FT_Short = Cshort

const FT_UShort = Cushort

const FT_Int = Cint

const FT_UInt = Cuint

const FT_Long = Clong

const FT_ULong = Culong

const FT_F2Dot14 = Cshort

const FT_F26Dot6 = Clong

const FT_Fixed = Clong

const FT_Error = Cint

const FT_Pointer = Ptr{Cvoid}

const FT_Offset = Csize_t

const FT_PtrDist = Cptrdiff_t

struct FT_UnitVector_
    x::FT_F2Dot14
    y::FT_F2Dot14
end

const FT_UnitVector = FT_UnitVector_

struct FT_Matrix_
    xx::FT_Fixed
    xy::FT_Fixed
    yx::FT_Fixed
    yy::FT_Fixed
end

const FT_Matrix = FT_Matrix_

struct FT_Data_
    pointer::Ptr{FT_Byte}
    length::FT_Int
end

const FT_Data = FT_Data_

# typedef void ( * FT_Generic_Finalizer ) ( void * object )
const FT_Generic_Finalizer = Ptr{Cvoid}

struct FT_Generic_
    data::Ptr{Cvoid}
    finalizer::FT_Generic_Finalizer
end

const FT_Generic = FT_Generic_

mutable struct __JL_FT_ListNodeRec_
end

function Base.unsafe_load(x::Ptr{__JL_FT_ListNodeRec_})
    unsafe_load(Ptr{FT_ListNodeRec_}(x))
end

function Base.getproperty(x::Ptr{__JL_FT_ListNodeRec_}, f::Symbol)
    getproperty(Ptr{FT_ListNodeRec_}(x), f)
end

function Base.setproperty!(x::Ptr{__JL_FT_ListNodeRec_}, f::Symbol, v)
    setproperty!(Ptr{FT_ListNodeRec_}(x), f, v)
end

const FT_ListNode = Ptr{__JL_FT_ListNodeRec_}

struct FT_ListRec_
    head::FT_ListNode
    tail::FT_ListNode
end

const FT_List = Ptr{FT_ListRec_}

struct FT_ListNodeRec_
    prev::FT_ListNode
    next::FT_ListNode
    data::Ptr{Cvoid}
end

const FT_ListNodeRec = FT_ListNodeRec_

const FT_ListRec = FT_ListRec_

@cenum __JL_Ctag_7::UInt32 begin
    FT_Mod_Err_Base = 0
    FT_Mod_Err_Autofit = 0
    FT_Mod_Err_BDF = 0
    FT_Mod_Err_Bzip2 = 0
    FT_Mod_Err_Cache = 0
    FT_Mod_Err_CFF = 0
    FT_Mod_Err_CID = 0
    FT_Mod_Err_Gzip = 0
    FT_Mod_Err_LZW = 0
    FT_Mod_Err_OTvalid = 0
    FT_Mod_Err_PCF = 0
    FT_Mod_Err_PFR = 0
    FT_Mod_Err_PSaux = 0
    FT_Mod_Err_PShinter = 0
    FT_Mod_Err_PSnames = 0
    FT_Mod_Err_Raster = 0
    FT_Mod_Err_SFNT = 0
    FT_Mod_Err_Smooth = 0
    FT_Mod_Err_TrueType = 0
    FT_Mod_Err_Type1 = 0
    FT_Mod_Err_Type42 = 0
    FT_Mod_Err_Winfonts = 0
    FT_Mod_Err_GXvalid = 0
    FT_Mod_Err_Max = 1
end

@cenum __JL_Ctag_8::UInt32 begin
    FT_Err_Ok = 0
    FT_Err_Cannot_Open_Resource = 1
    FT_Err_Unknown_File_Format = 2
    FT_Err_Invalid_File_Format = 3
    FT_Err_Invalid_Version = 4
    FT_Err_Lower_Module_Version = 5
    FT_Err_Invalid_Argument = 6
    FT_Err_Unimplemented_Feature = 7
    FT_Err_Invalid_Table = 8
    FT_Err_Invalid_Offset = 9
    FT_Err_Array_Too_Large = 10
    FT_Err_Missing_Module = 11
    FT_Err_Missing_Property = 12
    FT_Err_Invalid_Glyph_Index = 16
    FT_Err_Invalid_Character_Code = 17
    FT_Err_Invalid_Glyph_Format = 18
    FT_Err_Cannot_Render_Glyph = 19
    FT_Err_Invalid_Outline = 20
    FT_Err_Invalid_Composite = 21
    FT_Err_Too_Many_Hints = 22
    FT_Err_Invalid_Pixel_Size = 23
    FT_Err_Invalid_Handle = 32
    FT_Err_Invalid_Library_Handle = 33
    FT_Err_Invalid_Driver_Handle = 34
    FT_Err_Invalid_Face_Handle = 35
    FT_Err_Invalid_Size_Handle = 36
    FT_Err_Invalid_Slot_Handle = 37
    FT_Err_Invalid_CharMap_Handle = 38
    FT_Err_Invalid_Cache_Handle = 39
    FT_Err_Invalid_Stream_Handle = 40
    FT_Err_Too_Many_Drivers = 48
    FT_Err_Too_Many_Extensions = 49
    FT_Err_Out_Of_Memory = 64
    FT_Err_Unlisted_Object = 65
    FT_Err_Cannot_Open_Stream = 81
    FT_Err_Invalid_Stream_Seek = 82
    FT_Err_Invalid_Stream_Skip = 83
    FT_Err_Invalid_Stream_Read = 84
    FT_Err_Invalid_Stream_Operation = 85
    FT_Err_Invalid_Frame_Operation = 86
    FT_Err_Nested_Frame_Access = 87
    FT_Err_Invalid_Frame_Read = 88
    FT_Err_Raster_Uninitialized = 96
    FT_Err_Raster_Corrupted = 97
    FT_Err_Raster_Overflow = 98
    FT_Err_Raster_Negative_Height = 99
    FT_Err_Too_Many_Caches = 112
    FT_Err_Invalid_Opcode = 128
    FT_Err_Too_Few_Arguments = 129
    FT_Err_Stack_Overflow = 130
    FT_Err_Code_Overflow = 131
    FT_Err_Bad_Argument = 132
    FT_Err_Divide_By_Zero = 133
    FT_Err_Invalid_Reference = 134
    FT_Err_Debug_OpCode = 135
    FT_Err_ENDF_In_Exec_Stream = 136
    FT_Err_Nested_DEFS = 137
    FT_Err_Invalid_CodeRange = 138
    FT_Err_Execution_Too_Long = 139
    FT_Err_Too_Many_Function_Defs = 140
    FT_Err_Too_Many_Instruction_Defs = 141
    FT_Err_Table_Missing = 142
    FT_Err_Horiz_Header_Missing = 143
    FT_Err_Locations_Missing = 144
    FT_Err_Name_Table_Missing = 145
    FT_Err_CMap_Table_Missing = 146
    FT_Err_Hmtx_Table_Missing = 147
    FT_Err_Post_Table_Missing = 148
    FT_Err_Invalid_Horiz_Metrics = 149
    FT_Err_Invalid_CharMap_Format = 150
    FT_Err_Invalid_PPem = 151
    FT_Err_Invalid_Vert_Metrics = 152
    FT_Err_Could_Not_Find_Context = 153
    FT_Err_Invalid_Post_Table_Format = 154
    FT_Err_Invalid_Post_Table = 155
    FT_Err_DEF_In_Glyf_Bytecode = 156
    FT_Err_Missing_Bitmap = 157
    FT_Err_Syntax_Error = 160
    FT_Err_Stack_Underflow = 161
    FT_Err_Ignore = 162
    FT_Err_No_Unicode_Glyph_Name = 163
    FT_Err_Glyph_Too_Big = 164
    FT_Err_Missing_Startfont_Field = 176
    FT_Err_Missing_Font_Field = 177
    FT_Err_Missing_Size_Field = 178
    FT_Err_Missing_Fontboundingbox_Field = 179
    FT_Err_Missing_Chars_Field = 180
    FT_Err_Missing_Startchar_Field = 181
    FT_Err_Missing_Encoding_Field = 182
    FT_Err_Missing_Bbx_Field = 183
    FT_Err_Bbx_Too_Big = 184
    FT_Err_Corrupted_Font_Header = 185
    FT_Err_Corrupted_Font_Glyphs = 186
    FT_Err_Max = 187
end

function FT_Error_String(error_code)
    ccall((:FT_Error_String, libfreetype), Ptr{Cchar}, (FT_Error,), error_code)
end

struct FT_Glyph_Metrics_
    width::FT_Pos
    height::FT_Pos
    horiBearingX::FT_Pos
    horiBearingY::FT_Pos
    horiAdvance::FT_Pos
    vertBearingX::FT_Pos
    vertBearingY::FT_Pos
    vertAdvance::FT_Pos
end

const FT_Glyph_Metrics = FT_Glyph_Metrics_

struct FT_Bitmap_Size_
    height::FT_Short
    width::FT_Short
    size::FT_Pos
    x_ppem::FT_Pos
    y_ppem::FT_Pos
end

const FT_Bitmap_Size = FT_Bitmap_Size_

mutable struct FT_LibraryRec_ end

const FT_Library = Ptr{FT_LibraryRec_}

mutable struct FT_ModuleRec_ end

const FT_Module = Ptr{FT_ModuleRec_}

mutable struct FT_DriverRec_ end

const FT_Driver = Ptr{FT_DriverRec_}

mutable struct FT_RendererRec_ end

const FT_Renderer = Ptr{FT_RendererRec_}

mutable struct __JL_FT_FaceRec_
end

function Base.unsafe_load(x::Ptr{__JL_FT_FaceRec_})
    unsafe_load(Ptr{FT_FaceRec_}(x))
end

function Base.getproperty(x::Ptr{__JL_FT_FaceRec_}, f::Symbol)
    getproperty(Ptr{FT_FaceRec_}(x), f)
end

function Base.setproperty!(x::Ptr{__JL_FT_FaceRec_}, f::Symbol, v)
    setproperty!(Ptr{FT_FaceRec_}(x), f, v)
end

const FT_Face = Ptr{__JL_FT_FaceRec_}

struct FT_Size_Metrics_
    x_ppem::FT_UShort
    y_ppem::FT_UShort
    x_scale::FT_Fixed
    y_scale::FT_Fixed
    ascender::FT_Pos
    descender::FT_Pos
    height::FT_Pos
    max_advance::FT_Pos
end

const FT_Size_Metrics = FT_Size_Metrics_

mutable struct FT_Size_InternalRec_ end

const FT_Size_Internal = Ptr{FT_Size_InternalRec_}

struct FT_SizeRec_
    face::FT_Face
    generic::FT_Generic
    metrics::FT_Size_Metrics
    internal::FT_Size_Internal
end

const FT_Size = Ptr{FT_SizeRec_}

mutable struct __JL_FT_GlyphSlotRec_
end

function Base.unsafe_load(x::Ptr{__JL_FT_GlyphSlotRec_})
    unsafe_load(Ptr{FT_GlyphSlotRec_}(x))
end

function Base.getproperty(x::Ptr{__JL_FT_GlyphSlotRec_}, f::Symbol)
    getproperty(Ptr{FT_GlyphSlotRec_}(x), f)
end

function Base.setproperty!(x::Ptr{__JL_FT_GlyphSlotRec_}, f::Symbol, v)
    setproperty!(Ptr{FT_GlyphSlotRec_}(x), f, v)
end

const FT_GlyphSlot = Ptr{__JL_FT_GlyphSlotRec_}

mutable struct __JL_FT_CharMapRec_
end

function Base.unsafe_load(x::Ptr{__JL_FT_CharMapRec_})
    unsafe_load(Ptr{FT_CharMapRec_}(x))
end

function Base.getproperty(x::Ptr{__JL_FT_CharMapRec_}, f::Symbol)
    getproperty(Ptr{FT_CharMapRec_}(x), f)
end

function Base.setproperty!(x::Ptr{__JL_FT_CharMapRec_}, f::Symbol, v)
    setproperty!(Ptr{FT_CharMapRec_}(x), f, v)
end

const FT_CharMap = Ptr{__JL_FT_CharMapRec_}

@cenum FT_Encoding_::UInt32 begin
    FT_ENCODING_NONE = 0
    FT_ENCODING_MS_SYMBOL = 1937337698
    FT_ENCODING_UNICODE = 1970170211
    FT_ENCODING_SJIS = 1936353651
    FT_ENCODING_PRC = 1734484000
    FT_ENCODING_BIG5 = 1651074869
    FT_ENCODING_WANSUNG = 2002873971
    FT_ENCODING_JOHAB = 1785686113
    FT_ENCODING_GB2312 = 1734484000
    FT_ENCODING_MS_SJIS = 1936353651
    FT_ENCODING_MS_GB2312 = 1734484000
    FT_ENCODING_MS_BIG5 = 1651074869
    FT_ENCODING_MS_WANSUNG = 2002873971
    FT_ENCODING_MS_JOHAB = 1785686113
    FT_ENCODING_ADOBE_STANDARD = 1094995778
    FT_ENCODING_ADOBE_EXPERT = 1094992453
    FT_ENCODING_ADOBE_CUSTOM = 1094992451
    FT_ENCODING_ADOBE_LATIN_1 = 1818326065
    FT_ENCODING_OLD_LATIN_2 = 1818326066
    FT_ENCODING_APPLE_ROMAN = 1634889070
end

const FT_Encoding = FT_Encoding_

struct FT_CharMapRec_
    face::FT_Face
    encoding::FT_Encoding
    platform_id::FT_UShort
    encoding_id::FT_UShort
end

const FT_CharMapRec = FT_CharMapRec_

mutable struct FT_Face_InternalRec_ end

const FT_Face_Internal = Ptr{FT_Face_InternalRec_}

struct FT_FaceRec_
    num_faces::FT_Long
    face_index::FT_Long
    face_flags::FT_Long
    style_flags::FT_Long
    num_glyphs::FT_Long
    family_name::Ptr{FT_String}
    style_name::Ptr{FT_String}
    num_fixed_sizes::FT_Int
    available_sizes::Ptr{FT_Bitmap_Size}
    num_charmaps::FT_Int
    charmaps::Ptr{FT_CharMap}
    generic::FT_Generic
    bbox::FT_BBox
    units_per_EM::FT_UShort
    ascender::FT_Short
    descender::FT_Short
    height::FT_Short
    max_advance_width::FT_Short
    max_advance_height::FT_Short
    underline_position::FT_Short
    underline_thickness::FT_Short
    glyph::FT_GlyphSlot
    size::FT_Size
    charmap::FT_CharMap
    driver::FT_Driver
    memory::FT_Memory
    stream::FT_Stream
    sizes_list::FT_ListRec
    autohint::FT_Generic
    extensions::Ptr{Cvoid}
    internal::FT_Face_Internal
end

const FT_FaceRec = FT_FaceRec_

const FT_SizeRec = FT_SizeRec_

mutable struct FT_SubGlyphRec_ end

const FT_SubGlyph = Ptr{FT_SubGlyphRec_}

mutable struct FT_Slot_InternalRec_ end

const FT_Slot_Internal = Ptr{FT_Slot_InternalRec_}

struct FT_GlyphSlotRec_
    library::FT_Library
    face::FT_Face
    next::FT_GlyphSlot
    glyph_index::FT_UInt
    generic::FT_Generic
    metrics::FT_Glyph_Metrics
    linearHoriAdvance::FT_Fixed
    linearVertAdvance::FT_Fixed
    advance::FT_Vector
    format::FT_Glyph_Format
    bitmap::FT_Bitmap
    bitmap_left::FT_Int
    bitmap_top::FT_Int
    outline::FT_Outline
    num_subglyphs::FT_UInt
    subglyphs::FT_SubGlyph
    control_data::Ptr{Cvoid}
    control_len::Clong
    lsb_delta::FT_Pos
    rsb_delta::FT_Pos
    other::Ptr{Cvoid}
    internal::FT_Slot_Internal
end

const FT_GlyphSlotRec = FT_GlyphSlotRec_

function FT_Init_FreeType(alibrary)
    ccall((:FT_Init_FreeType, libfreetype), FT_Error, (Ptr{FT_Library},), alibrary)
end

function FT_Done_FreeType(library)
    ccall((:FT_Done_FreeType, libfreetype), FT_Error, (FT_Library,), library)
end

struct FT_Parameter_
    tag::FT_ULong
    data::FT_Pointer
end

const FT_Parameter = FT_Parameter_

struct FT_Open_Args_
    flags::FT_UInt
    memory_base::Ptr{FT_Byte}
    memory_size::FT_Long
    pathname::Ptr{FT_String}
    stream::FT_Stream
    driver::FT_Module
    num_params::FT_Int
    params::Ptr{FT_Parameter}
end

const FT_Open_Args = FT_Open_Args_

function FT_New_Face(library, filepathname, face_index, aface)
    ccall((:FT_New_Face, libfreetype), FT_Error, (FT_Library, Ptr{Cchar}, FT_Long, Ptr{FT_Face}), library, filepathname, face_index, aface)
end

function FT_New_Memory_Face(library, file_base, file_size, face_index, aface)
    ccall((:FT_New_Memory_Face, libfreetype), FT_Error, (FT_Library, Ptr{FT_Byte}, FT_Long, FT_Long, Ptr{FT_Face}), library, file_base, file_size, face_index, aface)
end

function FT_Open_Face(library, args, face_index, aface)
    ccall((:FT_Open_Face, libfreetype), FT_Error, (FT_Library, Ptr{FT_Open_Args}, FT_Long, Ptr{FT_Face}), library, args, face_index, aface)
end

function FT_Attach_File(face, filepathname)
    ccall((:FT_Attach_File, libfreetype), FT_Error, (FT_Face, Ptr{Cchar}), face, filepathname)
end

function FT_Attach_Stream(face, parameters)
    ccall((:FT_Attach_Stream, libfreetype), FT_Error, (FT_Face, Ptr{FT_Open_Args}), face, parameters)
end

function FT_Reference_Face(face)
    ccall((:FT_Reference_Face, libfreetype), FT_Error, (FT_Face,), face)
end

function FT_Done_Face(face)
    ccall((:FT_Done_Face, libfreetype), FT_Error, (FT_Face,), face)
end

function FT_Select_Size(face, strike_index)
    ccall((:FT_Select_Size, libfreetype), FT_Error, (FT_Face, FT_Int), face, strike_index)
end

@cenum FT_Size_Request_Type_::UInt32 begin
    FT_SIZE_REQUEST_TYPE_NOMINAL = 0
    FT_SIZE_REQUEST_TYPE_REAL_DIM = 1
    FT_SIZE_REQUEST_TYPE_BBOX = 2
    FT_SIZE_REQUEST_TYPE_CELL = 3
    FT_SIZE_REQUEST_TYPE_SCALES = 4
    FT_SIZE_REQUEST_TYPE_MAX = 5
end

const FT_Size_Request_Type = FT_Size_Request_Type_

struct FT_Size_RequestRec_
    type::FT_Size_Request_Type
    width::FT_Long
    height::FT_Long
    horiResolution::FT_UInt
    vertResolution::FT_UInt
end

const FT_Size_RequestRec = FT_Size_RequestRec_

const FT_Size_Request = Ptr{FT_Size_RequestRec_}

function FT_Request_Size(face, req)
    ccall((:FT_Request_Size, libfreetype), FT_Error, (FT_Face, FT_Size_Request), face, req)
end

function FT_Set_Char_Size(face, char_width, char_height, horz_resolution, vert_resolution)
    ccall((:FT_Set_Char_Size, libfreetype), FT_Error, (FT_Face, FT_F26Dot6, FT_F26Dot6, FT_UInt, FT_UInt), face, char_width, char_height, horz_resolution, vert_resolution)
end

function FT_Set_Pixel_Sizes(face, pixel_width, pixel_height)
    ccall((:FT_Set_Pixel_Sizes, libfreetype), FT_Error, (FT_Face, FT_UInt, FT_UInt), face, pixel_width, pixel_height)
end

function FT_Load_Glyph(face, glyph_index, load_flags)
    ccall((:FT_Load_Glyph, libfreetype), FT_Error, (FT_Face, FT_UInt, FT_Int32), face, glyph_index, load_flags)
end

function FT_Load_Char(face, char_code, load_flags)
    ccall((:FT_Load_Char, libfreetype), FT_Error, (FT_Face, FT_ULong, FT_Int32), face, char_code, load_flags)
end

function FT_Set_Transform(face, matrix, delta)
    ccall((:FT_Set_Transform, libfreetype), Cvoid, (FT_Face, Ptr{FT_Matrix}, Ptr{FT_Vector}), face, matrix, delta)
end

@cenum FT_Render_Mode_::UInt32 begin
    FT_RENDER_MODE_NORMAL = 0
    FT_RENDER_MODE_LIGHT = 1
    FT_RENDER_MODE_MONO = 2
    FT_RENDER_MODE_LCD = 3
    FT_RENDER_MODE_LCD_V = 4
    FT_RENDER_MODE_MAX = 5
end

const FT_Render_Mode = FT_Render_Mode_

function FT_Render_Glyph(slot, render_mode)
    ccall((:FT_Render_Glyph, libfreetype), FT_Error, (FT_GlyphSlot, FT_Render_Mode), slot, render_mode)
end

@cenum FT_Kerning_Mode_::UInt32 begin
    FT_KERNING_DEFAULT = 0
    FT_KERNING_UNFITTED = 1
    FT_KERNING_UNSCALED = 2
end

const FT_Kerning_Mode = FT_Kerning_Mode_

function FT_Get_Kerning(face, left_glyph, right_glyph, kern_mode, akerning)
    ccall((:FT_Get_Kerning, libfreetype), FT_Error, (FT_Face, FT_UInt, FT_UInt, FT_UInt, Ptr{FT_Vector}), face, left_glyph, right_glyph, kern_mode, akerning)
end

function FT_Get_Track_Kerning(face, point_size, degree, akerning)
    ccall((:FT_Get_Track_Kerning, libfreetype), FT_Error, (FT_Face, FT_Fixed, FT_Int, Ptr{FT_Fixed}), face, point_size, degree, akerning)
end

function FT_Get_Glyph_Name(face, glyph_index, buffer, buffer_max)
    ccall((:FT_Get_Glyph_Name, libfreetype), FT_Error, (FT_Face, FT_UInt, FT_Pointer, FT_UInt), face, glyph_index, buffer, buffer_max)
end

function FT_Get_Postscript_Name(face)
    ccall((:FT_Get_Postscript_Name, libfreetype), Ptr{Cchar}, (FT_Face,), face)
end

function FT_Select_Charmap(face, encoding)
    ccall((:FT_Select_Charmap, libfreetype), FT_Error, (FT_Face, FT_Encoding), face, encoding)
end

function FT_Set_Charmap(face, charmap)
    ccall((:FT_Set_Charmap, libfreetype), FT_Error, (FT_Face, FT_CharMap), face, charmap)
end

function FT_Get_Charmap_Index(charmap)
    ccall((:FT_Get_Charmap_Index, libfreetype), FT_Int, (FT_CharMap,), charmap)
end

function FT_Get_Char_Index(face, charcode)
    ccall((:FT_Get_Char_Index, libfreetype), FT_UInt, (FT_Face, FT_ULong), face, charcode)
end

function FT_Get_First_Char(face, agindex)
    ccall((:FT_Get_First_Char, libfreetype), FT_ULong, (FT_Face, Ptr{FT_UInt}), face, agindex)
end

function FT_Get_Next_Char(face, char_code, agindex)
    ccall((:FT_Get_Next_Char, libfreetype), FT_ULong, (FT_Face, FT_ULong, Ptr{FT_UInt}), face, char_code, agindex)
end

function FT_Face_Properties(face, num_properties, properties)
    ccall((:FT_Face_Properties, libfreetype), FT_Error, (FT_Face, FT_UInt, Ptr{FT_Parameter}), face, num_properties, properties)
end

function FT_Get_Name_Index(face, glyph_name)
    ccall((:FT_Get_Name_Index, libfreetype), FT_UInt, (FT_Face, Ptr{FT_String}), face, glyph_name)
end

function FT_Get_SubGlyph_Info(glyph, sub_index, p_index, p_flags, p_arg1, p_arg2, p_transform)
    ccall((:FT_Get_SubGlyph_Info, libfreetype), FT_Error, (FT_GlyphSlot, FT_UInt, Ptr{FT_Int}, Ptr{FT_UInt}, Ptr{FT_Int}, Ptr{FT_Int}, Ptr{FT_Matrix}), glyph, sub_index, p_index, p_flags, p_arg1, p_arg2, p_transform)
end

struct FT_LayerIterator_
    num_layers::FT_UInt
    layer::FT_UInt
    p::Ptr{FT_Byte}
end

const FT_LayerIterator = FT_LayerIterator_

function FT_Get_Color_Glyph_Layer(face, base_glyph, aglyph_index, acolor_index, iterator)
    ccall((:FT_Get_Color_Glyph_Layer, libfreetype), FT_Bool, (FT_Face, FT_UInt, Ptr{FT_UInt}, Ptr{FT_UInt}, Ptr{FT_LayerIterator}), face, base_glyph, aglyph_index, acolor_index, iterator)
end

function FT_Get_FSType_Flags(face)
    ccall((:FT_Get_FSType_Flags, libfreetype), FT_UShort, (FT_Face,), face)
end

function FT_Face_GetCharVariantIndex(face, charcode, variantSelector)
    ccall((:FT_Face_GetCharVariantIndex, libfreetype), FT_UInt, (FT_Face, FT_ULong, FT_ULong), face, charcode, variantSelector)
end

function FT_Face_GetCharVariantIsDefault(face, charcode, variantSelector)
    ccall((:FT_Face_GetCharVariantIsDefault, libfreetype), FT_Int, (FT_Face, FT_ULong, FT_ULong), face, charcode, variantSelector)
end

function FT_Face_GetVariantSelectors(face)
    ccall((:FT_Face_GetVariantSelectors, libfreetype), Ptr{FT_UInt32}, (FT_Face,), face)
end

function FT_Face_GetVariantsOfChar(face, charcode)
    ccall((:FT_Face_GetVariantsOfChar, libfreetype), Ptr{FT_UInt32}, (FT_Face, FT_ULong), face, charcode)
end

function FT_Face_GetCharsOfVariant(face, variantSelector)
    ccall((:FT_Face_GetCharsOfVariant, libfreetype), Ptr{FT_UInt32}, (FT_Face, FT_ULong), face, variantSelector)
end

function FT_MulDiv(a, b, c)
    ccall((:FT_MulDiv, libfreetype), FT_Long, (FT_Long, FT_Long, FT_Long), a, b, c)
end

function FT_MulFix(a, b)
    ccall((:FT_MulFix, libfreetype), FT_Long, (FT_Long, FT_Long), a, b)
end

function FT_DivFix(a, b)
    ccall((:FT_DivFix, libfreetype), FT_Long, (FT_Long, FT_Long), a, b)
end

function FT_RoundFix(a)
    ccall((:FT_RoundFix, libfreetype), FT_Fixed, (FT_Fixed,), a)
end

function FT_CeilFix(a)
    ccall((:FT_CeilFix, libfreetype), FT_Fixed, (FT_Fixed,), a)
end

function FT_FloorFix(a)
    ccall((:FT_FloorFix, libfreetype), FT_Fixed, (FT_Fixed,), a)
end

function FT_Vector_Transform(vector, matrix)
    ccall((:FT_Vector_Transform, libfreetype), Cvoid, (Ptr{FT_Vector}, Ptr{FT_Matrix}), vector, matrix)
end

function FT_Library_Version(library, amajor, aminor, apatch)
    ccall((:FT_Library_Version, libfreetype), Cvoid, (FT_Library, Ptr{FT_Int}, Ptr{FT_Int}, Ptr{FT_Int}), library, amajor, aminor, apatch)
end

function FT_Face_CheckTrueTypePatents(face)
    ccall((:FT_Face_CheckTrueTypePatents, libfreetype), FT_Bool, (FT_Face,), face)
end

function FT_Face_SetUnpatentedHinting(face, value)
    ccall((:FT_Face_SetUnpatentedHinting, libfreetype), FT_Bool, (FT_Face, FT_Bool), face, value)
end

function FT_Outline_Decompose(outline, func_interface, user)
    ccall((:FT_Outline_Decompose, libfreetype), FT_Error, (Ptr{FT_Outline}, Ptr{FT_Outline_Funcs}, Ptr{Cvoid}), outline, func_interface, user)
end

function FT_Outline_New(library, numPoints, numContours, anoutline)
    ccall((:FT_Outline_New, libfreetype), FT_Error, (FT_Library, FT_UInt, FT_Int, Ptr{FT_Outline}), library, numPoints, numContours, anoutline)
end

function FT_Outline_Done(library, outline)
    ccall((:FT_Outline_Done, libfreetype), FT_Error, (FT_Library, Ptr{FT_Outline}), library, outline)
end

function FT_Outline_Check(outline)
    ccall((:FT_Outline_Check, libfreetype), FT_Error, (Ptr{FT_Outline},), outline)
end

function FT_Outline_Get_CBox(outline, acbox)
    ccall((:FT_Outline_Get_CBox, libfreetype), Cvoid, (Ptr{FT_Outline}, Ptr{FT_BBox}), outline, acbox)
end

function FT_Outline_Translate(outline, xOffset, yOffset)
    ccall((:FT_Outline_Translate, libfreetype), Cvoid, (Ptr{FT_Outline}, FT_Pos, FT_Pos), outline, xOffset, yOffset)
end

function FT_Outline_Copy(source, target)
    ccall((:FT_Outline_Copy, libfreetype), FT_Error, (Ptr{FT_Outline}, Ptr{FT_Outline}), source, target)
end

function FT_Outline_Transform(outline, matrix)
    ccall((:FT_Outline_Transform, libfreetype), Cvoid, (Ptr{FT_Outline}, Ptr{FT_Matrix}), outline, matrix)
end

function FT_Outline_Embolden(outline, strength)
    ccall((:FT_Outline_Embolden, libfreetype), FT_Error, (Ptr{FT_Outline}, FT_Pos), outline, strength)
end

function FT_Outline_EmboldenXY(outline, xstrength, ystrength)
    ccall((:FT_Outline_EmboldenXY, libfreetype), FT_Error, (Ptr{FT_Outline}, FT_Pos, FT_Pos), outline, xstrength, ystrength)
end

function FT_Outline_Reverse(outline)
    ccall((:FT_Outline_Reverse, libfreetype), Cvoid, (Ptr{FT_Outline},), outline)
end

function FT_Outline_Get_Bitmap(library, outline, abitmap)
    ccall((:FT_Outline_Get_Bitmap, libfreetype), FT_Error, (FT_Library, Ptr{FT_Outline}, Ptr{FT_Bitmap}), library, outline, abitmap)
end

function FT_Outline_Render(library, outline, params)
    ccall((:FT_Outline_Render, libfreetype), FT_Error, (FT_Library, Ptr{FT_Outline}, Ptr{FT_Raster_Params}), library, outline, params)
end

@cenum FT_Orientation_::UInt32 begin
    FT_ORIENTATION_TRUETYPE = 0
    FT_ORIENTATION_POSTSCRIPT = 1
    FT_ORIENTATION_FILL_RIGHT = 0
    FT_ORIENTATION_FILL_LEFT = 1
    FT_ORIENTATION_NONE = 2
end

const FT_Orientation = FT_Orientation_

function FT_Outline_Get_Orientation(outline)
    ccall((:FT_Outline_Get_Orientation, libfreetype), FT_Orientation, (Ptr{FT_Outline},), outline)
end

const FT_CONFIG_OPTION_ENVIRONMENT_PROPERTIES = nothing

const FT_CONFIG_OPTION_INLINE_MULFIX = nothing

const FT_CONFIG_OPTION_USE_LZW = nothing

const FT_CONFIG_OPTION_USE_ZLIB = nothing

const FT_CONFIG_OPTION_SYSTEM_ZLIB = nothing

const FT_CONFIG_OPTION_USE_BZIP2 = nothing

const FT_CONFIG_OPTION_POSTSCRIPT_NAMES = nothing

const FT_CONFIG_OPTION_ADOBE_GLYPH_LIST = nothing

const FT_CONFIG_OPTION_MAC_FONTS = nothing

const FT_CONFIG_OPTION_GUESSING_EMBEDDED_RFORK = nothing

const FT_CONFIG_OPTION_INCREMENTAL = nothing

const FT_RENDER_POOL_SIZE = Clong(16384)

const FT_MAX_MODULES = 32

const TT_CONFIG_OPTION_EMBEDDED_BITMAPS = nothing

const TT_CONFIG_OPTION_COLOR_LAYERS = nothing

const TT_CONFIG_OPTION_POSTSCRIPT_NAMES = nothing

const TT_CONFIG_OPTION_SFNT_NAMES = nothing

const TT_CONFIG_CMAP_FORMAT_0 = nothing

const TT_CONFIG_CMAP_FORMAT_2 = nothing

const TT_CONFIG_CMAP_FORMAT_4 = nothing

const TT_CONFIG_CMAP_FORMAT_6 = nothing

const TT_CONFIG_CMAP_FORMAT_8 = nothing

const TT_CONFIG_CMAP_FORMAT_10 = nothing

const TT_CONFIG_CMAP_FORMAT_12 = nothing

const TT_CONFIG_CMAP_FORMAT_13 = nothing

const TT_CONFIG_CMAP_FORMAT_14 = nothing

const TT_CONFIG_OPTION_BYTECODE_INTERPRETER = nothing

const TT_CONFIG_OPTION_SUBPIXEL_HINTING = 2

const TT_CONFIG_OPTION_GX_VAR_SUPPORT = nothing

const TT_CONFIG_OPTION_BDF = nothing

const TT_CONFIG_OPTION_MAX_RUNNABLE_OPCODES = Clong(1000000)

const T1_MAX_DICT_DEPTH = 5

const T1_MAX_SUBRS_CALLS = 16

const T1_MAX_CHARSTRINGS_OPERANDS = 256

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_X1 = 500

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_Y1 = 400

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_X2 = 1000

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_Y2 = 275

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_X3 = 1667

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_Y3 = 275

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_X4 = 2333

const CFF_CONFIG_OPTION_DARKENING_PARAMETER_Y4 = 0

const AF_CONFIG_OPTION_CJK = nothing

const AF_CONFIG_OPTION_INDIC = nothing

const AF_CONFIG_OPTION_USE_WARPER = nothing

const TT_USE_BYTECODE_INTERPRETER = nothing

const TT_SUPPORT_SUBPIXEL_HINTING_MINIMAL = nothing

const ft_ptrdiff_t = ptrdiff_t

const FT_CHAR_BIT = CHAR_BIT

const FT_USHORT_MAX = USHRT_MAX

const FT_INT_MAX = INT_MAX

const FT_INT_MIN = INT_MIN

const FT_UINT_MAX = UINT_MAX

const FT_LONG_MIN = LONG_MIN

const FT_LONG_MAX = LONG_MAX

const FT_ULONG_MAX = ULONG_MAX

const FT_FILE = FILE

const FT_SIZEOF_INT = 32 ÷ FT_CHAR_BIT

const FT_SIZEOF_LONG = 64 ÷ FT_CHAR_BIT

const FT_LONG64 = nothing

const FT_INT64 = Clong

const FT_UINT64 = Culong

const FT_CALLBACK_TABLE_DEF = nothing

const ft_pixel_mode_none = FT_PIXEL_MODE_NONE

const ft_pixel_mode_mono = FT_PIXEL_MODE_MONO

const ft_pixel_mode_grays = FT_PIXEL_MODE_GRAY

const ft_pixel_mode_pal2 = FT_PIXEL_MODE_GRAY2

const ft_pixel_mode_pal4 = FT_PIXEL_MODE_GRAY4

const FT_OUTLINE_NONE = 0x00

const FT_OUTLINE_OWNER = 0x01

const FT_OUTLINE_EVEN_ODD_FILL = 0x02

const FT_OUTLINE_REVERSE_FILL = 0x04

const FT_OUTLINE_IGNORE_DROPOUTS = 0x08

const FT_OUTLINE_SMART_DROPOUTS = 0x10

const FT_OUTLINE_INCLUDE_STUBS = 0x20

const FT_OUTLINE_HIGH_PRECISION = 0x0100

const FT_OUTLINE_SINGLE_PASS = 0x0200

const ft_outline_none = FT_OUTLINE_NONE

const ft_outline_owner = FT_OUTLINE_OWNER

const ft_outline_even_odd_fill = FT_OUTLINE_EVEN_ODD_FILL

const ft_outline_reverse_fill = FT_OUTLINE_REVERSE_FILL

const ft_outline_ignore_dropouts = FT_OUTLINE_IGNORE_DROPOUTS

const ft_outline_high_precision = FT_OUTLINE_HIGH_PRECISION

const ft_outline_single_pass = FT_OUTLINE_SINGLE_PASS

const FT_CURVE_TAG_ON = 0x01

const FT_CURVE_TAG_CONIC = 0x00

const FT_CURVE_TAG_CUBIC = 0x02

const FT_CURVE_TAG_HAS_SCANMODE = 0x04

const FT_CURVE_TAG_TOUCH_X = 0x08

const FT_CURVE_TAG_TOUCH_Y = 0x10

const FT_CURVE_TAG_TOUCH_BOTH = FT_CURVE_TAG_TOUCH_X | FT_CURVE_TAG_TOUCH_Y

const FT_Curve_Tag_On = FT_CURVE_TAG_ON

const FT_Curve_Tag_Conic = FT_CURVE_TAG_CONIC

const FT_Curve_Tag_Cubic = FT_CURVE_TAG_CUBIC

const FT_Curve_Tag_Touch_X = FT_CURVE_TAG_TOUCH_X

const FT_Curve_Tag_Touch_Y = FT_CURVE_TAG_TOUCH_Y

const FT_Outline_MoveTo_Func = FT_Outline_MoveToFunc

const FT_Outline_LineTo_Func = FT_Outline_LineToFunc

const FT_Outline_ConicTo_Func = FT_Outline_ConicToFunc

const FT_Outline_CubicTo_Func = FT_Outline_CubicToFunc

const ft_glyph_format_none = FT_GLYPH_FORMAT_NONE

const ft_glyph_format_composite = FT_GLYPH_FORMAT_COMPOSITE

const ft_glyph_format_bitmap = FT_GLYPH_FORMAT_BITMAP

const ft_glyph_format_outline = FT_GLYPH_FORMAT_OUTLINE

const ft_glyph_format_plotter = FT_GLYPH_FORMAT_PLOTTER

const FT_Raster_Span_Func = FT_SpanFunc

const FT_RASTER_FLAG_DEFAULT = 0x00

const FT_RASTER_FLAG_AA = 0x01

const FT_RASTER_FLAG_DIRECT = 0x02

const FT_RASTER_FLAG_CLIP = 0x04

const ft_raster_flag_default = FT_RASTER_FLAG_DEFAULT

const ft_raster_flag_aa = FT_RASTER_FLAG_AA

const ft_raster_flag_direct = FT_RASTER_FLAG_DIRECT

const ft_raster_flag_clip = FT_RASTER_FLAG_CLIP

const FT_Raster_New_Func = FT_Raster_NewFunc

const FT_Raster_Done_Func = FT_Raster_DoneFunc

const FT_Raster_Reset_Func = FT_Raster_ResetFunc

const FT_Raster_Set_Mode_Func = FT_Raster_SetModeFunc

const FT_Raster_Render_Func = FT_Raster_RenderFunc

const FT_ERR_BASE = 0

const FT_INCLUDE_ERR_PROTOS = nothing

const FT_ERR_PROTOS_DEFINED = nothing

const ft_encoding_none = FT_ENCODING_NONE

const ft_encoding_unicode = FT_ENCODING_UNICODE

const ft_encoding_symbol = FT_ENCODING_MS_SYMBOL

const ft_encoding_latin_1 = FT_ENCODING_ADOBE_LATIN_1

const ft_encoding_latin_2 = FT_ENCODING_OLD_LATIN_2

const ft_encoding_sjis = FT_ENCODING_SJIS

const ft_encoding_gb2312 = FT_ENCODING_PRC

const ft_encoding_big5 = FT_ENCODING_BIG5

const ft_encoding_wansung = FT_ENCODING_WANSUNG

const ft_encoding_johab = FT_ENCODING_JOHAB

const ft_encoding_adobe_standard = FT_ENCODING_ADOBE_STANDARD

const ft_encoding_adobe_expert = FT_ENCODING_ADOBE_EXPERT

const ft_encoding_adobe_custom = FT_ENCODING_ADOBE_CUSTOM

const ft_encoding_apple_roman = FT_ENCODING_APPLE_ROMAN

const FT_FACE_FLAG_SCALABLE = Clong(1) << 0

const FT_FACE_FLAG_FIXED_SIZES = Clong(1) << 1

const FT_FACE_FLAG_FIXED_WIDTH = Clong(1) << 2

const FT_FACE_FLAG_SFNT = Clong(1) << 3

const FT_FACE_FLAG_HORIZONTAL = Clong(1) << 4

const FT_FACE_FLAG_VERTICAL = Clong(1) << 5

const FT_FACE_FLAG_KERNING = Clong(1) << 6

const FT_FACE_FLAG_FAST_GLYPHS = Clong(1) << 7

const FT_FACE_FLAG_MULTIPLE_MASTERS = Clong(1) << 8

const FT_FACE_FLAG_GLYPH_NAMES = Clong(1) << 9

const FT_FACE_FLAG_EXTERNAL_STREAM = Clong(1) << 10

const FT_FACE_FLAG_HINTER = Clong(1) << 11

const FT_FACE_FLAG_CID_KEYED = Clong(1) << 12

const FT_FACE_FLAG_TRICKY = Clong(1) << 13

const FT_FACE_FLAG_COLOR = Clong(1) << 14

const FT_FACE_FLAG_VARIATION = Clong(1) << 15

const FT_STYLE_FLAG_ITALIC = 1 << 0

const FT_STYLE_FLAG_BOLD = 1 << 1

const FT_OPEN_MEMORY = 0x01

const FT_OPEN_STREAM = 0x02

const FT_OPEN_PATHNAME = 0x04

const FT_OPEN_DRIVER = 0x08

const FT_OPEN_PARAMS = 0x10

const ft_open_memory = FT_OPEN_MEMORY

const ft_open_stream = FT_OPEN_STREAM

const ft_open_pathname = FT_OPEN_PATHNAME

const ft_open_driver = FT_OPEN_DRIVER

const ft_open_params = FT_OPEN_PARAMS

const FT_LOAD_DEFAULT = 0x00

const FT_LOAD_NO_SCALE = Clong(1) << 0

const FT_LOAD_NO_HINTING = Clong(1) << 1

const FT_LOAD_RENDER = Clong(1) << 2

const FT_LOAD_NO_BITMAP = Clong(1) << 3

const FT_LOAD_VERTICAL_LAYOUT = Clong(1) << 4

const FT_LOAD_FORCE_AUTOHINT = Clong(1) << 5

const FT_LOAD_CROP_BITMAP = Clong(1) << 6

const FT_LOAD_PEDANTIC = Clong(1) << 7

const FT_LOAD_IGNORE_GLOBAL_ADVANCE_WIDTH = Clong(1) << 9

const FT_LOAD_NO_RECURSE = Clong(1) << 10

const FT_LOAD_IGNORE_TRANSFORM = Clong(1) << 11

const FT_LOAD_MONOCHROME = Clong(1) << 12

const FT_LOAD_LINEAR_DESIGN = Clong(1) << 13

const FT_LOAD_NO_AUTOHINT = Clong(1) << 15

const FT_LOAD_COLOR = Clong(1) << 20

const FT_LOAD_COMPUTE_METRICS = Clong(1) << 21

const FT_LOAD_BITMAP_METRICS_ONLY = Clong(1) << 22

const FT_LOAD_ADVANCE_ONLY = Clong(1) << 8

const FT_LOAD_SBITS_ONLY = Clong(1) << 14

FT_LOAD_TARGET_(x) = FT_Int32(x & 15) << 16

const FT_LOAD_TARGET_NORMAL = FT_LOAD_TARGET_(FT_RENDER_MODE_NORMAL)

const FT_LOAD_TARGET_LIGHT = FT_LOAD_TARGET_(FT_RENDER_MODE_LIGHT)

const FT_LOAD_TARGET_MONO = FT_LOAD_TARGET_(FT_RENDER_MODE_MONO)

const FT_LOAD_TARGET_LCD = FT_LOAD_TARGET_(FT_RENDER_MODE_LCD)

const FT_LOAD_TARGET_LCD_V = FT_LOAD_TARGET_(FT_RENDER_MODE_LCD_V)

const ft_render_mode_normal = FT_RENDER_MODE_NORMAL

const ft_render_mode_mono = FT_RENDER_MODE_MONO

const ft_kerning_default = FT_KERNING_DEFAULT

const ft_kerning_unfitted = FT_KERNING_UNFITTED

const ft_kerning_unscaled = FT_KERNING_UNSCALED

const FT_SUBGLYPH_FLAG_ARGS_ARE_WORDS = 1

const FT_SUBGLYPH_FLAG_ARGS_ARE_XY_VALUES = 2

const FT_SUBGLYPH_FLAG_ROUND_XY_TO_GRID = 4

const FT_SUBGLYPH_FLAG_SCALE = 8

const FT_SUBGLYPH_FLAG_XY_SCALE = 0x40

const FT_SUBGLYPH_FLAG_2X2 = 0x80

const FT_SUBGLYPH_FLAG_USE_MY_METRICS = 0x0200

const FT_FSTYPE_INSTALLABLE_EMBEDDING = 0x0000

const FT_FSTYPE_RESTRICTED_LICENSE_EMBEDDING = 0x0002

const FT_FSTYPE_PREVIEW_AND_PRINT_EMBEDDING = 0x0004

const FT_FSTYPE_EDITABLE_EMBEDDING = 0x0008

const FT_FSTYPE_NO_SUBSETTING = 0x0100

const FT_FSTYPE_BITMAP_EMBEDDING_ONLY = 0x0200

const FREETYPE_MAJOR = 2

const FREETYPE_MINOR = 10

const FREETYPE_PATCH = 1

# exports
const PREFIXES = ["FREETYPE_", "FT_", "ft_", "TT_", "TTAG_", "CFF_", "T1_", "CID_", "PS_", "t1_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module

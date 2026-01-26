check_error(err, error_msg) = err == 0 || error("$error_msg with error: $err")

const LIBRARY_LOCK = ReentrantLock()
const FREE_FONT_LIBRARY = FT_Library[C_NULL]

function ft_init()
    @lock LIBRARY_LOCK begin
        if FREE_FONT_LIBRARY[1] != C_NULL
            error("Freetype already initalized. init() called two times?")
        end
        return FT_Init_FreeType(FREE_FONT_LIBRARY) == 0
    end
end

function ft_done()
    @lock LIBRARY_LOCK begin
        if FREE_FONT_LIBRARY[1] == C_NULL
            error("Library == CNULL. FreeTypeAbstraction.done() called before init(), or done called two times?")
        end
        err = FT_Done_FreeType(FREE_FONT_LIBRARY[1])
        FREE_FONT_LIBRARY[1] = C_NULL
        return err == 0
    end
end

function newface(facename, faceindex::Real=0, ftlib=FREE_FONT_LIBRARY)
    face = Ref{FT_Face}()
    err = @lock LIBRARY_LOCK FT_New_Face(ftlib[1], facename, Int32(faceindex), face)
    check_error(err, "Couldn't load font $facename")
    return face[]
end

function newface_mmapped(filepath, faceindex::Real=0, ftlib=FREE_FONT_LIBRARY)
    mmapped = open(filepath, "r") do io
        Mmap.mmap(io)
    end
    face = Ref{FT_Face}()
    err = @lock LIBRARY_LOCK FT_New_Memory_Face(ftlib[1], mmapped, length(mmapped), Int32(faceindex), face)
    check_error(err, "Couldn't load font \"$filepath\"")
    return face[], mmapped
end


struct FontExtent{T}
    vertical_bearing::Vec{2, T}
    horizontal_bearing::Vec{2, T}

    advance::Vec{2, T}
    scale::Vec{2, T}
end

hadvance(ext::FontExtent) = ext.advance[1]
vadvance(ext::FontExtent) = ext.advance[2]
inkwidth(ext::FontExtent) = ext.scale[1]
inkheight(ext::FontExtent) = ext.scale[2]
hbearing_ori_to_left(ext::FontExtent) = ext.horizontal_bearing[1]
hbearing_ori_to_top(ext::FontExtent) = ext.horizontal_bearing[2]
leftinkbound(ext::FontExtent) = hbearing_ori_to_left(ext)
rightinkbound(ext::FontExtent) = leftinkbound(ext) + inkwidth(ext)
bottominkbound(ext::FontExtent) = hbearing_ori_to_top(ext) - inkheight(ext)
topinkbound(ext::FontExtent) = hbearing_ori_to_top(ext)

BroadcastStyle(::Type{<: FontExtent}) = Style{FontExtent}()
BroadcastStyle(::Style{FontExtent}, x) = Style{FontExtent}()
BroadcastStyle(x, ::Style{FontExtent}) = Style{FontExtent}()

function broadcasted(op, f::FontExtent, scaling::Vec)
    return FontExtent(
        op.(f.vertical_bearing, scaling[1]),
        op.(f.horizontal_bearing, scaling[2]),
        op.(f.advance, scaling),
        op.(f.scale, scaling),
    )
end

function broadcasted(op, f::FontExtent)
    return FontExtent(
        op.(f.vertical_bearing),
        op.(f.horizontal_bearing),
        op.(f.advance),
        op.(f.scale),
    )
end

function broadcasted(op, ::Type{T}, f::FontExtent) where T
    return FontExtent(
        map(x-> op(T, x), f.vertical_bearing),
        map(x-> op(T, x), f.horizontal_bearing),
        map(x-> op(T, x), f.advance),
        map(x-> op(T, x), f.scale),
    )
end

function FontExtent(fontmetric::FreeType.FT_Glyph_Metrics, scale::T = 64.0) where T <: AbstractFloat
    return FontExtent(
        Vec{2, T}(fontmetric.vertBearingX, fontmetric.vertBearingY) ./ scale,
        Vec{2, T}(fontmetric.horiBearingX, fontmetric.horiBearingY) ./ scale,
        Vec{2, T}(fontmetric.horiAdvance, fontmetric.vertAdvance) ./ scale,
        Vec{2, T}(fontmetric.width, fontmetric.height) ./ scale
    )
end

function ==(x::FontExtent, y::FontExtent)
    return (
        x.vertical_bearing == y.vertical_bearing &&
        x.horizontal_bearing == y.horizontal_bearing &&
        x.advance == y.advance &&
        x.scale == y.scale
    )
end

function FontExtent(fontmetric::FreeType.FT_Glyph_Metrics, scale::Integer)
    return FontExtent(
        div.(Vec{2, Int}(fontmetric.vertBearingX, fontmetric.vertBearingY), scale),
        div.(Vec{2, Int}(fontmetric.horiBearingX, fontmetric.horiBearingY), scale),
        div.(Vec{2, Int}(fontmetric.horiAdvance, fontmetric.vertAdvance), scale),
        div.(Vec{2, Int}(fontmetric.width, fontmetric.height), scale)
    )
end

function bearing(extent::FontExtent{T}) where T
    # origin to SW corner of the horizontal metric
    # with the conventions of freetype.org/freetype2/docs/glyphs/glyphs-3.html
    return Vec2{T}(
        hbearing_ori_to_left(extent),
        hbearing_ori_to_top(extent) - inkheight(extent),
    )
end

function safe_free(face)
    @lock face.lock begin
        ptr = getfield(face, :ft_ptr)
        if ptr != C_NULL && FREE_FONT_LIBRARY[1] != C_NULL
            FT_Done_Face(face)
        end
    end
end

function boundingbox(extent::FontExtent{T}) where T
    return Rect2(bearing(extent), Vec2{T}(extent.scale))
end

mutable struct FTFont
    ft_ptr::FreeType.FT_Face
    use_cache::Bool
    extent_cache::Dict{UInt64, FontExtent{Float32}}
    lock::ReentrantLock # lock this for the duration of any FT operation on ft_ptr
    mmapped::Union{Nothing,Vector{UInt8}}
    fontname::String
    function FTFont(ft_ptr::FreeType.FT_Face, use_cache::Bool=true, mmapped = nothing)
        extent_cache = Dict{UInt64, FontExtent{Float32}}()
        face = new(ft_ptr, use_cache, extent_cache, ReentrantLock(), mmapped, "")
        face.fontname = "$(family_name(face)) $(style_name(face))"
        finalizer(safe_free, face)
        return face
    end
end

use_cache(face::FTFont) = getfield(face, :use_cache)
get_cache(face::FTFont) = getfield(face, :extent_cache)

function FTFont(path::String, use_cache::Bool=true)
    face, mmapped = newface_mmapped(path)
    FTFont(face, use_cache, mmapped)
end

# C interop
Base.cconvert(::Type{FreeType.FT_Face}, font::FTFont) = font

function Base.unsafe_convert(::Type{FreeType.FT_Face}, font::FTFont)
    return getfield(font, :ft_ptr)
end

Base.propertynames(font::FTFont) = fieldnames(FreeType.FT_FaceRec)

function Base.getproperty(font::FTFont, fieldname::Symbol)
    fieldname in fieldnames(FTFont) && return getfield(font, fieldname)
    @lock font.lock begin
        fontrect = unsafe_load(getfield(font, :ft_ptr))
        field = getfield(fontrect, fieldname)
        if field isa Ptr{FT_String}
            field == C_NULL && return ""
            return unsafe_string(field)
        else
            return field
        end
    end
end

function Base.show(io::IO, font::FTFont)
    print(io, "FTFont (family = $(font.family_name), style = $(font.style_name))")
end

# Allow broadcasting over fonts
Base.Broadcast.broadcastable(ft::FTFont) = Ref(ft)

function set_pixelsize(face::FTFont, size::Integer)
    @lock face.lock begin
        err = FT_Set_Pixel_Sizes(face, size, size)
        check_error(err, "Couldn't set pixelsize")
        return size
    end
end

function kerning(glyphspec1, glyphspec2, face::FTFont)
    i1 = glyph_index(face, glyphspec1)
    i2 = glyph_index(face, glyphspec2)
    kerning2d = Ref{FreeType.FT_Vector}()
    err = @lock face.lock begin
        FT_Get_Kerning(face, i1, i2, FreeType.FT_KERNING_DEFAULT, kerning2d)
    end
    # Can error if font has no kerning! Since that's somewhat expected, we just return 0
    err != 0 && return Vec2f(0)
    # 64 since metrics are in 1/64 units (units to 26.6 fractional pixels)
    divisor = 64
    return Vec2f(kerning2d[].x / divisor, kerning2d[].y / divisor)
end

function get_extent(face::FTFont, glyphspec)
    gi = glyph_index(face, glyphspec)
    if use_cache(face)
        lock(face.lock) do
            get!(get_cache(face), gi) do
                return internal_get_extent(face, gi)
            end
        end
    else
        return internal_get_extent(face, gi)
    end
end

glyph_index(face::FTFont, glyphname::String)::UInt64 = @lock face.lock FT_Get_Name_Index(face, glyphname)
glyph_index(face::FTFont, char::Char)::UInt64 = @lock face.lock FT_Get_Char_Index(face, char)
glyph_index(face::FTFont, int::Integer) = UInt64(int)

function internal_get_extent(face::FTFont, glyphspec)
    gi = glyph_index(face, glyphspec)
    #=
    Load chars without scaling. This leaves all glyph metrics that can be
    retrieved in font units, which can be normalized by dividing with the
    font's units_per_EM. This is more robust than relying on extents
    that are only valid with a specific pixelsize, because a font's
    pixelsize can be silently changed by third parties, such as Cairo.
    If that happens, all glyph metrics are incorrect. We avoid this by using the normalized space.
    =#
    err = @lock face.lock FT_Load_Glyph(face, gi, FT_LOAD_NO_SCALE)
    check_error(err, "Could not load glyph $(repr(glyphspec)) from $(face) to get extent.")
    # This gives us the font metrics in normalized units (0, 1), with negative
    # numbers interpreted as an offset
    metrics = @lock face.lock unsafe_load(face.glyph).metrics
    return FontExtent(metrics, Float32(face.units_per_EM))
end

descender(font) = font.descender / font.units_per_EM
ascender(font) = font.ascender / font.units_per_EM

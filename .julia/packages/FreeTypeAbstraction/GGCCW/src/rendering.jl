
function load_glyph(face::FTFont, glyph)
    gi = glyph_index(face, glyph)
    err = @lock face.lock FT_Load_Glyph(face, gi, FT_LOAD_RENDER)
    check_error(err, "Could not load glyph $(repr(glyph)) from $(face) to render.")
end

function loadglyph(face::FTFont, glyph, pixelsize::Integer)
    set_pixelsize(face, pixelsize)
    load_glyph(face, glyph)
    gl = @lock face.lock unsafe_load(face.glyph)
    @assert gl.format == FreeType.FT_GLYPH_FORMAT_BITMAP
    return gl
end

function renderface(face::FTFont, glyph, pixelsize::Integer)
    gl = loadglyph(face, glyph, pixelsize)
    return glyphbitmap(gl.bitmap), FontExtent(gl.metrics)
end

function extents(face::FTFont, glyph, pixelsize::Integer)
    return FontExtent(loadglyph(face, glyph, pixelsize).metrics)
end

function glyphbitmap(bitmap::FreeType.FT_Bitmap)
    @assert bitmap.pixel_mode == FreeType.FT_PIXEL_MODE_GRAY
    bmp = Matrix{UInt8}(undef, bitmap.width, bitmap.rows)
    row = bitmap.buffer
    if bitmap.pitch < 0
        row -= bitmap.pitch * (rbmpRec.rows - 1)
    end
    for r in 1:bitmap.rows
        src = unsafe_wrap(Array, row, bitmap.width)
        bmp[:, r] = src
        row += bitmap.pitch
    end
    return bmp
end

function one_or_typemax(::Type{T}) where {T<:Union{Real,Colorant}}
    return T<:Integer ? typemax(T) : oneunit(T)
end

"""
    renderstring!(img::AbstractMatrix, str::String, face, pixelsize, y0, x0;
    fcolor=one_or_typemax(T), bcolor=zero(T), halign=:hleft, valign=:vbaseline) -> Matrix

Render `str` into `img` using the font `face` of size `pixelsize` at coordinates `y0,x0`.
Uses the conventions of freetype.org/freetype2/docs/glyphs/glyphs-3.html

# Arguments
* `y0,x0`: origin is in upper left with positive `y` going down
* `fcolor`: foreground color; AbstractVector{T}, typemax(T) for T<:Integer, otherwise one(T)
* `gcolor`: background color; AbstractVector{T}, typemax(T) for T<:Integer, otherwise one(T)
* `bcolor`: canvas background color; set to `nothing` for transparent
* `halign`: :hleft, :hcenter, or :hright
* `valign`: :vtop, :vcenter, :vbaseline, or :vbottom
* `bbox_glyph`: glyph bounding box color (debugging)
* `bbox`: bounding box color (debugging)
* `gstr`: background string or array of chars (for background sizing)
* `incx`: extra x spacing
"""
function renderstring!(
        img::AbstractMatrix{T}, fstr::Union{AbstractVector{Char},String},
        face::FTFont, pixelsize::Union{Int, Tuple{Int, Int}}, y0, x0;
        fcolor::Union{AbstractVector{T},T} = one_or_typemax(T),
        gcolor::Union{AbstractVector{T},T,Nothing} = nothing,
        bcolor::Union{T,Nothing} = zero(T),
        halign::Symbol = :hleft,
        valign::Symbol = :vbaseline,
        bbox_glyph::Union{T,Nothing} = nothing,
        bbox::Union{T,Nothing} = nothing,
        gstr::Union{AbstractVector{Char},String,Nothing} = nothing,
        off_bg::Int = 0,
        incx::Int = 0,
    ) where T<:Union{Real,Colorant}

    if pixelsize isa Tuple
        @warn "using tuple for pixelsize is deprecated, please use one integer"
        pixelsize = pixelsize[1]
    end
    set_pixelsize(face, pixelsize)

    fstr = fstr isa AbstractVector ? fstr : collect(fstr)
    if gstr !== nothing
        gstr = gstr isa AbstractVector ? gstr : collect(gstr)
    end

    len = length(fstr)
    bitmaps = Vector{Matrix{UInt8}}(undef, len)
    metrics = Vector{FontExtent{Int}}(undef, len)

    y_min = y_max = sum_advance_x = 0  # y_min and y_max are w.r.t the baseline
    for (i, char) in enumerate(fstr)
        bitmap, metricf = renderface(face, char, pixelsize)
        metric = round.(Int, metricf)
        bitmaps[i] = bitmap
        metrics[i] = metric

        y_min = min(y_min, bottominkbound(metric))
        y_max = max(y_max, topinkbound(metric))
        sum_advance_x += hadvance(metric)
    end

    bitmap_max = bitmaps |> first |> eltype |> typemax
    imgh, imgw = size(img)

    # initial pen position
    px = x0 - (halign == :hright ? sum_advance_x : halign == :hcenter ? sum_advance_x >> 1 : 0)
    py = y0 + (
        valign == :vtop ? y_max : valign == :vbottom ? y_min :
        valign == :vcenter ? (y_max - y_min) >> 1 + y_min : 0
    )

    if bcolor !== nothing
        img[
            clamp(py - y_max, 1, imgh) : clamp(py - y_min, 1, imgh),
            clamp(px, 1, imgw) : clamp(px + sum_advance_x, 1, imgw)
        ] .= bcolor
    end

    local prev_char::Char
    for (i, char) in enumerate(fstr)
        bitmap = bitmaps[i]
        metric = metrics[i]
        bx, by = metric.horizontal_bearing
        ax, ay = metric.advance
        sx, sy = metric.scale

        if i == 1
            prev_char = char
        else
            kx, _ = map(x-> round(Int, x), kerning(prev_char, char, face))
            px += kx
        end

        # glyph origin
        oy = py - by
        ox = px + bx

        fcol = fcolor isa AbstractVector ? fcolor[i] : fcolor
        gcol = gcolor isa AbstractVector ? gcolor[i] : gcolor

        # trim parts of glyph images that are outside the destination
        row_lo, row_hi = 1 + max(0, -oy), sy - max(0, oy + sy - imgh)
        col_lo, col_hi = 1 + max(0, -ox), sx - max(0, ox + sx - imgw)

        if gcol === nothing
            for r in row_lo:row_hi, c in col_lo:col_hi
                (bm = bitmap[c, r]) == 0 && continue
                color = bm / bitmap_max * fcol
                img[oy + r, ox + c] = T <: Integer ? round(T, color) : T(color)
            end
        else
            if gstr !== nothing
                gmetric = round.(Int, extents(face, gstr[i], pixelsize))
                y_min = bottominkbound(gmetric)
                y_max = topinkbound(gmetric)
            end

            # fill background
            by1, by2 = py - y_max, py - y_min
            bx1, bx2 = px, px + ax
            r1, r2 = clamp(by1, 1, imgh), clamp(by2, 1, imgh)
            c1, c2 = clamp(bx1, 1, imgw), clamp(bx2, 1, imgw)
            for r in r1 + off_bg:r2 - off_bg, c in c1 + off_bg:c2 - off_bg
                img[r, c] = gcol
            end

            # render character by drawing the corresponding glyph
            for r in row_lo:row_hi, c in col_lo:col_hi
                (bm = bitmap[c, r]) == 0 && continue
                w1 = bm / bitmap_max
                color0 = w1 * fcol
                color1 = (1.0 - w1) * gcol
                img[oy + r, ox + c] = T <: Integer ? round(T, color0 + color1) : T(color0 + color1)
            end

            # draw background bounding box
            if bbox !== nothing && r2 > r1 && c2 > c1
                img[r1, c1:c2] .= bbox
                img[r2, c1:c2] .= bbox
                img[r1:r2, c1] .= bbox
                img[r1:r2, c2] .= bbox
            end
        end

        # draw glyph bounding box
        if bbox_glyph !== nothing
            # rect = boundingbox(metric)
            # (mc, mr), (Mc, Mr) = extrema(rect)
            # r1, r2 = clamp(py + mr, 1, imgh), clamp(py + Mr, 1, imgh)
            # c1, c2 = clamp(px + mc, 1, imgw), clamp(px + Mc, 1, imgw)
            # ^^^ should be equivalent to vvv: see JuliaGraphics/FreeTypeAbstraction.jl/pull/71
            r1, r2 = clamp(oy + row_lo, 1, imgh), clamp(oy + row_hi, 1, imgh)
            c1, c2 = clamp(ox + col_lo, 1, imgw), clamp(ox + col_hi, 1, imgw)
            if r2 > r1 && c2 > c1
                img[r1, c1:c2] .= bbox_glyph
                img[r2, c1:c2] .= bbox_glyph
                img[r1:r2, c1] .= bbox_glyph
                img[r1:r2, c2] .= bbox_glyph
            end
        end

        px += ax + incx
    end
    return img
end

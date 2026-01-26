module Graphics

using Colors
using NaNMath

import Base: +, -, *, /, &, fill
import LinearAlgebra: norm

"""
Graphics defines an API for drawing in two dimensions.

- Geometric primitives: `Vec2`, `Point`, `BoundingBox`

- Geometry API: `aspect_ratio`, `center`, `deform`, `diagonal`,
  `isinside`, `shift`, `height`, `width`, `xmin`, `xmax`, `ymin`,
  `ymax`, `xrange`, `yrange`

- 2d drawing contexts: `GraphicsDevice`, `GraphicsContext`, `creategc`, `getgc`

- Coordinate systems: `set_coordinates`, `reset_transform`, `rotate`,
  `scale`, `translate`, `user_to_device!`, `device_to_user!`,
  `user_to_device_distance!`, `device_to_user_distance!`,
  `user_to_device`, `device_to_user`

- Lines: `set_line_width`, `set_dash`

- Colors and painting (drawing attributes): `set_source`,
  `set_source_rgb`, `set_source_rgba`, `save`, `restore`

- Clipping: `clip`, `clip_preserve`, `reset_clip`, `inner_canvas`

- Paths: `move_to`, `line_to`, `rel_line_to`, `rel_move_to`,
  `new_path`, `new_sub_path`, `close_path`, `arc`

- High-level paths: `rectangle`, `circle`, `polygon`

- Fill and stroke: `fill`, `fill_preserve`, `paint`, `stroke`,
  `stroke_preserve`, `stroke_transformed`,
  `stroke_transformed_preserve`
"""
Graphics

export
    # Part 1. 2D Geometry
    Vec2, Point, BoundingBox,
    # limits in world coordinates
    isinside, xmin, xmax, ymin, ymax, center, xrange, yrange,
    aspect_ratio, with_aspect_ratio, diagonal, shift, deform,
    # TODO: more

    # TODO: 3D geometry

    # Part 2. 2D Drawing
    # device and context
    GraphicsDevice, GraphicsContext, creategc, getgc,
    # width, height are in user (world) coordinates for geometric objects,
    # but in device coordinates for GraphicsDevice, GraphicsContext, and
    # other concrete things like windows and widgets.
    width, height,

    # drawing attribute manipulation
    save, restore, set_line_width, set_dash, set_source_rgb, set_source_rgba,
    set_source,

    # coordinate systems
    reset_transform, set_coordinates, rotate, scale, translate, user_to_device!,
    device_to_user!, user_to_device_distance!, device_to_user_distance!,
    user_to_device, device_to_user,

    # clipping
    clip, clip_preserve, reset_clip, inner_canvas,

    # path primitives
    move_to, line_to, rel_line_to, rel_move_to, new_path, new_sub_path,
    close_path, arc,

    # fill and stroke
    fill, fill_preserve, paint, stroke, stroke_preserve,
    stroke_transformed, stroke_transformed_preserve,

    # derived path operations
    rectangle, circle, polygon

    # TODO: text drawing API

    # TODO: rendering pipeline API


# Part 1. geometric primitives

"""
    Vec2(x, y) -> v

Create a Cartesian representation `v` of a vector (or point) in two dimensions.
"""
struct Vec2
    x::Float64
    y::Float64
end

"""
    Point(x, y) -> p

Create a Cartesian representation `p` of a point in two dimensions.
`Point` is an alias of [`Vec2`](@ref).
"""
const Point = Vec2

(+)(a::Vec2, b::Vec2) = Vec2(a.x + b.x, a.y + b.y)
(-)(a::Vec2, b::Vec2) = Vec2(a.x - b.x, a.y - b.y)
(*)(p::Vec2, s::Real) = Vec2(p.x*s, p.y*s)
(/)(p::Vec2, s::Real) = Vec2(p.x/s, p.y/s)
(*)(s::Real, p::Vec2) = p*s

norm(p::Vec2) = hypot(p.x, p.y)

"""
    rotate(p::Vec2, angle::Real, o::Vec2 = Vec2(0, 0)) -> pnew::Vec2

Rotate `p` around `o` by `angle` (in radians).

!!! note
    The direction of rotation for positive angles is from the positive x-axis
    toward the positive y-axis. For example, the direction of rotation is
    "clockwise" when the x-axis is directed right and the y-axis is directed
    downward.

# Example
```jldoctest
julia> rotate(Vec2(2, 1), 0.5π, Vec2(1, 1))
Vec2(1.0, 2.0)
```
"""
function rotate(p::Vec2, angle::Real, o::Vec2 = Vec2(0., 0.))
    c = cos(angle)
    s = sin(angle)
    d = p - o
    Vec2(o.x + c*d.x - s*d.y, o.y + s*d.x + c*d.y)
end

"""
    BoundingBox(xmin, xmax, ymin, ymax) -> bb

Create a representation `bb` of a rectangular region, specifying the
coordinates of the horizontal (x) and vertical (y) edges.
"""
struct BoundingBox
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end

BoundingBox() = BoundingBox(NaN, NaN, NaN, NaN)

"""
    BoundingBox(p0::Point, points::Point...) -> bb

Compute the BoundingBox `bb` that minimally encloses all of the input points.
"""
function BoundingBox(p0::Point, points::Point...)
    xmin, xmax, ymin, ymax = p0.x, p0.x, p0.y, p0.y
    for p in points
        xmin = NaNMath.min(xmin, p.x)
        xmax = NaNMath.max(xmax, p.x)
        ymin = NaNMath.min(ymin, p.y)
        ymax = NaNMath.max(ymax, p.y)
    end
    return BoundingBox(xmin, xmax, ymin, ymax)
end

"""
    BoundingBox(bb0::BoundingBox, bboxes::BoundingBox...) -> bb

Compute the BoundingBox `bb` that minimally encloses all of the input boxes.
"""
function BoundingBox(bb0::BoundingBox, bboxes::BoundingBox...)
    xmin, xmax, ymin, ymax = bb0.xmin, bb0.xmax, bb0.ymin, bb0.ymax
    for bb in bboxes
        xmin = NaNMath.min(xmin, bb.xmin)
        xmax = NaNMath.max(xmax, bb.xmax)
        ymin = NaNMath.min(ymin, bb.ymin)
        ymax = NaNMath.max(ymax, bb.ymax)
    end
    return BoundingBox(xmin, xmax, ymin, ymax)
end

"""
    width(obj) -> w

Get the horizontal length of `obj`.
"""
width(bb::BoundingBox) = bb.xmax - bb.xmin # TODO: generalization with xmax()/xmin()

"""
    height(obj) -> h

Get the vertical length of `obj`.
"""
height(bb::BoundingBox) = bb.ymax - bb.ymin # TODO: generalization with ymax()/ymin()

"""
    diagonal(obj) -> diag

Get the diagonal length of `obj`. The fallback implementation of this function
returns the diagonal length of the bounding box of `obj`.
"""
diagonal(obj) = hypot(width(obj), height(obj))

"""
    aspect_ratio(bb::BoundingBox) -> r

Compute the ratio `r` of the height and width of `bb`.

# Example
```jldoctest
julia> bb = BoundingBox(Point(0, 0), Point(1920, 1080)); # landscape

julia> aspect_ratio(bb)
0.5625

julia> rationalize(ans)
9//16
```

"""
aspect_ratio(bb) = height(bb)/width(bb)

"""
    xmin(obj) -> xmin

Get the minimum x coordinate of the bounding box of `obj`.
"""
function xmin end
"""
    xmax(obj) -> xmax

Get the maximum x coordinate of the bounding box of `obj`.
"""
function xmax end
"""
    ymin(obj) -> ymin

Get the minimum x coordinate of the bounding box of `obj`.
"""
function ymin end
"""
    ymax(obj) -> ymax

Get the maximum y coordinate of the bounding box of `obj`.
"""
function ymax end

xmin(bb::BoundingBox) = bb.xmin
xmax(bb::BoundingBox) = bb.xmax
ymin(bb::BoundingBox) = bb.ymin
ymax(bb::BoundingBox) = bb.ymax

"""
    center(obj) -> p::Point

Compute the center coordinate of `obj`.

!!! note
    The fallback implementation of this function returns the center of bounding
    box, not the geometric center, or centroid.
"""
center(obj) = Point((xmin(obj)+xmax(obj))/2, (ymin(obj)+ymax(obj))/2)

"""
    xrange(obj) -> (xmin, xmax)

Get the horizontal range of the bounding box that minimally contains `obj`.
"""
xrange(obj) = xmin(obj), xmax(obj)

"""
    yrange(obj) -> (ymin, ymax)

Get the vertical range of the bounding box that minimally contains `obj`.
"""
yrange(obj) = ymin(obj), ymax(obj)

"""
    bb1 + bb2 -> bb

Compute the BoundingBox `bb` that minimally contains `bb1` and `bb2`
"""
function (+)(bb1::BoundingBox, bb2::BoundingBox)
    BoundingBox(NaNMath.min(bb1.xmin, bb2.xmin),
                NaNMath.max(bb1.xmax, bb2.xmax),
                NaNMath.min(bb1.ymin, bb2.ymin),
                NaNMath.max(bb1.ymax, bb2.ymax))
end

"""
    bb1 & bb2 -> bb

Compute the intersection of two BoundingBoxes.
"""
function (&)(bb1::BoundingBox, bb2::BoundingBox)
    BoundingBox(NaNMath.max(bb1.xmin, bb2.xmin),
                NaNMath.min(bb1.xmax, bb2.xmax),
                NaNMath.max(bb1.ymin, bb2.ymin),
                NaNMath.min(bb1.ymax, bb2.ymax))
end

"""
    deform(bb::BoundingBox, Δl, Δr, Δt, Δb) -> bbnew

Add `Δl` (left), `Δr` (right), `Δt` (top), and `Δb` (bottom) to the edges of a
`BoundingBox`. The sign of each value follows the positive direction of the
axis. The "top" and "bottom" are representations when the y-axis is directed
downward.

# Example
```jldoctest
julia> bb = BoundingBox(Point(1, 1), Point(10, 10));

julia> bbnew = deform(bb, 0.5, -1.0, -0.25, 1.0);

julia> xrange(bbnew), yrange(bbnew)
((1.5, 9.0), (0.75, 11.0))
```
"""
function deform(bb::BoundingBox, dl, dr, dt, db)
    BoundingBox(bb.xmin + dl, bb.xmax + dr, bb.ymin + dt, bb.ymax + db)
end

"""
    shift(bb::BoundingBox, Δx, Δy) -> bbnew

Shift center by `(Δx, Δy)`, keeping width & height fixed.
"""
function shift(bb::BoundingBox, dx, dy)
    BoundingBox(bb.xmin + dx, bb.xmax + dx, bb.ymin + dy, bb.ymax + dy)
end

"""
    s*bb -> bbnew
    bb*s -> bbnew

Scale width & height of `BoundingBox` `bb` by `s`, keeping center fixed.
"""
function (*)(bb::BoundingBox, s::Real)
    dw = 0.5*(s - 1)*width(bb)
    dh = 0.5*(s - 1)*height(bb)
    deform(bb, -dw, dw, -dh, dh)
end
(*)(s::Real, bb::BoundingBox) = bb*s

"""
    rotate(bb::BoundingBox, angle::Real, o::Point) -> bbnew

Rotate `bb` around `o` by `angle` (in radians), returning the `BoundingBox` that
encloses the vertices of the rotated box.
"""
function rotate(bb::BoundingBox, angle::Real, p::Point)
    a = rotate(Point(bb.xmin,bb.ymin), angle, p)
    b = rotate(Point(bb.xmax,bb.ymin), angle, p)
    c = rotate(Point(bb.xmin,bb.ymax), angle, p)
    d = rotate(Point(bb.xmax,bb.ymax), angle, p)
    BoundingBox(a, b, c, d)
end

function with_aspect_ratio(bb::BoundingBox, ratio::Real)
    if ratio < aspect_ratio(bb)
        dh = height(bb) - ratio * width(bb)
        return BoundingBox(bb.xmin, bb.xmax, bb.ymin + dh/2, bb.ymax - dh/2)
    else
        dw = width(bb) - height(bb) / ratio
        return BoundingBox(bb.xmin + dw/2, bb.xmax - dw/2, bb.ymin, bb.ymax)
    end
end

"""
    isinside(bb::BoundingBox, p::Point) -> tf::Bool
    isinside(bb::BoundingBox, x, y) -> tf::Bool

Determine whether the point lies within `bb`.
"""
isinside(bb::BoundingBox, x, y) = (bb.xmin <= x <= bb.xmax) && (bb.ymin <= y <= bb.ymax)
isinside(bb::BoundingBox, p::Point) = isinside(bb, p.x, p.y)


# Part 2. Drawing

macro mustimplement(sig)
    fname = sig.args[1]
    arg1 = sig.args[2]
    if isa(arg1,Expr)
        arg1 = arg1.args[1]
    end
    :($(esc(sig)) = error(typeof($(esc(arg1))),
                          " must implement ", $(Expr(:quote,fname))))
end

"""
    GraphicDevice

An abstract graphics output device; can create [`GraphicsContext`](@ref)s.
"""
abstract type GraphicsDevice end

@mustimplement width(gd::GraphicsDevice)
@mustimplement height(gd::GraphicsDevice)
"""
    creategc(gd::GraphicsDevice) -> gc::GraphicContext

Create a new `GraphicContext`.
"""
@mustimplement creategc(gd::GraphicsDevice)
xmin(g::GraphicsDevice) = 0 # FIXME: use the same type as `xmax`
xmax(g::GraphicsDevice) = width(g)
ymin(g::GraphicsDevice) = 0 # FIXME: use the same type as `ymax`
ymax(g::GraphicsDevice) = height(g)

"""
    GraphicsContext

An abstract object that can actually be drawn to.
"""
abstract type GraphicsContext end

@mustimplement width(gc::GraphicsContext)
@mustimplement height(gc::GraphicsContext)

"""
    getgc(obj) -> gc::GraphicContext

Get a `GraphicsContext` from something that might be drawable.
"""
getgc(gc::GraphicsContext) = gc


# transformations

"""
    inner_canvas(gc::GraphicsContext, device::BoundingBox, user::BoundingBox)
    inner_canvas(gc::GraphicsContext, x, y, w, h, l, r, t, b)

Create a rectangular drawing area inside `device` (represented in
device-coordinates), giving it user-coordinates `user`. Any drawing
that occurs outside this box is clipped.

`x`, `y`, `w`, and `h` are an alternative parametrization of `device`,
and `l`, `r`, `t`, `b` parametrize `user`.

See also: [`set_coordinates`](@ref).
"""
inner_canvas(gc::GraphicsContext, device::BoundingBox, user::BoundingBox) =
    inner_canvas(gc,
                 device.xmin, device.ymin, width(device), height(device),
                 user.xmin, user.xmax, user.ymin, user.ymax)

function inner_canvas(gc::GraphicsContext, x, y, w, h, l, r, t, b)
    reset_transform(gc)
    reset_clip(gc)
    rectangle(gc, x, y, w, h)
    clip(gc)
    _set_coordinates(gc, x, y, w, h, l, r, t, b)
end

function set_coordinates(gc::GraphicsContext, x, y, w, h, l, r, t, b)
    reset_transform(gc)
    _set_coordinates(gc, x, y, w, h, l, r, t, b)
end

function _set_coordinates(gc::GraphicsContext, x, y, w, h, l, r, t, b)
    if (r-l) != w || (b-t) != h || l != x || t != y
        # note: Cairo assigns integer pixel-space coordinates to the grid
        # points between sample locations, not to the centers of pixels.
        xs = w/(r-l)
        ys = h/(b-t)
        scale(gc, xs, ys)

        xcent = (l+r)/2
        ycent = (t+b)/2
        translate(gc, -xcent + (w/2 + x)/xs, -ycent + (h/2 + y)/ys)
    end
    gc
end

"""
    set_coordinates(gc::GraphicsContext, device::BoundingBox, user::BoundingBox)
    set_coordinates(gc::GraphicsContext, user::BoundingBox)

Set the device->user coordinate transformation of `c` so that
`device`, expressed in "device coordinates" (pixels), is equivalent to
`user` as expressed in "user coordinates". If `device` is omitted, it
defaults to the full span of `gc`,
`BoundingBox(0, width(gc), 0, height(gc))`.
"""
set_coordinates(gc::GraphicsContext, device::BoundingBox, user::BoundingBox) =
    set_coordinates(gc,
                    device.xmin, device.ymin, width(device), height(device),
                    user.xmin, user.xmax, user.ymin, user.ymax)
set_coordinates(gc::GraphicsContext, user::BoundingBox) =
    set_coordinates(gc, BoundingBox(0, width(gc), 0, height(gc)), user)

"""
    save(gc::GraphicsContext)

Save the copy of current context `gc`. The context is saved onto an internal
stack, for example.

See also: [`restore`](@ref)

!!! warning
    The function name `save` conflicts with the `save` in the FileIO package,
    which is likely to be used with the Graphics package.
"""
@mustimplement save(gc::GraphicsContext)

"""
    restore(gc::GraphicsContext)

Restore the context saved by a preceding [`save`](@ref). The state at the time
of the call is overwritten with the restored context, and the restored context
is removed from an internal stack, for example.
"""
@mustimplement restore(gc::GraphicsContext)

"""
    reset_transform(gc::GraphicsContext)

Reset the current transformation.
"""
@mustimplement reset_transform(gc::GraphicsContext)

"""
    rotate(gc::GraphicsContext, angle)

Rotate the user-space axes by `angle` (in radians). The rotation takes places
after any existing transformation.

See also: [`rotate(p::Vec2, angle::Real, o::Vec2)`](@ref).
"""
@mustimplement rotate(gc::GraphicsContext, ::Real)

"""
    scale(gc::GraphicsContext, sx, sy)

Scale the user-space x-axis and y-axis by `sx` and `sy` respectively. The
scaling takes places after any existing transformation.
"""
@mustimplement scale(gc::GraphicsContext, ::Real, ::Real)

"""
    translate(gc::GraphicsContext, Δx, Δy)

Translate the user-space origin by `(Δx, Δy)`. The translation takes places
after any existing transformation.
"""
@mustimplement translate(gc::GraphicsContext, ::Real, ::Real)

"""
    user_to_device!(gc::GraphicsContext, c::Vector{Float64})

Transform a coordinate `c` from the user space to the device space.

See also: [`user_to_device`](@ref), [`device_to_user!`](@ref)
"""
user_to_device!(gc::GraphicsContext, c::Vector{Float64}) = c # TODO: generalization

"""
    device_to_user!(gc::GraphicsContext, c::Vector{Float64})

Transform a coordinate `c` from the device space to the user space.

See also: [`device_to_user`](@ref), [`user_to_device!`](@ref)
"""
device_to_user!(gc::GraphicsContext, c::Vector{Float64}) = c # TODO: generalization

"""
    user_to_device_distance!(gc::GraphicsContext, d::Vector{Float64})

Transform a distance vector `d` from the user space to the device space. This
function is similar to the [`device_to_user!`](@ref) except that the translation
components will be cancelled.
"""
user_to_device_distance!(gc::GraphicsContext, c::Vector{Float64}) = c # TODO: generalization

"""
    device_to_user_distance!(gc::GraphicsContext, d::Vector{Float64})

Transform a distance vector `d` from the device space to the user space. This
function is similar to the [`user_to_device!`](@ref) except that the translation
components will be cancelled.
"""
device_to_user_distance!(gc::GraphicsContext, c::Vector{Float64}) = c # TODO: generalization

const d2ubuf = zeros(2)

"""
    user_to_device(gc::GraphicsContext, x, y) -> (xd, yd)

Transform a user space coordinate `(x, y)` to the device space coordinate
`(xd, yd)`.

See also: [`user_to_device!`](@ref), [`device_to_user`](@ref)
"""
function user_to_device(gc::GraphicsContext, x::Real, y::Real) # FIXME: avoid use of global variable
    d2ubuf[1] = x
    d2ubuf[2] = y
    user_to_device!(gc, d2ubuf)
    d2ubuf[1], d2ubuf[2]
end

"""
    device_to_user(gc::GraphicsContext, x, y) -> (xu, yu)

Transform a device space coordinate `(x, y)` to the user space coordinate
`(xu, yu)`.

See also: [`device_to_user!`](@ref), [`user_to_device`](@ref)
"""
function device_to_user(gc::GraphicsContext, x::Real, y::Real) # FIXME: avoid use of global variable
    d2ubuf[1] = x
    d2ubuf[2] = y
    device_to_user!(gc, d2ubuf)
    d2ubuf[1], d2ubuf[2]
end

# drawing and properties

"""
    set_line_width(gc::GraphicsContext, w)

Set the current line width (in device-depended units). The actual width and
aspect-ratio on the screen may be affected by the transformation.
"""
@mustimplement set_line_width(gc::GraphicsContext, ::Real)

"""
    set_dash(gc::GraphicsContext, dashes::Vector{Float64}, offset)

Set the dash pattern. The `dashes` is a `Vector` of positive lengths. The
odd-numbered elements represent the length of the "on" state, and
the even-numbered elements represent the length of the "off" (blank) state. The
`offset` specifies an offset at which the stroke begins.
If `dashes` is empty, dashing is disabled.
"""
@mustimplement set_dash(gc::GraphicsContext, ::Vector{Float64}, ::Real)

"""
    set_source_rgb(gc::GraphicsContext, r, g, b)

Set the source pattern to an opaque color `RGB(r, g, b)`. The color components
are in the range `[0, 1]`, not `[0, 255]`.

See also: [`set_source_rgba`](@ref).
"""
@mustimplement set_source_rgb(gc::GraphicsContext, ::Real, ::Real, ::Real)

"""
    set_source_rgba(gc::GraphicsContext, r, g, b, a)

Set the source pattern to a transparent color `RGBA(r, g, b, a)`. The color and
alpha components are in the range `[0, 1]`, not `[0, 255]`.

See also: [`set_source_rgb`](@ref).
"""
@mustimplement set_source_rgba(gc::GraphicsContext, ::Real, ::Real, ::Real, ::Real)

"""
    set_source(gc::GraphicsContext, src)

Set the source pattern to `src`. If the `src` is a `Color`, it is used as the
source color.
"""
@mustimplement set_source(gc::GraphicsContext, src)

function set_source(gc::GraphicsContext, c::Color)
    rgb = convert(RGB, c)
    set_source_rgb(gc, rgb.r, rgb.g, rgb.b)
end

"""
    clip(gc::GraphicsContext)

Set a new clipping region based on the current path and fill within the context.

See also: [`reset_clip`](@ref), [`clip_preserve`](@ref).
"""
@mustimplement clip(gc::GraphicsContext)

"""
    clip_preserve(gc::GraphicsContext)

Set a new clipping region based on the current path and fill within the context.
Unlike the [`clip`](@ref) function, this function preserves the path within the
context.
"""
@mustimplement clip_preserve(gc::GraphicsContext)

"""
    reset_clip(gc::GraphicsContext)

Remove the clipping region set by the [`clip`](@ref) or [`clip_preserve`](@ref)
function.
"""
@mustimplement reset_clip(gc::GraphicsContext)


"""
    move_to(gc::GraphicsContext, x, y)

Begin a new sub-path. The current point will be moved to `(x, y)`.

See also: [`rel_move_to`](@ref).
"""
@mustimplement move_to(gc::GraphicsContext, ::Real, ::Real)

"""
    line_to(gc::GraphicsContext, x, y)

Add a line to the current path from the current point to the position `(x, y)`.
The current point will be moved to `(x, y)`.

See also: [`rel_line_to`](@ref).
"""
@mustimplement line_to(gc::GraphicsContext, ::Real, ::Real)

"""
    rel_move_to(gc::GraphicsContext, Δx, Δy)

Begin a new sub-path. The current point will be moved to `(x + Δx, y + Δy)`,
where `(x, y)` is the previous current point.

See also: [`move_to`](@ref).
"""
@mustimplement rel_move_to(gc::GraphicsContext, ::Real, ::Real)

"""
    rel_line_to(gc::GraphicsContext, Δx, Δy)

Add a line to the current path from the current point of `(x, y)` to the
position `(x + Δx, y + Δy)`. The current point will be moved to
`(x + Δx, y + Δy)`.

See also: [`line_to`](@ref).
"""
@mustimplement rel_line_to(gc::GraphicsContext, ::Real, ::Real)

"""
    arc(gc::GraphicsContext, xc, yc, radius, angle1, angle2)

Add a circular arc with the specified `radius` to the current path. The arc is
centered at `(xc, yc)`, begins at `angle1` and ends at `angle2`. The `angle1`
and `angle2` are in radians. The arc will be drawn in the direction of
increasing angles.
"""
@mustimplement arc(gc::GraphicsContext, ::Real, ::Real, ::Real, ::Real, ::Real)

"""
    close_path(gc::GraphicsContext)

Add a line to the current path from the current point to the beginning of the
current sub-path and closes the sub-path. The current point will be changed to
the joined point.

There is a difference between closing a subpath and drawing a line to the
equivalent coordinate. This difference might be visualized as a difference in
drawing stroke endpoints.
"""
@mustimplement close_path(gc::GraphicsContext)

"""
    new_path(gc::GraphicsContext)

Clear the current path.

!!! note
    Depending on the backend, the new path may actually begin when a path
    element is added.
"""
@mustimplement new_path(gc::GraphicsContext)

"""
    new_sub_path(gc::GraphicsContext)

Begin a new sub-path. The current point will be cleared.

See also: [`move_to`](@ref).
"""
@mustimplement new_sub_path(gc::GraphicsContext)


"""
    fill(gc::GraphicsContext)

Fill the current path according to the current fill rule. The current path will
be cleared from the context.

See also: [`fill_preserve`](@ref).
"""
@mustimplement fill(gc::GraphicsContext)

"""
    fill_preserve(gc::GraphicsContext)

Fill the current path according to the current fill rule. Unlike the
[`fill`](@ref) function, this function preserves the current path within the
context.
"""
@mustimplement fill_preserve(gc::GraphicsContext)

"""
    stroke(gc::GraphicsContext)

Stroke the current path according to the current stroke style. The current path
will be cleared from the context.

!!! note
    The `stroke` function ignores the current transformation. If you want to
    apply the current transformation to the stroke, use
    [`stroke_transformed`](@ref).

See also: [`stroke_preserve`](@ref).
"""
@mustimplement stroke(gc::GraphicsContext)

"""
    stroke_preserve(gc::GraphicsContext)

Stroke the current path according to the current stroke style. Unlike the
[`stroke`](@ref) function, this function preserves the current path within
the context.

!!! note
    The `stroke_preserve` function ignores the current transformation. If you
    want to apply the current transformation to the stroke, use
    [`stroke_transformed_preserve`](@ref).

"""
@mustimplement stroke_preserve(gc::GraphicsContext)

"""
    paint(gc::GraphicsContext)

Paint the current source everywhere within the current clipping region.
"""
@mustimplement paint(gc::GraphicsContext)

"""
    stroke_transformed(gc::GraphicsContext)

Stroke the current path according to the current stroke style and the current
transformation.

See also: [`stroke`](@ref).
"""
stroke_transformed(gc::GraphicsContext) = stroke(gc)

"""
    stroke_transformed_preserve(gc::GraphicsContext)

Stroke the current path according to the current stroke style and the current
transformation. Unlike the [`stroke_transformed`](@ref) function, this function
preserves the current path within the context.

See also: [`stroke_preserve`](@ref).
"""
stroke_transformed_preserve(gc::GraphicsContext) = stroke_preserve(gc)

# generic path functions

"""
    rectangle(gc::GraphicsContext, x, y, width, height)
    rectangle(gc::GraphicsContext, user::BoundingBox)

Add a sub-path rectangle to the current path. The `x` and `y` specify the
(typically the upper left) corner coordinate of the rectangle, and the `width`
and `height` specify the size.

You can also specify the position and size by a `BoundingBox` in user-space
coordinate.
"""
function rectangle(gc::GraphicsContext, x::Real, y::Real, width::Real, height::Real)
    move_to(gc, x, y)
    rel_line_to(gc, width, 0)
    rel_line_to(gc, 0, height)
    rel_line_to(gc, -width, 0)
    close_path(gc)
end
rectangle(gc::GraphicsContext, user::BoundingBox) =
    rectangle(gc, user.xmin, user.ymin, width(user), height(user))

"""
    circle(gc::GraphicsContext, x, y, r)

Add a sub-path circle to the current path. The `x` and `y` specify the center
coordinate of the circle, and the `r` specifies the radius.
"""
circle(gc::GraphicsContext, x::Real, y::Real, r::Real) =
    arc(gc, x, y, r, 0., 2pi)

"""
    polygon(gc::GraphicsContext, verts::Matrix, idx::Vector)

Add a closed sub-path polygon with the given vertices. The `verts` is a
collection of vertex coordinates in the following matrix form:
```julia
[x1 x2 x3 ... xn;
 y1 y2 y3 ... yn]
```
The `idx` is a vector of vertex indices, i.e. the matrix column numbers.

!!! tip
    You can reuse the vertex coordinates by specifying the same index in `idx`.
    This is useful when drawing meshes.
"""
function polygon(gc::GraphicsContext, verts::Matrix, idx::Vector)
    move_to(gc, verts[1,idx[1]], verts[2,idx[1]])
    for i=2:length(idx)
        n = idx[i]
        line_to(gc, verts[1,n], verts[2,n])
    end
    close_path(gc)
end

"""
    polygon(gc::GraphicsContext, points::AbstractVector)

Add a closed sub-path polygon with the given vertices to .
"""
function polygon(gc::GraphicsContext, points::AbstractVector)
    move_to(gc, points[1].x, points[1].y)
    for i in 2:length(points)
        line_to(gc, points[i].x, points[i].y)
    end
    close_path(gc)
end

end # module

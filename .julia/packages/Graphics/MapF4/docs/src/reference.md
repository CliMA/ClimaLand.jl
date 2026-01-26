# Reference

## Geometric primitives

```@docs
Vec2
Point
BoundingBox
```

## Geometry API

```@docs
aspect_ratio
center
deform
diagonal
isinside
shift
height
width
xmin
xmax
ymin
ymax
xrange
yrange
```

## 2d drawing contexts

```@docs
GraphicsDevice
GraphicsContext
creategc
getgc
```

## Coordinate systems

```@docs
set_coordinates
reset_transform
rotate
scale
translate
user_to_device!
device_to_user!
user_to_device_distance!
device_to_user_distance!
user_to_device
device_to_user
```

## Lines

```@docs
set_line_width
set_dash
```

## Colors and painting (drawing attributes)

```@docs
set_source
set_source_rgb
set_source_rgba
save
restore
```

## Clipping

```@docs
clip
clip_preserve
reset_clip
inner_canvas
```

## Paths

```@docs
move_to
line_to
rel_line_to
rel_move_to
new_path
new_sub_path
close_path
arc
```

## High-level paths

```@docs
rectangle
circle
polygon
```

## Fill and stroke

```@docs
fill
fill_preserve
paint
stroke
stroke_preserve
stroke_transformed
stroke_transformed_preserve
```

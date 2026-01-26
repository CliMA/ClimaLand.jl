# Tutorial

## TL;DR: If you read nothing else, read this

```@meta
CurrentModule = Contour
```

The most common use case for this package is plotting iso lines. Here's a
complete example that lets you do that, while showing off all of the most
important features of the package:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
plot(args...;kwargs...) = nothing # hide
for cl in levels(contours(x,y,z))
    lvl = level(cl) # the z-value of this contour level
    for line in lines(cl)
        xs, ys = coordinates(line) # coordinates of this line segment
        plot(xs, ys, color=lvl) # pseuod-code; use whatever plotting package you prefer
    end
end
```

## Preface: some test data...

The Contour module expects input data to be on a Cartesian grid, and supports
both uniform and non-uniform grid spacings.  For the following examples, `x` and
`y` are 1D sorted arrays that contain the grid coordinates, and `z` is a matrix
arranged such that `z[xi,yi]` correspond to the location `(x[xi], y[yi])`.

Let's consider the function ``z(x,y) = x^2 + y^2``:

```julia
x = -3:0.01:3
y = -4:0.02:5

z = [Float64((xi^2 + yi^2)) for xi in x, yi in y]
nothing # hide
```

`x` and `y` don't have to be evenly spaced - they can just as well be (sorted)
arrays of coordinate values.

## Example: plotting isolines

Usually, you'll start by calling [`contours`](@ref):

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
c = contours(x,y,z)
nothing # hide
```

The package is designed so that you shouldn't have to worry about the types of
the outputs - instead, there are functions that let you extract the data you
need. So, instead of simply returning a `Vector{ContourLevel}`, we return a
special object which supports the [`levels`](@ref) function. `levels` in turn
returns an iterable, where each item represents a contour level:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
c = contours(x,y,z) # hide
for cl in levels(c)
    # do something
end
```

On each level (`cl` in the snippet above) there are two pieces of information
that can be of interest. You find the ``z``-value of the isoline with the
[`level`](@ref) function, while [`lines`](@ref) yields an iterable collection
of line segments (remember that there might be more than one isoline for a given
``z``-value):

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
c = contours(x,y,z) # hide
cl = first(levels(c)) # hide
level(cl) # the z-value of the current isoline collection
lines(cl) # an iterable collection of isolines
nothing # hide
```

This contour level only had one line. An isoline is represented as a sequence of
vertices, which either starts and ends at the boundaries of the data set, or
closes on itself, in which case the first and last points are equal.

The ``x``- and ``y``-coordinates of an isoline can be extracted using the
[`coordinates`](@ref) or [`vertices`](@ref) functions:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
c = contours(x,y,z) # hide
cl = first(levels(c)) # hide
l = first(lines(cl))
xs, ys = coordinates(l)
x_y_pairs = vertices(l)
nothing # hide
```

Now we understand all the parts of the plotting example at the top:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
plot(args...; kwargs...) = nothing # hide
for cl in levels(contours(x,y,z))
    lvl = level(cl) # the z-value of this contour level
    for line in lines(cl)
        xs, ys = coordinates(line) # coordinates of this line segment
        plot(xs, ys, color=lvl) # pseuod-code; use whatever plotting package you prefer
    end
end
```

## Affecting the choice of contour levels

There are several ways to affect the choice of contour levels.

First, you can specify them manually:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
contours(x, y, z, [2,3])
nothing # hide
```

You can also just specify the number of levels you want, and let the package
choose them:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
contours(x, y, z, 2)
nothing # hide
```

The package uses [`Contour.contourlevels`](@ref) to choose the levels, so it's
entirely possible to investigate what levels would be traced without doing any
plotting:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
Contour.contourlevels(z, 4)
nothing # hide
```

If you only want a single contour level, use the [`contour`](@ref) function
directly - its fourth parameter is the ``z``-value at which to trace the isolines:

```@example
using Contour # hide
x = -3:0.01:3 # hide
y = -4:0.02:5 # hide
z = [Float64((xi^2 + yi^2)) for xi in x, yi in y] # hide
contour(x, y, z, 2.3)
nothing # hide
```

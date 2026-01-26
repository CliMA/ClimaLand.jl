# Lazy transformation of values

In image display and input/output, it is sometimes necessary to
transform the value (or the type) of individual pixels.  For example,
if you want to view an image with an unconventional range (e.g., -1000
to 1000, for which the normal range 0=black to 1=white will not be
very useful), then those values might need to be transformed before
display. Likewise, if try to save an image to disk that contains some
out-of-range or NaN values, you are likely to experience an error
unless the values are put in a range that makes sense for the specific
file format.

There are several approaches to handling this problem. One is to
compute a new image with scaled values, and for many users this may be
the simplest option.  However, particularly with large images (or
movies) this can present a performance problem.  In such cases, it's
better to separate the concept of the "map" (transformation) function
from the image (array) itself. (Here it's worth mentioning the
[MappedArrays](https://github.com/JuliaArrays/MappedArrays.jl)
package, which allows you to express lazy transformations on values
for an entire array.)

ImageCore contains several such transformation functions that are
frequently useful when working with images. Some of these functions
operate directly on values:

- [`clamp01`](@ref)
- [`clamp01nan`](@ref)

These two functions force the returned value to lie between 0 and 1,
or each color channel to lie between 0 and 1 for color
images. (`clamp01nan` forces `NaN` to 0, whereas `clamp01` does not
handle `NaN`.)

A simple application of these functions is in saving images, where you
may have some out-of-range values but don't care if they get
truncated:

```julia
img01 = clamp01nan.(img)
```

`img01` is safe to save to an image file, whereas trying to save `img`
might possibly result in an error (depending on the contents of
`img`).

Other functions require parameters:

- [`scaleminmax`](@ref)
- [`scalesigned`](@ref)
- [`colorsigned`](@ref)

These return a function rather than a value; that function can then
be applied to pixels of the image.  For example:

```julia
julia> f = scaleminmax(-10, 10)
(::#9) (generic function with 1 method)

julia> f(10)
1.0

julia> f(-10)
0.0

julia> f(5)
0.75
```

It's worth noting that you can combine these: for example, you can
combine `scalesigned` and `colorsigned` to map real values to linear
colormaps. For example, suppose we want to visualize some data,
mapping negative values to green hues and positive values to magenta
hues. Let's say the negative values are a bit more compressed, so
we're going to map -5 to pure green and +20 to pure magenta. We can
achieve this easily with the following:

```julia
julia> sc = scalesigned(-5, 0, 20)  # maps [-5, 0, 20] -> [-1, 0, 1]
(::#15) (generic function with 1 method)

julia> col = colorsigned()          # maps -1 -> green, +1->magenta
(::#17) (generic function with 1 method)

julia> f = x->col(sc(x))            # combine the two
(::#1) (generic function with 1 method)

julia> f(-5)
RGB{N0f8}(0.0,1.0,0.0)

julia> f(20)
RGB{N0f8}(1.0,0.0,1.0)

julia> f(0)
RGB{N0f8}(1.0,1.0,1.0)

julia> f(10)
RGB{N0f8}(1.0,0.502,1.0)
```

Finally, [`takemap`](@ref) exists to automatically set the parameters of
certain functions from the image itself.  For example,

    takemap(scaleminmax, A)

will return a function that scales the minimum value of `A` to 0 and
the maximum value of `A` to 1.

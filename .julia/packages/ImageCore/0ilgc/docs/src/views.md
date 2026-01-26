# Views

## View types defined in ImageCore

It is quite possible that the default representation of images will
satisfy most or all of your needs. However, to enhance flexibility in
working with image data, it is possible to leverage several different
kinds of "views." Generically, a view is an *interpretation* of array
data, one that may change the apparent meaning of the array but which
shares the same underlying storage: change an element of the view, and
you also change the original array. Views can facilitate processing images
of immense size without making copies, and writing algorithms in the
most convenient format often without having to worry about the
potential cost of converting from one format to another.

To illustrate views, it's helpful to begin with a very simple image:

```jldoctest
julia> using Colors

julia> img = [RGB(1,0,0) RGB(0,1,0);
              RGB(0,0,1) RGB(0,0,0)]
2×2 Matrix{RGB{FixedPointNumbers.N0f8}}:
 RGB(1.0, 0.0, 0.0)  RGB(0.0, 1.0, 0.0)
 RGB(0.0, 0.0, 1.0)  RGB(0.0, 0.0, 0.0)
```

which displays as

![rgbk](assets/rgbk.png)

```@meta
DocTestSetup = quote
    using Colors, ImageCore
    img = [RGB(1,0,0) RGB(0,1,0);
           RGB(0,0,1) RGB(0,0,0)]
    v = channelview(img)
    r = rawview(v)
end
```

Most commonly, it's convenient that all dimensions of this array
correspond to pixel indices: you don't need to worry about some
dimensions of the array corresponding to "color channels" and other
the spatial location, and you're guaranteed to get the entire pixel
contents when you access that location.

That said, occasionally there are reasons to want to treat `RGB` as a
3-component vector.  That's motivation for introducing our first view:

```jldoctest
julia> v = channelview(img)
3×2×2 reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}}) with eltype N0f8:
[:, :, 1] =
 1.0  0.0
 0.0  0.0
 0.0  1.0

[:, :, 2] =
 0.0  0.0
 1.0  0.0
 0.0  0.0
```

`channelview` does exactly what the name suggests: provide a view of
the array using separate channels for the color components.

To access the underlying representation of the `N0f8` numbers, there's
another view called `rawview`:

```jldoctest
julia> r = rawview(v)
3×2×2 rawview(reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}})) with eltype UInt8:
[:, :, 1] =
 0xff  0x00
 0x00  0x00
 0x00  0xff

[:, :, 2] =
 0x00  0x00
 0xff  0x00
 0x00  0x00
```

Let's make a change in one of the entries:

```jldoctest
julia> r[3,1,1] = 128
128
```

If we display `img`, now we get this:

![mgbk](assets/mgbk.png)

You can see that the first pixel has taken on a magenta hue, which is
a mixture of red and blue.  Why does this happen? Let's look at the
array values themselves:

```@meta
DocTestSetup = quote
    using Colors, ImageCore
    img = [RGB(1,0,0) RGB(0,1,0);
           RGB(0,0,1) RGB(0,0,0)]
    v = channelview(img)
    r = rawview(v)
    r[3,1,1] = 128
end
```

```jldoctest
julia> r
3×2×2 rawview(reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}})) with eltype UInt8:
[:, :, 1] =
 0xff  0x00
 0x00  0x00
 0x80  0xff

[:, :, 2] =
 0x00  0x00
 0xff  0x00
 0x00  0x00

julia> v
3×2×2 reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}}) with eltype N0f8:
[:, :, 1] =
 1.0    0.0
 0.0    0.0
 0.502  1.0

[:, :, 2] =
 0.0  0.0
 1.0  0.0
 0.0  0.0

julia> img
2×2 Matrix{RGB{N0f8}}:
 RGB(1.0, 0.0, 0.502)  RGB(0.0, 1.0, 0.0)
 RGB(0.0, 0.0, 1.0)    RGB(0.0, 0.0, 0.0)
```

The hexadecimal representation of 128 is 0x80; this is approximately
halfway to 255, and as a consequence the `N0f8` representation is
very near 0.5.  You can see the same change is reflected in `r`, `v`,
and `img`: there is only one underlying array, `img`, and the two
views simply reference it.

Maybe you're used to having the color channel be the last dimension,
rather than the first. We can achieve that using `PermutedDimsArray`:

```jldoctest
julia> p = PermutedDimsArray(v, (2,3,1))
2×2×3 PermutedDimsArray(reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}}), (2, 3, 1)) with eltype N0f8:
[:, :, 1] =
 1.0  0.0
 0.0  0.0

[:, :, 2] =
 0.0  1.0
 0.0  0.0

[:, :, 3] =
 0.502  0.0
 1.0    0.0

julia> p[1,2,:] .= 0.25
3-element view(PermutedDimsArray(reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}}), (2, 3, 1)), 1, 2, :) with eltype N0f8:
 0.251N0f8
 0.251N0f8
 0.251N0f8

julia> p
2×2×3 PermutedDimsArray(reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}}), (2, 3, 1)) with eltype N0f8:
[:, :, 1] =
 1.0  0.251
 0.0  0.0

[:, :, 2] =
 0.0  0.251
 0.0  0.0

[:, :, 3] =
 0.502  0.251
 1.0    0.0

julia> v
3×2×2 reinterpret(reshape, N0f8, ::Matrix{RGB{N0f8}}) with eltype N0f8:
[:, :, 1] =
 1.0    0.0
 0.0    0.0
 0.502  1.0

[:, :, 2] =
 0.251  0.0
 0.251  0.0
 0.251  0.0

julia> img
2×2 Matrix{RGB{N0f8}}:
 RGB(1.0, 0.0, 0.502)  RGB(0.251, 0.251, 0.251)
 RGB(0.0, 0.0, 1.0)    RGB(0.0, 0.0, 0.0)
```

Once again, `p` is a view, and as a consequence changing it leads to
changes in all the coupled arrays and views.

Finally, you can combine multiple arrays into a "virtual" multichannel
array. We'll use the
[lighthouse](http://juliaimages.github.io/TestImages.jl/images/lighthouse.png)
image:

```@example
using ImageCore, TestImages, Colors
img = testimage("lighthouse")
# Split out into separate channels
cv = channelview(img)
# Recombine the channels, filling in 0 for the middle (green) channel
rb = colorview(RGB, cv[1,:,:], zeroarray, cv[3,:,:])
```

`zeroarray` is a constant which serves as a placeholder to create a
(virtual) all-zeros array of size that matches the other arguments.

`rb` looks like this:

![redblue](assets/redblue.png)

In this case, we could have done the same thing somewhat more simply
with `cv[2,:,:] .= 0` and then visualize `img`. However, more generally
you can apply this to independent arrays which may not allow you to
set values to 0. In IJulia,

![linspace1](assets/linspace1.png)

The error comes from the fact that `img1d` does not store values
separately from the `LinSpace` objects used to create it, and
`LinSpace` (which uses a compact representation of a range, storing
just the endpoints and the number of values) does not allow you to set
specific values. However, if you need to set individual values, you
can make a `copy`:

![linspace2](assets/linspace2.png)

The fact that no storage is allocated by `colorview`
is very convenient in certain situations, particularly when processing
large images.

`colorview`'s ability to combine multiple grayscale images is based on
another view, `StackedView`, which you can also use directly.

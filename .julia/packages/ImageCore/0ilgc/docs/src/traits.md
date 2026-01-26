# Traits

ImageCore supports several "traits" that are sometimes useful in
viewing or analyzing images. Many of these traits become much more
powerful if you are using add-on packages like ImagesAxes, which
allows you to give "physical meaning" to the different axes of your
image.  Readers are encouraged to view the documentation for ImageAxes
to gain a better appreciation of how to exploit these traits.  When
using plain arrays to represent images, most of the traits default to
"trivial" outcomes.

Let's illustrate with a couple of examples:

```julia
julia> using Colors, ImageCore

julia> img = rand(RGB{N0f8}, 680, 480);

julia> pixelspacing(img)
(1,1)
```

`pixelspacing` returns the spacing between adjacent pixels along each
axis. Using ImagesAxes, you can even use physical units to encode this
information, which might be important for microscopy or biomedical imaging.

```@meta
DocTestSetup = quote
    using Colors, ImageCore
    img = rand(RGB{N0f8}, 680, 480);
end
```

Another simple trait is `coords_spatial`:

```julia
julia> coords_spatial(img)
(1,2)
```

This trait indicates that both dimensions 1 and 2 are "spatial
dimensions," meaning they correspond to physical space. This trait
again becomes more interesting with ImagesAxes, where you can denote
that some axes correspond to time (e.g., for a movie).

A full list of traits is presented in the reference section.

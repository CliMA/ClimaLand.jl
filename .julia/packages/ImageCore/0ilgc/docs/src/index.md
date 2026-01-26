# ImageCore.jl

ImageCore is the lowest-level component of the system of packages
designed to support image processing and computer vision. Its main
role is to simplify "conversions" between different image
representations through different "view" types, and to provide some
useful low-level functions (including "traits") that simplify image
display, input/output, and the writing of algorithms.

Some of the key features and functionalities provided by ImageCore.jl include:

- Storage type transformations: ImageCore.jl provides various methods like `float32`, `n0f8`, etc which can be used to transform raw storage type of the image without changing the color space. 
- Image view conversion: It provides various methods like `rawview`, `channelview`, `colorview` , `stackedview` which help in conversions between various types of representations. Default representation of images is likely all that is needed most of the time but these functions provided more flexibility while working with image data.
- Pixel level transformations: Several methods like `clamp01` , `clamp01nan`, `scaleminmax`, etc are available that help with truncating image pixel values within certain valid limits and helps to avoid IO errors.
- Image Traits: These can prove to be useful to users when trying to add more meaning to different axes of image and to add dimension specific information.


ImageCore often acts as the "common core" dependency for working with any kind of image in the JuliaImages ecosystem, but itself provides almost none of the standard image-processing algorithms. These algorithms are found in other packages in the JuliaImages ecosystem.

If you want to use ImageCore.jl, you can add it to your project environment using the package manager. Open a Julia REPL and use the following command:

```jl
] add ImageCore
```


If you're just getting started with images in Julia, it's recommended
that you see the
[introductory documentation](http://juliaimages.github.io/latest/). In
particular, this document assumes that you understand how Julia
represents color through the used of fixed-point numbers.

```@contents
Pages = ["views.md", "map.md", "traits.md", "reference.md"]
```

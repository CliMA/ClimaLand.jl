# # Writing TIFFs

# This page is a tutorial for saving TIFFs using TiffImages.jl and covers some common
# use cases

#md # ```@contents
#md # Pages = ["writing.md"]
#md # Depth = 5
#md # ```

# You might want to write TIFFs to disk too. Now this can be done quite simply
# with TiffImages.jl. Say you have some AbstractArray type that you want to save, here
# we'll call it `data`:

using Random
using Images # for nice inline images

Random.seed!(123)
data = rand(RGB{N0f8}, 10, 10)

#md # !!! note
#md #     TiffImages.jl only works with AbstractArrays with `eltype`s of `<:ColorOrTuple` because
#md #     the writer needs to know how to represent the image data on disk. Make sure to
#md #     convert your `AbstractArrays` using before passing them. See the
#md #     [common strategies](#Strategies-for-saving-common-types) section
#md #     below for tips.

# ## Simple cases

# In most simple cases, all you need to do is use the `save` function

using TiffImages
TiffImages.save("test.tif", data)

# That's it! TiffImages will convert your data into its own internal file type
# and then rapidly write it to disk. See the "Incremental writing" section of
# [Lazy TIFFs](@ref) for building a TIFF piece by piece.

# ## Complex cases
# If you need more fine-grained control over what tags are included when the
# image is written, this section is for you!

# ### Converting to `TiffImages.jl`'s TIFF type
# Next lets convert `data` to a TIFF type

using TiffImages
img = TiffImages.DenseTaggedImage(data)

# Wait nothing happened! Hang with me, lets take a closer look at our new object
# using the `dump` command. We can see that there's now new information
# associated with our data! TiffImages.jl usually represents TIFF images as simply the
# data and associated tags that describe the data

dump(img; maxdepth=1)

# The tags are organized as a vector of what are called [Image File
# Directories](https://www.awaresystems.be/imaging/tiff/faq.html) (IFDs). For a
# simple 2D image like what we have, the IFDs will be stored a vector of
# length=1. For 3D images, the length of the IFDs vector will equal the length
# of the image in the third dimension.

# Lets take a look at what tags there are:

ifd = ifds(img) # returns a single IFD since our data is 2D

# ### Manipulating TIFF Tags
#
# These are some of the most basic tags that are required by the TIFF spec. We
# can even update it to add our own custom tags

ifd[TiffImages.IMAGEDESCRIPTION] = "This is very important data"
ifd

# We can even add tags that aren't in the standard set in
# [`TiffImages.TiffTag`](@ref) as long as they are a `UInt16`

ifd[UInt16(34735)] = UInt16[1, 2, 3]
ifd

# We can also delete tags if we decide we don't want them:

delete!(ifd, TiffImages.IMAGEDESCRIPTION)
ifd

#md # !!! warning
#md #     Careful with `delete!`, if any of core tags are deleted, TiffImages.jl and
#md #     other readers might fail to read the file


# ### Saving to disk
#
# Once you're happy with your TIFF object, you can write it to disk as follows:

TiffImages.save("test.tif", img)

# And to just double check, we can load it right back in

TiffImages.load("test.tif")

# ## Strategies for saving common types

# The general strategy for saving arrays will differ a bit depending on the
# type. The key step is the convert or reinterpret the arrays so that the
# elements are subtypes of `ColorOrTuple`

# #### Unsigned Integers

# Say you want to save a 3D array of small integers as grayscale values.

data2 = rand(UInt8.(1:255), 5, 10)
eltype(data2)

# You can't directly save the `data2` since TiffImages.jl needs some color information
# to properly save the file. You can use
# [`reinterpret`](https://docs.julialang.org/en/v1/base/arrays/#Base.reinterpret)
# to accomplish this:

grays = reinterpret(Gray{N0f8}, data2)
img2 = TiffImages.DenseTaggedImage(grays)

# Here the data are first reinterpreted as `N0f8`s, which is a
# [`FixedPointNumber`](https://github.com/JuliaMath/FixedPointNumbers.jl) then
# wrapped with a Gray type that marks this as a grayscale image. TiffImages.jl uses
# this information to update the TIFF tags

# #### Floating point numbers

# With RGB we can reinterpret the first dimension of a 3D array as the 3
# different color components (red, green, and blue):

data = rand(Float64, 3, 5, 10);
colors = dropdims(reinterpret(RGB{eltype(data)}, data), dims=1) # drop first dimension
img3 = TiffImages.DenseTaggedImage(colors)

# Here we dropped the first dimension since it was collapsed into the RGB type
# when we ran the `reinterpret` command.

# #### Signed integers

# Say you want to save data that has negative integer values. In that case, you
# can't use `N0f8`, etc because those only worked for unsigned integers. You
# have to instead use `Q0f63`, etc, which is a different kind of fixed point
# number that uses one bit for the sign info (that's why it's `Q0f63`, not `Q0f64`!)

data = rand(-100:100, 5, 5)
#--------------------------
img4 = TiffImages.DenseTaggedImage(reinterpret(Gray{Q0f63}, data))
println(ifds(img4))

# As you can see the `SAMPLEFORMATS` and `BITSPERSAMPLE` tags correctly updated
# to show that this TIFF contains signed integers and 64-bit data, respectively.

#md # !!! warning
#md #     Currently, several of the display libraries struggle with showing
#md #     `Colorant`s backed by a signed type so you might run into errors, but
#md #     the data will still save properly

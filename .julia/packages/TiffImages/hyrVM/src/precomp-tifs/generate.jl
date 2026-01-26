# this file must be run manually whenever the set of images for precompilation should be changed
using TiffImages
using ColorTypes
using FixedPointNumbers

types = [N0f8, N0f16, Float32, Float64]
colors = [Gray, RGB]
sizes = [(2, 2), (2, 2, 2)]

images = joinpath(@__DIR__, "images")
rm(images, recursive = true, force = true)
mkdir(images)

# generate tifs to load during precompilation, we can't generate them during precompilation
# because we don't want to spend the time calling `save` there (see issue #180)
files = map(enumerate(Iterators.product(types, colors, sizes))) do (i, (t, c, s))
    file = joinpath(@__DIR__, "images", "$(t)_$(c)_$(s).tif")
    TiffImages.save(file, rand(c{t}, s))
end
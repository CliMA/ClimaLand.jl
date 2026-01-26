# Reference

## List of view types

With that as an introduction, let's list all the view types supported
by this package.  `channelview` and `colorview` are opposite
transformations, as are `rawview` and `normedview`. `channelview` and
`colorview` typically create objects of type `ChannelView` and
`ColorView`, respectively, unless they are "undoing" a previous view
of the opposite type.

```@docs
channelview
colorview
rawview
normedview
StackedView
```

## List of value-transformations (map functions)

```@docs
clamp01
clamp01!
clamp01nan
clamp01nan!
scaleminmax
scalesigned
colorsigned
takemap
```

## List of storage-type transformations

```@docs
float32
float64
n0f8
n6f10
n4f12
n2f14
n0f16
```

## List of traits

```@docs
pixelspacing
spacedirections
sdims
coords_spatial
size_spatial
indices_spatial
nimages
assert_timedim_last
```

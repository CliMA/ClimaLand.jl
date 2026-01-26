# Extending TiffImages.jl

If you want to extend `TiffImages.jl` to add support for more features or change how
TIFF data is loaded, you have come to right place.

### Types

```@docs
TiffImages.TiffFile
TiffImages.IFD
TiffImages.Tag
TiffImages.Iterable
TiffImages.RemoteData
```

### Functions

```@docs
TiffImages.load!
TiffImages.sizeof
```
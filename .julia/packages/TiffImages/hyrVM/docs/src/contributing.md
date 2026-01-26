# Contributing

Supporting all TIFFs is non-trivial and I would greatly appreciate any help from
the community in identifying edge cases.

## Add edge case TIFFs

There is incredible diversity in the TIFF ecosystem so much so that there is a
backronym "Thousand Incompatible File Formats" to describe it. I have tried to
establish a good baseline test set of TIFFs that should guarantee that `TiffImages.jl`
should "just work tm" for most people, but if you have a TIFF that you run into
that breaks `TiffImages.jl` please do the following:

1. create a pull request against the
   [example TIFF](https://github.com/tlnagy/exampletiffs) repo adding the file.
   The smaller the file, the better.
2. update the README table with license information, etc.
3. Open an issue against `TiffImages.jl` with the error message and the expected result
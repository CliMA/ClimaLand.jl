# ---
# title: Loading a remote TIFF
# description: This demo shows how to load a TIFF over the network
# author: Santiago Pelufo, Tamas Nagy
# cover: assets/coffee.png
# ---

# The following code can be used to load remote TIFFs without saving them to the
# disk. 

using TiffImages, Downloads, FileIO

function load_tiff_without_saving(url)
  buffer = IOBuffer()
  Downloads.download(url, buffer)
  bufstream = TiffImages.getstream(format"TIFF", buffer)
  TiffImages.load(read(bufstream, TiffFile))
end

# We'll load an example from
# [`tlnagy/exampletiffs`](https://github.com/tlnagy/exampletiffs) 

coffeepath = "https://github.com/tlnagy/exampletiffs/raw/11516d288c4b03a258aa3027705b0e9d2ce2b5de/coffee.tif"
coffeeimg = load_tiff_without_saving(coffeepath)

# No disks involved!

save("assets/coffee.png", collect(coffeeimg[:, 70:470])) #hide
nothing #hide
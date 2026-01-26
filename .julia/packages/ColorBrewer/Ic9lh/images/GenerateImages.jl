using Colors, ColorBrewer

function do_palette(name="Blues", size::Int=9)
  width, height = "1cm", "0.5cm"
  path = "$colorbrewer_git_path/images"
  svg_filename = "$path/$name.svg"
  png_filename = "$path/$name.png"
  f = open(svg_filename, "w")
  @printf f "%s\n" "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"$(size)cm\" height=\"0.5cm\">"
  x = 0
  for c in palette(name, size)
    r, g, b = int(round(c.r*255)), int(round(c.g*255)), int(round(c.b*255))
    @printf f "%s\n" "<rect x=\"$(x)cm\" width=\"$width\" height=\"$height\" fill=\"rgb($r, $g, $b)\"/>"
    x += 1
  end
  @printf f "</svg>"
  close(f)
  run(`svg2png $svg_filename $png_filename`)
  run(`rm $svg_filename`)
end

function main()
  for name in ["Blues","Oranges","Greens","Reds","Purples","Greys","OrRd","GnBu","PuBu","PuRd","BuPu","BuGn","YlGn","RdPu","YlOrBr","YlGnBu","YlOrRd","PuBuGn"]
    do_palette(name, 9)
  end
  for name in ["Spectral","RdYlGn","RdBu","PiYG","PRGn","RdYlBu","BrBG","RdGy","PuOr"]
    do_palette(name, 11)
  end
  do_palette("Set1", 9)
  do_palette("Set2", 8)
  do_palette("Set3", 12)
  do_palette("Accent", 8)
  do_palette("Dark2", 8)
  do_palette("Paired", 12)
  do_palette("Pastel1", 9)
  do_palette("Pastel2", 8)
end

nothing

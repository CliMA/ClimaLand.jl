packer = GuillotinePacker(1024, 1024)
rects = [Rect(0, 0, rand(10:100), rand(10:100)) for i in 1:300]
push!(packer, rects)
push!(packer, Rect(0, 0, 100, 100))

# using Makie
# linesegments(Rectf(0, 0, 1024, 1024))
# poly!(packer.used_rectangles, color = (:red, 0.1), strokewidth=1, strokecolor=:black)

root = RectanglePacker(Rect2(0,0,1024,1024))
push!(root, Rect2(0,0,20,20))
push!(root, Rect2(0,0,20,20))
push!(root, [Rect2(0,0, rand(5:50), rand(5:50)) for i=1:20])


# function get_rectangles(packer::Nothing, rectangles=Rect2i[])
#     return rectangles
# end
#
# function get_rectangles(packer, rectangles)
#     push!(rectangles, packer.area)
#     get_rectangles(packer.children.left, rectangles)
#     get_rectangles(packer.children.right, rectangles)
#     return rectangles
# end
# function get_rectangles(packer)
#     rectangles = Rect2i[]
#     get_rectangles(packer.children.left, rectangles)
#     get_rectangles(packer.children.right, rectangles)
#     return rectangles
# end
#
# rectangles = get_rectangles(root)
# using Makie
# linesegments(Rectf(0, 0, 1024, 1024))
# poly!(rectangles, color = (:red, 0.1), strokewidth=1, strokecolor=:black)
#

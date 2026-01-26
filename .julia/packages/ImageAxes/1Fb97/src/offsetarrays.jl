if isdefined(OffsetArrays, :centered)
    # Requires OffsetArrays v1.9
    # https://github.com/JuliaArrays/OffsetArrays.jl/pull/242
    # Originally from ImageFiltering https://github.com/JuliaImages/ImageFiltering.jl/pull/214
    OffsetArrays.centered(ax::Axis{name}) where name = Axis{name}(OffsetArrays.centered(ax.val))
    OffsetArrays.centered(a::AxisArray) = AxisArray(OffsetArrays.centered(a.data), OffsetArrays.centered.(AxisArrays.axes(a)))
end

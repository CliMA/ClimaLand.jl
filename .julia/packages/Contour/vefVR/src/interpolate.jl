# Given the row and column indices of the lower left
# vertex, add the location where the contour level
# crosses the specified edge.
function interpolate(x, y, z::AbstractMatrix, h::Number, ind, edge::UInt8, ::Type{VT}) where {VT}
    xi, yi = ind
    @inbounds if edge == W
        y_interp = y[yi] + (y[yi + 1] - y[yi]) * (h - z[xi, yi]) / (z[xi, yi + 1] - z[xi, yi])
        x_interp = x[xi]
    elseif edge == E
        y_interp = y[yi] + (y[yi + 1] - y[yi]) * (h - z[xi + 1, yi]) / (z[xi + 1, yi + 1] - z[xi + 1, yi])
        x_interp = x[xi + 1]
    elseif edge == N
        y_interp = y[yi + 1]
        x_interp = x[xi] + (x[xi + 1] - x[xi]) * (h - z[xi, yi + 1]) / (z[xi + 1, yi + 1] - z[xi, yi + 1])
    elseif edge == S
        y_interp = y[yi]
        x_interp = x[xi] + (x[xi + 1] - x[xi]) * (h - z[xi, yi]) / (z[xi + 1, yi] - z[xi, yi])
    end

    return VT <: Tuple ? (x_interp, y_interp) : VT(x_interp, y_interp)
end

function interpolate(x::AbstractRange, y::AbstractRange, z::AbstractMatrix, h::Number, ind, edge::UInt8, ::Type{VT}) where {VT}
    xi, yi = ind
    @inbounds if edge == W
        y_interp = y[yi] + step(y) * (h - z[xi, yi]) / (z[xi, yi + 1] - z[xi, yi])
        x_interp = x[xi]
    elseif edge == E
        y_interp = y[yi] + step(y) * (h - z[xi + 1, yi]) / (z[xi + 1, yi + 1] - z[xi + 1, yi])
        x_interp = x[xi + 1]
    elseif edge == N
        y_interp = y[yi + 1]
        x_interp = x[xi] + step(x) * (h - z[xi, yi + 1]) / (z[xi + 1, yi + 1] - z[xi, yi + 1])
    elseif edge == S
        y_interp = y[yi]
        x_interp = x[xi] + step(x) * (h - z[xi, yi]) / (z[xi + 1, yi] - z[xi, yi])
    end

    return VT <: Tuple ? (x_interp, y_interp) : VT(x_interp, y_interp)
end

function interpolate(x::AbstractMatrix, y::AbstractMatrix, z::AbstractMatrix, h::Number, ind, edge::UInt8, ::Type{VT}) where {VT}
    xi, yi = ind
    @inbounds if edge == W
        Δ = [y[xi,  yi+1] - y[xi,  yi  ], x[xi,  yi+1] - x[xi,  yi  ]].*(h - z[xi,  yi  ])/(z[xi,  yi+1] - z[xi,  yi  ])
        y_interp = y[xi,yi] + Δ[1]
        x_interp = x[xi,yi] + Δ[2]
    elseif edge == E
        Δ = [y[xi+1,yi+1] - y[xi+1,yi  ], x[xi+1,yi+1] - x[xi+1,yi  ]].*(h - z[xi+1,yi  ])/(z[xi+1,yi+1] - z[xi+1,yi  ])
        y_interp = y[xi+1,yi] + Δ[1]
        x_interp = x[xi+1,yi] + Δ[2]
    elseif edge == N
        Δ = [y[xi+1,yi+1] - y[xi,  yi+1], x[xi+1,yi+1] - x[xi,  yi+1]].*(h - z[xi,  yi+1])/(z[xi+1,yi+1] - z[xi,  yi+1])
        y_interp = y[xi,yi+1] + Δ[1]
        x_interp = x[xi,yi+1] + Δ[2]
    elseif edge == S
        Δ = [y[xi+1,yi  ] - y[xi,  yi  ], x[xi+1,yi  ] - x[xi,  yi  ]].*(h - z[xi,  yi  ])/(z[xi+1,yi  ] - z[xi,  yi  ])
        y_interp = y[xi,yi] + Δ[1]
        x_interp = x[xi,yi] + Δ[2]
    end

    return VT <: Tuple ? (x_interp, y_interp) : VT(x_interp, y_interp)
end

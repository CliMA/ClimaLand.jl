# Originally from Images.jl
restrict(img::AxisArray, ::Tuple{}) = img

function restrict(img::AxisArray, region::Dims)
    inregiont = ntuple(i->in(i, region), ndims(img))
    rdata = restrict(img.data, region)
    AxisArray(rdata, map(modax, AxisArrays.axes(img), axes(rdata), inregiont))
end

function restrict(img::AxisArray, ::Type{Ax}) where Ax
    dim = axisdim(img, Ax)
    A = restrict(img.data, dim)
    AxisArray(A, replace_axis(modax(img[Ax], axes(A)[dim]), AxisArrays.axes(img)))
end

replace_axis(newax, axs) = _replace_axis(newax, axnametype(newax), axs...)
@inline _replace_axis(newax, ::Type{Ax}, ax::Ax, axs...) where {Ax} = (newax, _replace_axis(newax, Ax, axs...)...)
@inline _replace_axis(newax, ::Type{Ax}, ax, axs...) where {Ax} = (ax, _replace_axis(newax, Ax, axs...)...)
_replace_axis(newax, ::Type{Ax}) where {Ax} = ()

axnametype(ax::Axis{name}) where {name} = Axis{name}

ofaxes(::Base.OneTo, r) = r
ofaxes(axs::Any, r) = OffsetArray(r, axs)

function modax(ax, Aax)
    v = ax.val
    if iseven(length(v))
        return ax(ofaxes(Aax, range(first(v)-step(v)/2, stop=last(v)+step(v)/2, length=length(v)÷2 + 1)))
    else
        return ax(ofaxes(Aax, range(first(v)/1, stop=last(v)/1, length=length(v)÷2 + 1)))
    end
end

modax(ax, Aax, inregion::Bool) = inregion ? modax(ax, Aax) : ax

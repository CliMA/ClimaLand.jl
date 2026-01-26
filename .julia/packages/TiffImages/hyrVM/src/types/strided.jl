struct StridedTaggedImage{O <: Unsigned, AA <: AbstractMatrix} <: AbstractTIFF{Any, 1}
    data::Vector{AA}
    ifds::Vector{IFD{O}}
end

Base.size(t::StridedTaggedImage) = size(t.data)
Base.getindex(img::StridedTaggedImage, i) = getindex(img.data, i)
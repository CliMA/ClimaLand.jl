# have to import this or @deprecate doesn't work
import Base: getindex, setindex!, delete!, haskey, get, copy!, getproperty, setproperty!

@deprecate(ImageMeta(data::AbstractArray{T,N}, props::AbstractDict{String,Any}) where {T,N},
           ImageMeta(data, to_dict(props)))

@deprecate(ImageMeta(data::AbstractArray{T,N}, props::Dict{String,Any}) where {T,N},
           ImageMeta(data, to_dict(props)))

@deprecate(ImageMeta(data::AbstractArray, props::Dict{<:AbstractString}),
           ImageMeta(data, to_dict(props)))

@deprecate(copy!(imgdest::ImageMeta, imgsrc::ImageMeta, prop1::AbstractString, props::AbstractString...),
           copy!(imgdest, imgsrc, Symbol(prop1), Symbol.(props)...))

@deprecate(delete!(img::ImageMeta, propname::AbstractString),
           delete!(img, Symbol(propname)))

@deprecate(haskey(img::ImageMeta, k::AbstractString),
           hasproperty(img, Symbol(k)))

@deprecate(hasproperty(img::ImageMeta, k::AbstractString),
           hasproperty(img, Symbol(k)))

@deprecate(get(img::ImageMeta, k::AbstractString, default),
           get(img, Symbol(k), default))

@deprecate(getindex(img::ImageMeta, propname::AbstractString),
           getproperty(img, Symbol(propname)))

@deprecate(getproperty(img::ImageMeta, propname::AbstractString),
           getproperty(img, Symbol(propname)))

@deprecate(setproperty!(img::ImageMeta, propname::AbstractString, val),
           setproperty!(img::ImageMeta, Symbol(propname), val))

@deprecate(setindex!(img::ImageMeta, X, propname::AbstractString),
           setproperty!(img, Symbol(propname), X))

@deprecate(data(img::ImageMeta), arraydata(img))

function to_dict(dold::AbstractDict{String})
    dnew = Dict{Symbol,Any}()
    for (k, v) in dold
        dnew[Symbol(k)] = v
    end
    return dnew
end

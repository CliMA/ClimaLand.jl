using ImageCore: AbstractGray, TransparentGray, TransparentRGB

# Specialize the low-level broadcasting machinery for ImageMeta.
# This make it work for all operators, rather than specializing
# each operator individually.
Base.BroadcastStyle(::Type{<:ImageMeta}) = Broadcast.ArrayStyle{ImageMeta}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ImageMeta}}, ::Type{ElType}) where ElType
    Mt = imagemeta(bc.args...)
    if length(Mt) > 0
        M = Mt[1]
        if length(Mt) > 1
            for i = 2:length(Mt)
                Mt[i] == M || return ambigop(:Broadcast)
            end
        end
        # Use the properties field of img to create the output
        ret = ImageMeta(similar(Array{ElType}, axes(bc)), properties(M))
    else
        ret = ImageMeta(similar(Array{ElType}, axes(bc)))
    end
    ret
end
# Select all the ImageMeta arrays
@inline imagemeta(As...) = _imagemeta((), As...)
_imagemeta(out) = out
@inline _imagemeta(out, A::ImageMeta, As...) = _imagemeta((out..., A), As...)
@inline _imagemeta(out, A, As...) = _imagemeta(out, As...)


(-)(img::ImageMeta) = shareproperties(img, -arraydata(img))

batch2 = (:+, :-)

for op in batch2
    @eval begin
        ($op)(img::ImageMeta, B::BitArray) = shareproperties(img, ($op)(arraydata(img), B))
        ($op)(B::BitArray, img::ImageMeta) = shareproperties(img, ($op)(B, arraydata(img)))
        ($op)(img::ImageMeta, B::ImageMeta) = ambigop(Symbol($op))
        ($op)(img::ImageMeta, B::AbstractArray) = shareproperties(img, ($op)(arraydata(img), B))
        ($op)(B::AbstractArray, img::ImageMeta) = shareproperties(img, ($op)(B, arraydata(img)))
    end

    for CV in (:AbstractGray, :TransparentGray, :AbstractRGB, :TransparentRGB)
        @eval begin
            ($op)(img::ImageMeta{CV}, n::$CV) where {CV<:$CV} =
                shareproperties(img, ($op)(arraydata(img), n))
            ($op)(n::$CV, img::ImageMeta{CV}) where {CV<:$CV} =
                shareproperties(img, ($op)(n, arraydata(img)))
        end
    end
end

(/)(img::ImageMeta, n::Number) = shareproperties(img, arraydata(img)/n)

ambigop(s::Symbol) = error("$s with two ImageMeta arrays: dictionary choice is ambiguous")

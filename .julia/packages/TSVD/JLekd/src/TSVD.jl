# Copyright 2015-2018 Andreas Noack
module TSVD

    export tsvd

    import Base: *, hcat, maximum, size
    import LinearAlgebra: axpy!, mul!

    using LinearAlgebra
    using LinearAlgebra: BlasComplex, BlasFloat, BlasInt, BlasReal
    using Adapt: adapt

    include("common.jl")
    include("eig.jl")
    include("svd.jl")

end

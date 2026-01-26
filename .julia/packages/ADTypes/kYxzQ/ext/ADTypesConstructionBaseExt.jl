module ADTypesConstructionBaseExt

using ADTypes: AutoEnzyme, AutoForwardDiff, AutoPolyesterForwardDiff
using ConstructionBase: ConstructionBase

struct InternalAutoEnzymeReconstructor{A} end

InternalAutoEnzymeReconstructor{A}(mode::M) where {M, A} = AutoEnzyme{M, A}(mode)

function ConstructionBase.constructorof(::Type{<:AutoEnzyme{M, A}}) where {M, A}
    return InternalAutoEnzymeReconstructor{A}
end

function ConstructionBase.constructorof(::Type{<:AutoForwardDiff{chunksize}}) where {chunksize}
    return AutoForwardDiff{chunksize}
end

function ConstructionBase.constructorof(::Type{<:AutoPolyesterForwardDiff{chunksize}}) where {chunksize}
    return AutoPolyesterForwardDiff{chunksize}
end

end

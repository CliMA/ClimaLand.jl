module ADTypesEnzymeCoreExt

using ADTypes: ADTypes, AutoEnzyme
using EnzymeCore: EnzymeCore

ADTypes.mode(::AutoEnzyme{M}) where {M <: EnzymeCore.ForwardMode} = ADTypes.ForwardMode()
ADTypes.mode(::AutoEnzyme{M}) where {M <: EnzymeCore.ReverseMode} = ADTypes.ReverseMode()

end

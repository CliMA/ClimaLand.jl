module SciMLOperatorsStaticArraysCoreExt

import SciMLOperators
import StaticArraysCore

function Base.copyto!(L::SciMLOperators.MatrixOperator,
        rhs::Base.Broadcast.Broadcasted{<:StaticArraysCore.StaticArrayStyle})
    (copyto!(L.A, rhs); L)
end

end #module

module SciMLOperatorsSparseArraysExt

import SciMLOperators
import SparseArrays: sparse, issparse
import SparseArrays
import LinearAlgebra

SparseArrays.sparse(L::SciMLOperators.MatrixOperator) = sparse(L.A)
SparseArrays.issparse(L::SciMLOperators.MatrixOperator) = issparse(L.A)
SparseArrays.sparse(L::SciMLOperators.ScaledOperator) = L.λ * sparse(L.L)
SparseArrays.sparse(L::SciMLOperators.AddedOperator) = sum(sparse, L.ops)
SparseArrays.sparse(L::SciMLOperators.ComposedOperator) = prod(sparse, L.ops)
function SparseArrays.sparse(L::SciMLOperators.TensorProductOperator)
    LinearAlgebra.kron(sparse.(AbstractMatrix, L.ops)...)
end

end

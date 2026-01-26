"""
    SciMLBase.allows_arbitrary_number_types(alg::T)
        where {T <: ClimaTimeSteppers.RosenbrockAlgorithm}

Return `true`. Enable RosenbrockAlgorithms to run with `ClimaUtilities.ITime`.
"""
function SciMLBase.allows_arbitrary_number_types(alg::T) where {T <: ClimaTimeSteppers.RosenbrockAlgorithm}
    true
end

"""
    SciMLBase.allows_arbitrary_number_types(alg::T)
        where {T <: ClimaTimeSteppers.IMEXAlgorithm}

Return `true`. Enable IMEXAlgorithms to run with `ClimaUtilities.ITime`.
"""
function SciMLBase.allows_arbitrary_number_types(alg::T) where {T <: ClimaTimeSteppers.IMEXAlgorithm}
    true
end

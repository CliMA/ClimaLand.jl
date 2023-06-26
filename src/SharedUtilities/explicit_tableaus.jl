### This is only here temporarily -
### soon it will be moved to ClimaTimeSteppers

import ClimaTimeSteppers: ERKAlgorithmName, ExplicitTableau
using StaticArrays

"""
    RK4

The RK4 algorithm, a Runge-Kutta method with 4 stages and
4th order accuracy.
"""
struct RK4 <: ERKAlgorithmName end
function ExplicitTableau(::RK4)
    return ExplicitTableau(;
        a = @SArray([
            0 0 0 0
            1/2 0 0 0
            0 1/2 0 0
            0 0 1 0
        ]),
        b = @SArray([1 / 6, 1 / 3, 1 / 3, 1 / 6])
    )
end
